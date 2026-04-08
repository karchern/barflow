import groovy.json.JsonSlurper

include { get_comparison_status as validate_comparisons_pre_qc } from './../utils'
include { get_comparison_status as validate_comparisons_post_qc } from './../utils'
include {
    buildComparisonList
} from './../utils'

// Workflow helper: load and validate comparisons, then build the derived channels
// Returns a list: [comparisons_ch, comparisons_ch, comparisons_status_ch]
def check_and_validate_comparisons_post_qc(input_channel, comparisons) {
    // build comparisons channel (each element is a list of comparison tuples)
    def comparisons_built_ch = input_channel.
        toList()
        .map { tuples ->
            buildComparisonList(tuples, comparisons)
        }

    // prepare inputs for merging: (name, treat_ids, treat_paths, ctrl_ids, ctrl_paths, good_barcodes_file)
    def fitness_analysis_input_ch = comparisons_built_ch
        .flatMap { it }
        .map { name, treat_list, ctrl_list, good_barcodes_file, status, status_detail ->
            def treat_ids   = treat_list.collect { sid, p -> sid }
            def treat_paths = treat_list.collect { sid, p -> p }
            def ctrl_ids    = ctrl_list.collect { sid, p -> sid }
            def ctrl_paths  = ctrl_list.collect { sid, p -> p }
            tuple(name, treat_ids, treat_paths, ctrl_ids, ctrl_paths, good_barcodes_file)
        }

    // build a list of comparison status tuples and materialize to a single channel (toList)
    def comparisons_status_list = comparisons_built_ch
        .flatMap { it }
        .map { name, treat_list, ctrl_list, good_barcodes_file, status, status_detail ->
            tuple(name, status, status_detail)
        }
        .toList()

    return [fitness_analysis_input_ch, comparisons_status_list]
}

def load_json(String path) {
    def jsonFile = file(path)
    def jsonData = new JsonSlurper().parseText(jsonFile.text)
    return jsonData
}

workflow barseq_qc_wf {

    take:
    all_counts_with_good_barcodes_and_contigs_and_position_ch   // emits: [ sample_id, counts_path, sample_goodbarcodes_contig_map, library]

    main:

    def comparisons = load_json(params.comparisons)

    all_counts_with_good_barcodes_and_contigs_and_position_ch.
         map { sample_id, counts_path, sample_goodbarcodes_contig_map, library -> tuple(sample_id, counts_path) }
         .set { all_counts } // retain this and pass to most downstream functions to ensure everything is working as intended - at some point this could be improved

    def comps = check_and_validate_comparisons_post_qc(all_counts, comparisons)
    def fitness_analysis_input_ch = comps[0]
    def comparisons_status_list = comps[1]

    // Pre-QC comparison validation (materialized list is acceptable here)
    comparison_validation_log_pre_qc = validate_comparisons_pre_qc(comparisons_status_list, "before_barseq_qc")        

    // Run per-sample BarSeq QC using the channel of counts
    // all_counts_with_good_barcodes_and_contigs_and_position_ch also contains the sample_goodbarcodes_contig_map and library info - so retain the info
    barseq_qc(all_counts_with_good_barcodes_and_contigs_and_position_ch, params.minimum_read_sum_for_qc, params.minimum_median_barcode_count)

    // merge barcode_count_sample_metrics
    barseq_qc.out.
         // extract sample_id from [sample_id, counts_path]
        map { sample_id, metrics_path, qc_passed, median_of_medians, ptr_diag_image, ptr_diag_image_top50 -> metrics_path}
        .toList()
        .set { metrics_paths }
    
    // Based on the barseq_qc output, build a channel of sample_ids that passed QC. We will use this to filter the comparisons.
    sample_ids_passed_ch = barseq_qc.out
        .filter { sample_id, metrics_path, qc_passed, median_of_medians, ptr_diag_image, ptr_diag_image_top50 -> qc_passed.text.trim() == '1' }
        .map { sample_id, metrics_path, qc_passed, median_of_medians, ptr_diag_image, ptr_diag_image_top50 -> sample_id }
        .map { sid -> tuple(sid) }

    // Join the original counts CHANNEL with the passed-samples channel to get only passing samples
    all_counts_PASSED_ch = all_counts.join(sample_ids_passed_ch)

    // After barseq qc and corresponding sample filtering, we need to check and filter the comparisons based on which samples passed QC.
    // Pass the CHANNEL into the comparison-check helper (it expects channel operations)
    def comps_post_filter = check_and_validate_comparisons_post_qc(all_counts_PASSED_ch, comparisons)
    def fitness_analysis_input_post_filter_ch = comps_post_filter[0]
    def comparisons_status_post_filter_list = comps_post_filter[1]

    comparison_validation_log_post_qc = validate_comparisons_post_qc(comparisons_status_post_filter_list, "after_barseq_qc")

    // Publish merged BarSeq sample QC metrics
    collate_barseq_qc_results(metrics_paths)

    comparisons_status_post_filter_PASSED_list = comparisons_status_post_filter_list
        .flatMap { it }
        .filter { name, status, detail -> status != 'ALL_SAMPLES_MISSING' } 
        .map { name, status, detail -> name }

    // Attention: Comparisons with completely missing samples have already been filtered out in filter_ch_comparison_post_filter
    // so this join removes comparisons that have all samples missing, but keeps comparisons that have at least some samples passing QC.
    // This is important to consider downstream.
    fitness_analysis_input = fitness_analysis_input_post_filter_ch
        .join(comparisons_status_post_filter_PASSED_list)

    print_summary_table_p(
        all_counts,
        sample_ids_passed_ch,
        comparisons_status_list,
        comparisons_status_post_filter_PASSED_list
    )

    emit:
    fitness_analysis_input
    comparison_validation_log_pre_qc
    comparison_validation_log_post_qc

}



workflow print_summary_table_p {

    take:
    all_counts_ch
    sample_qcs_PASSED_CH
    comparisons_status_ch
    filter_ch_comparison_post_filter

    main:
    def total_samples_ch = all_counts_ch.toList().map { all_samples_list ->
        all_samples_list.size()
    }

    def samples_after_qc_ch = sample_qcs_PASSED_CH.count()

    def comps_before_ch = comparisons_status_ch.map { x ->
        x.size()
    }

    def comps_after_ch = filter_ch_comparison_post_filter.count()

    print_summary_table(
        total_samples_ch,
        samples_after_qc_ch,
        comps_before_ch,
        comps_after_ch
    )
}

process collate_barseq_qc_results {
    
    label 'python_basic'

    publishDir "${params.outdir}/barseq_qc_collated/", mode: 'copy', overwrite: true

    input:
    path metrics_file_paths

    output:
    path("all_samples.barcode_metrics.tsv"), emit: sample_wise_qc_metrics

    script:
    """
    #echo -e "sample_id,total_reads,n_barcodes,max_count,mean_count,median_count," > all_samples.barcode_metrics.tsv
    #Take header from first file, saver than this shit above :)
    header=\$(head -n 1 ${metrics_file_paths[0]})
    echo "\$header" > all_samples.barcode_metrics.tsv
    for f in *.csv; do
        tail -n +2 "\$f" >> all_samples.barcode_metrics.tsv
    done
    """

}

process barseq_qc {

    label 'python_basic_quick_process'

    publishDir "${params.outdir}/barseq_qc/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(counts_path), path(sample_goodbarcodes_contig_position_map), val(library)
    val(minimum_read_sum_for_qc)
    val(minimum_median_barcode_count)

    output:
    tuple val(sample_id), path("${sample_id}.barcode_metrics.csv"), path("${sample_id}.passed_qc.txt"), path("${sample_id}.median_of_means_over_genomes.csv"), path("${sample_id}_median_sliding_window.png"), path("${sample_id}_median_sliding_window_top50pct_sites.png"), emit: metrics

    """
    barseq_qc.py \
        --sample_id ${sample_id} \
        --input_counts ${counts_path} \
        --sample_goodbarcodes_contig_position_map ${sample_goodbarcodes_contig_position_map} \
        --output_barcode_metrics ${sample_id}.barcode_metrics.csv \
        --output_median_of_means_over_genomes ${sample_id}.median_of_means_over_genomes.csv \
        --output_passed ${sample_id}.passed_qc.txt \
        --min_read_sum_for_qc ${minimum_read_sum_for_qc} \
        --min_median_barcode_count ${minimum_median_barcode_count}
    """
}

process print_summary_table {

    cache false
    debug true

    input:
    val total_samples
    val samples_after_qc
    val comps_before
    val comps_after

    exec:
    println """
        ============================================
                    Pre-mbarq numbers
        ============================================
        | ${"Metric".padRight(28)} | ${"Number".padLeft(8)} |
        --------------------------------------------
        | ${"BarSeq samples before QC".padRight(28)} | ${String.valueOf(total_samples).padLeft(8)} |
        | ${"BarSeq samples after QC".padRight(28)} | ${String.valueOf(samples_after_qc).padLeft(8)} |
        | ${"Comparisons before filtering".padRight(28)} | ${String.valueOf(comps_before).padLeft(8)} |
        | ${"Comparisons after filtering".padRight(28)} | ${String.valueOf(comps_after).padLeft(8)} |
        ============================================
        """.stripIndent()
}