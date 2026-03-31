import groovy.json.JsonSlurper

include { get_comparison_status as validate_comparisons } from './../utils'
include { get_comparison_status as filter_comparisons } from './../utils'
include {
    buildComparisonList
} from './../utils'

// Workflow helper: load and validate comparisons, then build the derived channels
// Returns a list: [comparisons_ch, comparisons_ch, comparisons_status_ch]
def check_and_filter_comparisons(input_list, comparisons) {
    // build comparisons channel (each element is a list of comparison tuples)
    def comparisons_built_ch = input_list
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
    all_counts

    main:

    def comparisons = load_json(params.comparisons)
    all_counts_list = all_counts.toList() // materialize to list for multiple passes

    def comps = check_and_filter_comparisons(all_counts_list, comparisons)
    def fitness_analysis_input_ch = comps[0]
    def comparisons_status_list = comps[1]

    // check_and_filter_comparisons will check which comparisons have all their samples passing QC, and filter out the ones that don't. 
    // It will also build a list of comparison status tuples that we can materialize to a channel and publish as a tsv file.
    // Before barseq qc, all comparisons should be present (no filtering based on barseq QC yet), but we can still check their status and also need build the list of comparison status tuples.
    validate_comparisons(comparisons_status_list, "before_barseq_qc")        

    barseq_qc(all_counts_list.flatMap { it }, params.minimum_read_sum_for_qc)

    // merge barcode_count_sample_metrics
    barseq_qc.out.
         // extract sample_id from [sample_id, counts_path]
        map { sample_id, metrics_path, qc_passed -> metrics_path }
        .toList()
        .set { metrics_paths }
    
    // Based on the barseq_qc output, build a channel of sample_ids that passed QC. We will use this to filter the comparisons.
    sample_ids_passed_ch = barseq_qc.out
        .filter { sample_id, metrics_path, qc_passed -> qc_passed.text.trim() == '1' }
        .map { sample_id, metrics_path, qc_passed -> sample_id }
        .map { sid -> tuple(sid) }

    all_counts_PASSED_ch = all_counts_list
    .flatMap { it }
    .join(sample_ids_passed_ch)

    // After barseq qc and corresponding sample filtering, we need to check and filter the comparisons based on which samples passed QC.
    // If a comparison has some samples removed, we keep the comparison but log the missing samples in the comparison status. 
    // If all samples of a comparison are removed, we filter out the comparison entirely.
    def comps_post_filter = check_and_filter_comparisons(all_counts_PASSED_ch.toList(), comparisons)
    //def _ = comps_post_filter[0]
    def fitness_analysis_input_post_filter_ch = comps_post_filter[0]
    def comparisons_status_post_filter_list = comps_post_filter[1]

    filter_comparisons(comparisons_status_post_filter_list, "after_barseq_qc")

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
        all_counts_list,
        sample_ids_passed_ch,
        comparisons_status_list,
        comparisons_status_post_filter_PASSED_list
    )

    emit:
    fitness_analysis_input
}



workflow print_summary_table_p {

    take:
    all_counts_list_ch
    sample_qcs_PASSED_CH
    comparisons_status_ch
    filter_ch_comparison_post_filter

    main:
    def total_samples_ch = all_counts_list_ch.map { all_samples_list ->
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

    publishDir "${params.outdir}/logs/", mode: 'copy', overwrite: true

    input:
    path metrics_file_paths

    output:
    path("all_samples.barcode_metrics.tsv"), emit: sample_wise_qc_metrics

    script:
    """
    echo -e "sample_id,total_reads,n_barcodes,max_count,mean_count,median_count" > all_samples.barcode_metrics.tsv
    for f in *.csv; do
        tail -n +2 "\$f" >> all_samples.barcode_metrics.tsv
    done
    """

}

process barseq_qc {

    label 'python_basic'

    input:
    tuple val(sample_id), path(counts_path)
    val(minimum_read_sum_for_qc)

    output:
    tuple val(sample_id), path("${sample_id}.barcode_metrics.csv"), path("${sample_id}.passed_qc.txt"), emit: metrics

    """
    barseq_qc.py \
        --sample_id ${sample_id} \
        --input ${counts_path} \
        --output ${sample_id}.barcode_metrics.csv \
        --output_passed ${sample_id}.passed_qc.txt \
        --min_read_sum_for_qc ${minimum_read_sum_for_qc}
    
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