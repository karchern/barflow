#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.json.JsonSlurper

include { barcode_counter } from './modules/count_barcodes'
include {
    buildComparisonList
    createSampleInputChannelAndDecideIfToRun2Fast2Q
    print_summary_table
} from './modules/utils.nf'
include { fitness_analysis } from './modules/fitness_analysis'
include { get_comparison_status as get_comparison_status_before_barseq_qc } from './modules/get_comparison_status'
include { get_comparison_status as get_comparison_status_after_barseq_qc } from './modules/get_comparison_status'

//params.singularity  = params.singularity == true
// raw CLI param (may or may not exist)
def cli_twofast2q = params.containsKey('2fast2q_folder') ? params['2fast2q_folder'] : null
if( !params.twofast2q_folder && cli_twofast2q )
    params.twofast2q_folder = cli_twofast2q

def mbarqConfig = [
    normalization: params.mbarq_normalization,
]

def filterConfig = [
    lowly_abundant_barcode_cutoff : params.lowly_abundant_barcode_cutoff as Integer,
    filter_on_what                : params.filter_on_what as String,
    remove_all_0_barcodes         : params.remove_all_0_barcodes as Boolean
]
// Sanity check for existing 2fast2q folder
if( params.twofast2q_folder ) {

    def dir = file(params.twofast2q_folder)

    if( !dir.exists() ) {
        log.error "Provided --twofast2q_folder '${params.twofast2q_folder}' does not exist."
        System.exit(1)
    }
    if( !dir.isDirectory() ) {
        log.error "Provided --twofast2q_folder '${params.twofast2q_folder}' is not a directory."
        System.exit(1)
    }

    // list files ending in '2fast2q' (top level only)
    def files = dir.listFiles()?.findAll { f ->
        f.isFile() && f.name.endsWith('2fast2q')
    } ?: []

    if( !files ) {
        log.error "Directory '${params.twofast2q_folder}' contains no files ending in '2fast2q'."
        System.exit(1)
    }

    log.info "Found ${files.size()} file(s) ending in '2fast2q' in '${params.twofast2q_folder}\nWill try to build mbarq comparisons from those :)'."
}

// sanity check for filter_on_what
def allowedFilterValues = [
    'both',
    'treatments',
    'controls'
] as Set

if( !allowedFilterValues.contains(params.filter_on_what as String) ) {
    log.error "Invalid value for --filter_on_what: '${params.filter_on_what}'. Allowed values: ${allowedFilterValues.join(', ')}"
    System.exit(1)
}

// Workflow helper: load and validate comparisons, then build the derived channels
// Returns a list: [comparisons_ch, comparisons_ch, comparisons_status_ch]
def check_and_filter_comparisons(all_counts_list_ch, comparisons) {
    // build comparisons channel (each element is a list of comparison tuples)
    def comparisons_ch = all_counts_list_ch
        .map { tuples ->
            buildComparisonList(tuples, comparisons)
        }

    // prepare inputs for merging: (name, treat_ids, treat_paths, ctrl_ids, ctrl_paths, good_barcodes_file)
    def merge_inputs_ch = comparisons_ch
        .flatMap { it }
        .map { name, treat_list, ctrl_list, good_barcodes_file, status, status_detail ->
            def treat_ids   = treat_list.collect { sid, p -> sid }
            def treat_paths = treat_list.collect { sid, p -> p }
            def ctrl_ids    = ctrl_list.collect { sid, p -> sid }
            def ctrl_paths  = ctrl_list.collect { sid, p -> p }
            tuple(name, treat_ids, treat_paths, ctrl_ids, ctrl_paths, good_barcodes_file)
        }

    // build a list of comparison status tuples and materialize to a single channel (toList)
    def comparisons_status_ch = comparisons_ch
        .flatMap { it }
        .map { name, treat_list, ctrl_list, good_barcodes_file, status, status_detail ->
            tuple(name, status, status_detail)
        }
        .toList()

    return [comparisons_ch, merge_inputs_ch, comparisons_status_ch]
}

workflow {

    sample_goodbarcodes_library_map = Channel.fromPath(params.sample_goodbarcodes_library_map)
                                      .splitCsv(header:false)
                                      .map { row -> tuple(row[0], file(row[1]), row[2]) } //TODO: This can be moved into the conditional statement - but I have a pipeline running so don't want to mess things up

    input_info = createSampleInputChannelAndDecideIfToRun2Fast2Q(
        params.samplesheet,
        params.twofast2q_folder
    )

    def reads_ch   = input_info.reads_ch
    def run_counts = input_info.run_counts



    if( run_counts ) {
        // join read_ch and sample_goodbarcodes_library_map on first field to get 3-element tuples
        reads_ch = reads_ch
            .join(sample_goodbarcodes_library_map)
        barcode_counter(reads_ch)
        all_counts_list_ch = barcode_counter.out.result.toList()
    }
    else {
        all_counts_list_ch = reads_ch.toList()
    }

    if( params.stop_after_barcode_extraction as boolean ) {
        // force materialization so that barcode_counter (or its cache) is actually used
        all_counts_list_ch
            .view { tuples ->
                log.warn """
                ================================
                STOP AFTER BARCODE EXTRACTION
                ================================
                Parameter: --stop_after_barcode_extraction true

                Barcode extraction and count matrix generation have completed.
                No downstream merging or mbarq analysis will be performed.
                Disable this behavior by omitting the parameter or setting:
                    --stop_after_barcode_extraction false
                ================================
                """.stripIndent()
            }
    // just return from workflow body: pipeline ends successfully
    return
    }

    comparisons = load_json(params.comparisons)

    def comps = check_and_filter_comparisons(all_counts_list_ch, comparisons)
    def _ = comps[0]
    def merge_inputs_ch = comps[1]
    def comparisons_status_ch = comps[2]

    // check_and_filter_comparisons will check which comparisons have all their samples passing QC, and filter out the ones that don't. 
    // It will also build a list of comparison status tuples that we can materialize to a channel and publish as a tsv file.
    // Before barseq qc, all comparisons should be present (no filtering based on barseq QC yet), but we can still check their status and also need build the list of comparison status tuples.
    get_comparison_status_before_barseq_qc(comparisons_status_ch, "before_barseq_qc")        

    barseq_qc(
        all_counts_list_ch
        .flatMap { it }, 
    params.minimum_read_sum_for_qc
    )

    // merge barcode_count_sample_metrics
    barseq_qc.out.
         // extract sample_id from [sample_id, counts_path]
        map { sample_id, metrics_path, qc_passed -> metrics_path }
        .toList()
        .set { metrics_paths }
    
    // Based on the barseq_qc output, build a channel of sample_ids that passed QC. We will use this to filter the comparisons.
    sample_qcs_PASSED_CH = barseq_qc.out
        .filter { sample_id, metrics_path, qc_passed -> qc_passed.text.trim() == '1' }
        .map { sample_id, metrics_path, qc_passed -> sample_id }
        .map { sid -> tuple(sid) }

    all_counts_list_PASSED_ch = all_counts_list_ch
    .flatMap { it }
    .join(sample_qcs_PASSED_CH)

    // Publish merged BarSeq sample QC metrics
    collate_barseq_qc_results(metrics_paths)

    // After barseq qc and corresponding sample filtering, we need to check and filter the comparisons based on which samples passed QC.
    // If a comparison has some samples removed, we keep the comparison but log the missing samples in the comparison status. 
    // If all samples of a comparison are removed, we filter out the comparison entirely.
    def comps_post_filter = check_and_filter_comparisons(all_counts_list_PASSED_ch.toList(), comparisons)
    def merge_inputs_post_filter_ch = comps_post_filter[1]
    def comparisons_status_post_filter_ch = comps_post_filter[2]

    get_comparison_status_after_barseq_qc(comparisons_status_post_filter_ch, "after_barseq_qc")

    // TODO: Add some more logging statements
    filter_ch_comparison_post_filter = comparisons_status_post_filter_ch
        .flatMap { it }
        .filter { name, status, detail -> status != 'ALL_SAMPLES_MISSING' } 
        .map { name, status, detail -> name }

    // Attention: Comparisons with completely missing samples have already been filtered out in filter_ch_comparison_post_filter
    // so this join removes comparisons that have all samples missing, but keeps comparisons that have at least some samples passing QC.
    // This is important to consider downstream.
    merge_inputs_post_filtered_PASSED_ch = merge_inputs_post_filter_ch
        .join(filter_ch_comparison_post_filter)

    //
    // SUMMARY METRICS FOR print_summary_table
    //

    print_summary_table(
        all_counts_list_ch,
        sample_qcs_PASSED_CH,
        comparisons_status_ch,
        filter_ch_comparison_post_filter
    )
        
    // only run fitness analysis on samples that passed QC.
    // adapt comparisons accordingly (see above)
    fitness_analysis(
        merge_inputs_post_filtered_PASSED_ch,
        filterConfig,
        mbarqConfig,
    )

}

def load_json(
    String path
) {

    def jsonFile = file(path)
    def jsonData = new JsonSlurper().parseText(jsonFile.text)

    return jsonData
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

