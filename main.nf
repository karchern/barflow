#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.json.JsonSlurper

include { create_counts } from './modules/create_barcode_count_matrix'
include {
    buildSampleIndex
    resolveSelectors
    buildComparisonList
    createSampleInputChannelAndDecideIfToRun2Fast2Q
} from './modules/utils.nf'
include { merge_and_analyze } from './modules/merge_and_mbarq'

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

workflow {

    // load comparisons json
    comparisons = load_json(params.comparisons)

    sample_goodbarcodesfile_map = Channel.fromPath(params.sample_goodbarcodesfile_map)
                                      .splitCsv(header:false)
                                      .map { row -> tuple(row[0], row[1]) } //TODO: This can be moved into the conditional statement - but I have a pipeline running so don't want to mess things up

    input_info = createSampleInputChannelAndDecideIfToRun2Fast2Q(
        params.samplesheet,
        params.twofast2q_folder
    )

    def reads_ch   = input_info.reads_ch
    def run_counts = input_info.run_counts



    if( run_counts ) {
        // join read_ch and sample_goodbarcodesfile_map on first field to get 3-element tuples
        reads_ch = reads_ch
            .join(sample_goodbarcodesfile_map)
        create_counts(reads_ch)
        all_counts_list_ch = create_counts.out.result.toList()
    }
    else {
        all_counts_list_ch = reads_ch.toList()
    }

    if( params.stop_after_barcode_extraction as boolean ) {
        // force materialization so that create_counts (or its cache) is actually used
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

    // downstream as before
    all_counts_list_ch
        .map { tuples ->
            buildComparisonList(tuples, comparisons)
        }
        .set { comparisons_ch }

    comparisons_ch
        .flatMap { it }
        .map { name, treat_list, ctrl_list, good_barcodes_file ->
            def treat_ids   = treat_list.collect { sid, p -> sid }
            def treat_paths = treat_list.collect { sid, p -> p }
            def ctrl_ids    = ctrl_list.collect { sid, p -> sid }
            def ctrl_paths  = ctrl_list.collect { sid, p -> p }
            tuple(name, treat_ids, treat_paths, ctrl_ids, ctrl_paths, good_barcodes_file)
        }
        .set { merge_inputs_ch }

    merge_and_analyze(
        merge_inputs_ch,
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

