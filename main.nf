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

params.singularity  = params.singularity == true
// raw CLI param (may or may not exist)
def cli_twofast2q = params.containsKey('2fast2q_folder') ? params['2fast2q_folder'] : null

// unified internal params, with explicit defaults
params.samplesheet      = params.containsKey('samplesheet')      ? params.samplesheet      : null
params.twofast2q_folder = params.containsKey('twofast2q_folder') ? params.twofast2q_folder : cli_twofast2q
params.mbarq_normalization = params.containsKey('mbarq_normalization') ? params.mbarq_normalization : 'median'
params.remove_all_0_barcodes = params.containsKey('remove_all_0_barcodes') ? params.remove_all_0_barcodes : false
params.lowly_abundant_barcode_cutoff = params.containsKey('lowly_abundant_barcode_cutoff') ? params.lowly_abundant_barcode_cutoff : 0
params.filter_on_what = params.containsKey('filter_on_what') ? params.filter_on_what : 'treatments_and_controls' // allowed values: 'treatments_and_controls', 'treatments_only', 'controls_only'

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
    'treatments_and_controls',
    'treatments_only',
    'controls_only'
] as Set

if( !allowedFilterValues.contains(params.filter_on_what as String) ) {
    log.error "Invalid value for --filter_on_what: '${params.filter_on_what}'. Allowed values: ${allowedFilterValues.join(', ')}"
    System.exit(1)
}

workflow {

    // load comparisons json
    comparisons = load_json(params.comparisons)



    // good barcodes
    good_barcodes_ch = Channel.fromPath(params.good_barcodes_csv)
                              .map { it }
                              .first()




    input_info = createSampleInputChannelAndDecideIfToRun2Fast2Q(
        params.samplesheet,
        params.twofast2q_folder,
        good_barcodes_ch
    )

    def reads_ch   = input_info.reads_ch
    def run_counts = input_info.run_counts

    if( run_counts ) {
        create_counts(reads_ch, good_barcodes_ch)
        all_counts_list_ch = create_counts.out.result.toList()
    }
    else {
        all_counts_list_ch = reads_ch.toList()
    }

    // downstream as before
    all_counts_list_ch
        .map { tuples ->
            buildComparisonList(tuples, comparisons)
        }
        .set { comparisons_ch }

    comparisons_ch
        .flatMap { it }
        .map { name, treat_list, ctrl_list ->
            def treat_ids   = treat_list.collect { sid, p -> sid }
            def treat_paths = treat_list.collect { sid, p -> p }
            def ctrl_ids    = ctrl_list.collect { sid, p -> sid }
            def ctrl_paths  = ctrl_list.collect { sid, p -> p }
            tuple(name, treat_ids, treat_paths, ctrl_ids, ctrl_paths)
        }
        .set { merge_inputs_ch }

    merge_and_analyze(
        merge_inputs_ch,
        good_barcodes_ch,
        params.mbarq_normalization,
        params.lowly_abundant_barcode_cutoff,
        params.filter_on_what,
        params.remove_all_0_barcodes
    )
}




def load_json(
    String path
) {

    def jsonFile = file(path)
    def jsonData = new JsonSlurper().parseText(jsonFile.text)

    return jsonData
}

