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
        good_barcodes_ch
    )
}




def load_json(
    String path
) {

    def jsonFile = file(path)
    def jsonData = new JsonSlurper().parseText(jsonFile.text)

    return jsonData
}

