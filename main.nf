#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.json.JsonSlurper

include { create_counts } from './modules/create_barcode_count_matrix'
include {
    buildSampleIndex
    resolveSelectors
    buildComparisonList
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




    // decide input to create_counts
    if ( params.samplesheet && !params.twofast2q_folder ) {

        // current behaviour

        Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:false)
        .map { row ->
            def filename = file(row[0]).name.replaceFirst(/(?i)\.(fastq|fq)(\.gz)?|_sequence\.txt\.gz$/, '')
            tuple(filename, row[0])
        }
        .set { reads_ch_raw }

        reads_ch = reads_ch_raw
        create_counts(reads_ch, good_barcodes_ch)

    }
    else if ( !params.samplesheet && params.twofast2q_folder ) {

        Channel
            .fromPath("${params.twofast2q_folder}/*2fast2q")
            .map { Path p ->
                // remove trailing ".2fast2q" (case-insensitive)
                def sample_id = p.name.replaceFirst(/(?i)\.2fast2q$/, '')
                tuple(sample_id, p)
            }
            .toList()
            .set { all_counts_list_ch }
    }
    else {
        log.warn "\u001b[33mExactly one of --samplesheet or --twofast2q_folder must be set (not both, not neither).\u001b[0m"
        System.exit(1)
    }

    // if we ran create_counts, collect its output into all_counts_list_ch
    if ( params.samplesheet && !params.twofast2q_folder ) {
        create_counts.out.result
            .toList()
            .set { all_counts_list_ch }
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

