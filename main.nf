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

params.outdir       = params.outdir ?: './results'
params.singularity  = params.singularity == true
params.samplesheet = params.samplesheet ?: './samplesheet.csv'

workflow {

    // load comparisons json
    comparisons = load_json(params.comparisons)

    Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:false)
    .map { row ->
        def filename = file(row[0]).name.replaceFirst(/(?i)\.(fastq|fq)(\.gz)?|_sequence\.txt\.gz$/, '')
        tuple(filename, row[0])
    }
    .set { reads_ch }

    //reads_ch.view { "A: All input FASTQ files: ${it}" }

    good_barcodes_ch = Channel.value(file(params.good_barcodes_csv))

    // 1) Run 2FAST2Q
    create_counts(reads_ch, good_barcodes_ch)

    //create_counts.out.result.view { "Counts output from create_counts: ${it}" }

    // collect into a sample index: [ sample_id : Path ]
    create_counts.out.result
        .toList()
        .set { all_counts_list_ch }

    //all_counts_list_ch.view { "All counts list: ${it}" }

    all_counts_list_ch
        .map { tuples ->
            buildComparisonList(tuples, comparisons)
        }
        .set { comparisons_ch } 

    comparisons_ch
        .flatMap { it }      // turn List<cmp> into stream of cmp
        .map { name, treat_list, ctrl_list ->
            def treat_ids   = treat_list.collect { sid, p -> sid }
            def treat_paths = treat_list.collect { sid, p -> p }
            def ctrl_ids    = ctrl_list.collect { sid, p -> sid }
            def ctrl_paths  = ctrl_list.collect { sid, p -> p }
            tuple(name, treat_ids, treat_paths, ctrl_ids, ctrl_paths)
        }        
        .set { merge_inputs_ch }

    merge_inputs_ch.view()

    // comparisons_ch
    //     .map { cmp ->
    //         def name = cmp.name
    //         def treatPaths   = cmp.treatments.collect { sid, p -> p }
    //         def controlPaths = cmp.controls.collect { sid, p -> p }
    //         tuple(name, treatPaths, controlPaths)
    //     }
    //     .set { merge_inputs_ch }

    // comparisons_ch
    //     .flatMap { list_of_tuples -> list_of_tuples }  // flatten list into stream
    //     .set { comparisons_ch2 }

    // comparisons_ch
    // .map { cmp ->
    //     def name = cmp[0]
    //     def treat = cmp[1]
    //     def control = cmp[2]
    //     tuple(name, treat, control)
    // }
    // .set { merge_inputs_ch }

    // merge_inputs_ch.view()

    merge_and_analyze(
        merge_inputs_ch
    )
    
}



def load_json(
    String path
) {

    def jsonFile = file(path)
    def jsonData = new JsonSlurper().parseText(jsonFile.text)

    return jsonData
}

