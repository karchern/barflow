#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { barcode_counter } from './modules/count_barcodes'
include { barseq_qc_wf } from './modules/barseq_qc'
include {
    createSampleInputChannelAndDecideIfToRun2Fast2Q
    stop_after_barcode_extraction_and_warn
} from './modules/utils'
include { fitness_analysis } from './modules/fitness_analysis'

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

    sample_goodbarcodes_library_map = Channel.fromPath(params.sample_goodbarcodes_library_map)
                                      .splitCsv(header:false)
                                      .map { row ->
                                          def mapPath = row[1] as String
                                          def resolvedMapPath = mapPath.startsWith('/') ? mapPath : "${projectDir}/${mapPath}"
                                          tuple(row[0], file(resolvedMapPath), row[2])
                                      } //TODO: This can be moved into the conditional statement - but I have a pipeline running so don't want to mess things up

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
        all_counts_with_good_barcodes_and_contigs_and_position_ch = barcode_counter.out.result.join(sample_goodbarcodes_library_map)
    }
    else {
        all_counts_with_good_barcodes_and_contigs_and_position_ch = reads_ch.join(sample_goodbarcodes_library_map)
    }

    // use utility helper to handle the "stop after barcode extraction" behavior
    if( stop_after_barcode_extraction_and_warn(all_counts_with_good_barcodes_and_contigs_and_position_ch) ) {
        // just return from workflow body: pipeline ends successfully
        return
    }

    barseq_qc_wf(
        all_counts_with_good_barcodes_and_contigs_and_position_ch
    )
        
    // only run fitness analysis on samples that passed QC.
    // adapt comparisons accordingly (see above)
    fitness_analysis(
        barseq_qc_wf.out.fitness_analysis_input,
        filterConfig,
        mbarqConfig,
    )

    ch_pre_mbarq = fitness_analysis.out.pre_mbarq_qc_data_aggregated_ch
    ch_post_mbarq = fitness_analysis.out.post_mbarq_qc_data_aggregated_ch

    res = ch_pre_mbarq.join(ch_post_mbarq)
               .map { name, pre_m, post_m -> tuple(name, pre_m, post_m) }

    create_master_comparison_log_proc(
        res
    )

    concat_mbarq_qc_results(
        create_master_comparison_log_proc.out.master_log_ch.collect()
    )

    get_final_master_log(
            barseq_qc_wf.out.comparison_validation_log_pre_qc,
            barseq_qc_wf.out.comparison_validation_log_post_qc,
            concat_mbarq_qc_results.out.all_comparisons_master_log
    )

}

process create_master_comparison_log_proc {
    tag { comparison_name }
    label 'r_basic'
    //publishDir "${params.outdir}/comparison_master_logs", mode: 'copy', overwrite: true

    input:
    tuple val(comparison_name), path(pre_mbarq_log), path(post_mbarq_log)

    output:
    path("${comparison_name}.master_log.tsv"), emit: master_log_ch

    script:
    """
    join_master_logs.R "${comparison_name}" "${pre_mbarq_log}" "${post_mbarq_log}" "${comparison_name}.master_log.tsv"
    """
}

process concat_mbarq_qc_results {
    label 'r_basic'
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true

    input:
    path master_files

    output:
    path("all_comparisons.master_log.tsv"), emit: all_comparisons_master_log

    script:
    """
    # Get column names from the first file (assuming all files have the same columns)
    header=\$(head -n 1 ${master_files[0]})
    echo "\$header" > all_comparisons.master_log.tsv
    # Concatenate the rest files, skipping the header of each subsequent file
    for f in ${master_files.join(' ')}; do cat \${f} | tail -n 1 >> all_comparisons.master_log.tsv; done
    """
}

process get_final_master_log {
    label 'r_basic'
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true

    input:
    path pre_qc_logs
    path post_qc_logs
    path barseq_qc_logs

    output:
    path("final_master_log.tsv")

    script:
    """
    join_master_logs_2.R "${pre_qc_logs}" "${post_qc_logs}" "${barseq_qc_logs}" "final_master_log.tsv"
    """
}
