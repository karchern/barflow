process pre_mbarq_qc_process {


    tag { comparison_name }
    label 'r_basic'
    publishDir "${params.outdir}/qc_${called_from_where}/${comparison_name}", mode: 'copy'


    input:
    tuple val(comparison_name),
          path(merged_matrices_path),
          path(mbarq_meta_path)
    val(called_from_where)


    output:
    // step log only
    tuple val(comparison_name),
          path("${comparison_name}.mbarq_qc.log.txt"),
          emit: mbarq_qc_log_ch
    
    tuple val(comparison_name),
          path("${comparison_name}.${called_from_where}.mbarq.qc.tsv"),
          emit: mbarq_qc_data_ch



    script:
    """
    do_control_based_qc.r  ${merged_matrices_path} ${mbarq_meta_path} "${comparison_name}.${called_from_where}.mbarq.qc.tsv"
    touch "${comparison_name}.mbarq_qc.log.txt"
    """
}