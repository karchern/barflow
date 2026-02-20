nextflow.enable.dsl=2

process merge_barcode_matrices_process {

    tag { comparison_name }
    publishDir "${params.outdir}/comparisons", mode: 'copy'

    input:
    tuple val(comparison_name),
          val(treat_ids),   path(treat_paths),
          val(ctrl_ids),    path(ctrl_paths)

    output:
    tuple val(comparison_name),
        path("${comparison_name}.merged.barcode.matrices.tsv"),
        path("treatments.tsv"),
        path("controls.tsv"),
        emit: merged_tuple_ch
    

    script:
    """

    paste <(printf "%s\n" ${treat_ids.join(' ')}) <(printf "%s\n" ${treat_paths.join(' ')}) > treatments.tsv
    paste <(printf "%s\n" ${ctrl_ids.join(' ')} ) <(printf "%s\n" ${ctrl_paths.join(' ')} ) > controls.tsv
    touch ${comparison_name}.merged.barcode.matrices.tsv

    join_barcode_matrices.r treatments.tsv controls.tsv ${comparison_name}.merged.barcode.matrices.tsv
    """
}

process mbarq_analysis_process {

    tag { comparison_name }
    publishDir "${params.outdir}/comparisons", mode: 'copy'

    input:
    tuple val(comparison_name),
          path(merged_matrices_path),
          path(treatments_path),
          path(controls_path)

    output:
    path "${comparison_name}.mbarq.results.tsv", emit: mbarq_results

    script:
    """
    #mbarq.r ${merged_matrices_path} ${treatments_path} ${controls_path} > ${comparison_name}.mbarq.results.tsv
    echo "hello :)" > ${comparison_name}.mbarq.results.tsv
    """
}

workflow merge_and_analyze {

  take:
  comparisons_ch

  main:
  merged_result = merge_barcode_matrices_process(comparisons_ch)
  mbarq_results = mbarq_analysis_process(merged_result.merged_tuple_ch)

  emit:
  mbarq_results
}
