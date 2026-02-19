nextflow.enable.dsl=2

process merge_barcode_matrices_process {

    tag { comparison_name }
    publishDir "${params.outdir}/comparisons", mode: 'copy'

    input:
    tuple val(comparison_name),
          val(treat_ids),   path(treat_paths),
          val(ctrl_ids),    path(ctrl_paths)

    output:
    path "${comparison_name}.merged.barcode.matrices.tsv", emit: merged_matrices_ch
    path "treatments.tsv", emit: treatments_ch
    path "controls.tsv", emit: controls_ch

    script:
    """

    paste <(printf "%s\n" ${treat_ids.join(' ')}) <(printf "%s\n" ${treat_paths.join(' ')}) > treatments.tsv
    paste <(printf "%s\n" ${ctrl_ids.join(' ')} ) <(printf "%s\n" ${ctrl_paths.join(' ')} ) > controls.tsv
    touch ${comparison_name}.merged.barcode.matrices.tsv

    join_barcode_matrices.r treatments.tsv controls.tsv ${comparison_name}.merged.barcode.matrices.tsv
    """
}

workflow merge_and_analyze {

  take:
  comparisons_ch

  main:
  merged_result = merge_barcode_matrices_process(comparisons_ch)
  mbarq_results = mbarq_analysis(merged_result.merged_matrices_ch, merged_result.treatments_ch, merged_result.controls_ch)

  emit:
  mbarq_results
}
