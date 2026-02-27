nextflow.enable.dsl=2

process merge_barcode_matrices_process {

    tag { comparison_name }
    publishDir "${params.outdir}/comparisons", mode: 'copy'

    input:
    tuple val(comparison_name),
          val(treat_ids),   path(treat_paths),
          val(ctrl_ids),    path(ctrl_paths) 
    val(good_barcodes_csv_path)

    output:
    tuple val(comparison_name), val(treat_ids), emit: treat_ids_ch
    tuple val(comparison_name), val(ctrl_ids), emit: ctrl_ids_ch
    tuple val(comparison_name), path("${comparison_name}.merged.barcode.matrices.csv"), emit: merged_matrices_ch
    tuple val(comparison_name), path("treatments.tsv"), emit: treatments_ch
    tuple val(comparison_name), path("controls.tsv"), emit: controls_ch

    script:
    """

    paste <(printf "%s\n" ${treat_ids.join(' ')}) <(printf "%s\n" ${treat_paths.join(' ')}) > treatments.tsv
    paste <(printf "%s\n" ${ctrl_ids.join(' ')} ) <(printf "%s\n" ${ctrl_paths.join(' ')} ) > controls.tsv
    touch ${comparison_name}.merged.barcode.matrices.tsv

    join_barcode_matrices.r treatments.tsv controls.tsv ${good_barcodes_csv_path} ${comparison_name}.merged.barcode.matrices.csv
    """
}

process create_metadata_file_for_mbarq_process {

    tag { comparison_name }
    publishDir "${params.outdir}/comparisons", mode: 'copy'

    input:
    tuple val(comparison_name), 
    val(treat_ids),
    val(ctrl_ids)

    output:
    tuple val(comparison_name),
        path("${comparison_name}.mbarq.meta.csv"), 
        emit: mbarq_results
    

    script:
    """
    create_mbarq_meta_matrix.py \
    --treat-ids ${treat_ids.join(' ')} \
    --ctrl-ids  ${ctrl_ids.join(' ')} \
    --batch batch1 \
    --output ${comparison_name}.mbarq.meta.csv    
    """
}

process run_mbarq_process {
    
    tag { comparison_name }
    publishDir "${params.outdir}/mbarq", mode: 'copy'

    input:
    tuple val(comparison_name),
          path(merged_matrices_path),
          path(mbarq_meta_path)
    

    output:
    tuple val(comparison_name),
        path("*"),
        emit: mbarq_results
    

    script:
    """
    mbarq analyze -i ${merged_matrices_path} \
    -s ${mbarq_meta_path} \
    --treatment_column treatment --baseline control
    """
}

// process add_locus_tags_to_barcode_matrices_process {

//     tag { comparison_name }
//     publishDir "${params.outdir}/comparisons", mode: 'copy'

//     input:
//     tuple val(comparison_name), path(merged_matrices_path)
//     path good_barcodes_csv_path

//     output:
//     tuple val(comparison_name), path("${comparison_name}.merged.barcode.matrices.with_locustag.tsv"), emit: merged_matrices_ch

//     script:
//     """
//     add_locus_tag_information_to_merged_matrices.r \
//         --barcode-matrix ${merged_matrices_path} \
//         --good-barcodes ${good_barcodes_csv_path} \
//         --output ${comparison_name}.merged.barcode.matrices.with_locustag.tsv
//     """
// }

workflow merge_and_analyze {

  take:
  comparisons_ch
  good_barcodes_ch

  main:
  mr = merge_barcode_matrices_process(comparisons_ch, good_barcodes_ch)
  ch_joined = mr.treat_ids_ch.join(mr.ctrl_ids_ch) 
  mbarq_meta = create_metadata_file_for_mbarq_process(ch_joined)

  mbarq_results = run_mbarq_process(
    mr.merged_matrices_ch.join(mbarq_meta)
  )

  emit:
  mbarq_meta
}
