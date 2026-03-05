nextflow.enable.dsl=2

process merge_barcode_matrices_process {

    tag { comparison_name }
    label 'r_basic'
    publishDir "${params.outdir}/merged_barcode_matrices", mode: 'copy'

    input:
    tuple val(comparison_name),
          val(treat_ids),   path(treat_paths),
          val(ctrl_ids),    path(ctrl_paths) 
    path(good_barcodes_csv_path)

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
    publishDir "${params.outdir}/mbarq/${comparison_name}", mode: 'copy'

    input:
    tuple val(comparison_name),
          path(merged_matrices_path),
          path(mbarq_meta_path)
    val(mbarq_config)
    

    output:
    tuple val(comparison_name),
        path("*"),
        emit: mbarq_results
    

    script:
    """
    mbarq analyze -i ${merged_matrices_path} \
    -s ${mbarq_meta_path} \
    --treatment_column treatment --baseline control --norm_method ${mbarq_config.mbarq_normalization}
    """
}

process filter_barcodes_in_merged_matrices {

    tag { comparison_name }
    label 'r_basic'
    publishDir "${params.outdir}/filtered_barcode_matrices", mode: 'copy'

    input:
    tuple val(comparison_name),
          path(merged_matrices_path),
          path(mbarq_meta_path)
    val(filter_config)

    output:
    tuple val(comparison_name),
          path("${comparison_name}.filtered.merged.barcode.matrices.csv"),
          emit: merged_matrices_filtered_ch

    script:
    """
    echo "lowly_abundant_barcode_cutoff=${filter_config.lowly_abundant_barcode_cutoff}"
    echo "filter_on_what=${filter_config.filter_on_what}"
    echo "remove_all_0_barcodes=${filter_config.remove_all_0_barcodes}"

    filter_merged_barcode_matrix.r \
        ${merged_matrices_path} \
        ${mbarq_meta_path} \
        ${filter_config.filter_on_what} \
        ${filter_config.lowly_abundant_barcode_cutoff} \
        ${filter_config.remove_all_0_barcodes} \
        ${comparison_name}.filtered.merged.barcode.matrices.csv
    """
}

workflow merge_and_analyze {

  take:
    comparisons_ch
    good_barcodes_ch
    filter_config
    mbarq_config

  main:
    mr = merge_barcode_matrices_process(comparisons_ch, good_barcodes_ch)

    ch_joined  = mr.treat_ids_ch.join(mr.ctrl_ids_ch)
    mbarq_meta = create_metadata_file_for_mbarq_process(ch_joined)

    merged_with_meta = mr.merged_matrices_ch.join(mbarq_meta)

    mr_filtered = filter_barcodes_in_merged_matrices(
        merged_with_meta,
        filter_config
    )

    mbarq_results = run_mbarq_process(
        mr_filtered.merged_matrices_filtered_ch.join(mbarq_meta),
        mbarq_config
    )

  emit:
    mbarq_meta
}
