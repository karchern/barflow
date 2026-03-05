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
    tuple val(comparison_name),
          path("${comparison_name}.master.log.tx"),
          emit: log_ch

    script:
    """
    paste <(printf "%s\n" ${treat_ids.join(' ')}) <(printf "%s\n" ${treat_paths.join(' ')}) > treatments.tsv
    paste <(printf "%s\n" ${ctrl_ids.join(' ')} ) <(printf "%s\n" ${ctrl_paths.join(' ')} ) > controls.tsv
    touch ${comparison_name}.merged.barcode.matrices.tsv

    echo "### merge_barcode_matrices_process for ${comparison_name} ###" > ${comparison_name}.master.log.tx

    join_barcode_matrices.r treatments.tsv controls.tsv ${good_barcodes_csv_path} ${comparison_name}.merged.barcode.matrices.csv

    # create initial log
    {
        echo "comparison_name=${comparison_name}"
        echo "treatments.tsv:\n\$(cat treatments.tsv)"
        echo "controls.tsv:\n\$(cat controls.tsv)"
    } >> ${comparison_name}.master.log.tx
    """
}


process create_metadata_file_for_mbarq_process {

    tag { comparison_name }
    publishDir "${params.outdir}/comparisons", mode: 'copy'

    input:
    tuple val(comparison_name),
          val(treat_ids),
          val(ctrl_ids),
          path(log_path)

    output:
    // pure data for downstream
    tuple val(comparison_name),
          path("${comparison_name}.mbarq.meta.csv"),
          emit: mbarq_meta_ch

    // log stream (same log_path, updated)
    tuple val(comparison_name),
          path(log_path),
          emit: log_ch

    script:
    """

    echo "### create_metadata_file_for_mbarq_process for ${comparison_name} ###" >> ${log_path}

    create_mbarq_meta_matrix.py \
      --treat-ids ${treat_ids.join(' ')} \
      --ctrl-ids  ${ctrl_ids.join(' ')} \
      --batch batch1 \
      --output ${comparison_name}.mbarq.meta.csv

    # append note to log (optional)
    echo "created meta: ${comparison_name}.mbarq.meta.csv" >> ${log_path}
    """
}


process filter_barcodes_in_merged_matrices {

    tag { comparison_name }
    label 'r_basic'
    publishDir "${params.outdir}/filtered_barcode_matrices", mode: 'copy'

    input:
    tuple val(comparison_name),
          path(merged_matrices_path),
          path(mbarq_meta_path),
          path(log_path)
    val(filter_config)

    output:
    // filtered matrix only
    tuple val(comparison_name),
          path("${comparison_name}.filtered.merged.barcode.matrices.csv"),
          emit: filtered_matrices_ch

    // log stream (updated)
    tuple val(comparison_name),
          path(log_path),
          emit: log_ch

    script:
    """

    echo "### filter_barcodes_in_merged_matrices for ${comparison_name} ###" >> ${log_path}

    echo "lowly_abundant_barcode_cutoff=${filter_config.lowly_abundant_barcode_cutoff}" >> ${log_path}
    echo "filter_on_what=${filter_config.filter_on_what}" >> ${log_path}
    echo "remove_all_0_barcodes=${filter_config.remove_all_0_barcodes}" >> ${log_path}

    filter_merged_barcode_matrix.r \
        ${merged_matrices_path} \
        ${mbarq_meta_path} \
        ${filter_config.filter_on_what} \
        ${filter_config.lowly_abundant_barcode_cutoff} \
        ${filter_config.remove_all_0_barcodes} \
        ${comparison_name}.filtered.merged.barcode.matrices.csv

    echo "filtered matrix: ${comparison_name}.filtered.merged.barcode.matrices.csv" >> ${log_path}
    """
}


process run_mbarq_process {

    tag { comparison_name }
    publishDir "${params.outdir}/mbarq/${comparison_name}", mode: 'copy'

    input:
    tuple val(comparison_name),
          path(merged_matrices_path),
          path(mbarq_meta_path),
          path(log_path)
    val(mbarq_config)

    output:
    tuple val(comparison_name),
          path("${comparison_name}.filtered.merged.barcode.matrices_treated_vs_control.log"),
          emit: mbarq_treated_vs_control_log

    tuple val(comparison_name),
          path("${comparison_name}.filtered.merged.barcode.matrices_rra_results.csv"),
          emit: mbarq_rra_results

    tuple val(comparison_name),
          path("${comparison_name}.filtered.merged.barcode.matrices_barcodes_results.csv"),
          emit: mbarq_barcodes_results

    tuple val(comparison_name),
          path("${comparison_name}.filtered.merged.barcode.matrices_Experiment.log"),
          emit: mbarq_experiment_log

    tuple val(comparison_name),
          path(log_path),
          emit: mbarq_log          

    script:
    """

    echo "Running mbarq for comparison ${comparison_name} with normalization method ${mbarq_config.normalization}" >> ${log_path}

    mbarq analyze -i ${merged_matrices_path} \
      -s ${mbarq_meta_path} \
      --treatment_column treatment --baseline control --norm_method ${mbarq_config.normalization}

    echo "mbarq finished for comparison ${comparison_name}" >> ${log_path}

    # ensure log is in the output folder
    #cp ${log_path} ${comparison_name}.master.log.tx
    """
}

process plot_valcano_plots {

    tag { comparison_name }
    publishDir "${params.outdir}/mbarq/${comparison_name}", mode: 'copy'

    input:
    tuple val(comparison_name),
          path(mbarq_results_path),

    output:
    tuple val(comparison_name),
            path("${comparison_name}.volcano_plot.pdf"),
            emit: volcano_plot_data_ch  

    script:
    """

    echo "Running mbarq for comparison ${comparison_name} with normalization method ${mbarq_config.normalization}" >> ${log_path}

    mbarq analyze -i ${merged_matrices_path} \
      -s ${mbarq_meta_path} \
      --treatment_column treatment --baseline control --norm_method ${mbarq_config.normalization}

    echo "mbarq finished for comparison ${comparison_name}" >> ${log_path}

    # ensure log is in the output folder
    #cp ${log_path} ${comparison_name}.master.log.tx
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

    // treat + ctrl + initial log
    ch_joined = mr.treat_ids_ch
                  .join(mr.ctrl_ids_ch)
                  .join(mr.log_ch)
    // => [comparison_name, treat_ids, ctrl_ids, log_path]

    meta_step     = create_metadata_file_for_mbarq_process(ch_joined)
    mbarq_meta_ch = meta_step.mbarq_meta_ch     // (comparison_name, meta_path)
    log_after_meta = meta_step.log_ch          // (comparison_name, log_path)

    // merged matrices + meta + log
    merged_with_meta = mr.merged_matrices_ch
                         .join(mbarq_meta_ch)
                         .join(log_after_meta)
    // => [comparison_name, merged_matrices_path, meta_path, log_path]

    filtered_step    = filter_barcodes_in_merged_matrices(
        merged_with_meta,
        filter_config
    )
    filtered_mats_ch = filtered_step.filtered_matrices_ch   // (comparison_name, filtered_matrix)
    log_after_filter = filtered_step.log_ch                 // (comparison_name, log_path)

    // join filtered matrix + meta, then add log
    to_mbarq_data = filtered_mats_ch.join(mbarq_meta_ch)
    // => [comparison_name, filtered_matrix_path, meta_path]

    to_mbarq = to_mbarq_data.join(log_after_filter)
    // => [comparison_name, filtered_matrix_path, meta_path, log_path]

    mbarq_results = run_mbarq_process(
        to_mbarq,
        mbarq_config
    )

  emit:
    mbarq_meta_ch
}
