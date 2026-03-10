nextflow.enable.dsl=2


/*
 * PROCESSES
 */

process merge_barcode_matrices_process {

    tag { comparison_name }
    label 'r_basic'
    publishDir "${params.outdir}/merged_barcode_matrices", mode: 'copy'

    input:
    tuple val(comparison_name),
          val(treat_ids),   path(treat_paths),
          val(ctrl_ids),    path(ctrl_paths), 
          path(good_barcodes_file)

    output:
    tuple val(comparison_name), val(treat_ids), emit: treat_ids_ch
    tuple val(comparison_name), val(ctrl_ids), emit: ctrl_ids_ch
    tuple val(comparison_name), path("${comparison_name}.merged.barcode.matrices.csv"), emit: merged_matrices_ch
    tuple val(comparison_name), path("treatments.tsv"), emit: treatments_ch
    tuple val(comparison_name), path("controls.tsv"), emit: controls_ch

    // step-specific log
    tuple val(comparison_name),
          path("${comparison_name}.merge.log.txt"),
          emit: merge_log_ch

    script:
    """
    paste <(printf "%s\n" ${treat_ids.join(' ')}) <(printf "%s\n" ${treat_paths.join(' ')}) > treatments.tsv
    paste <(printf "%s\n" ${ctrl_ids.join(' ')} ) <(printf "%s\n" ${ctrl_paths.join(' ')} ) > controls.tsv
    touch ${comparison_name}.merged.barcode.matrices.tsv

    join_barcode_matrices.r treatments.tsv controls.tsv ${good_barcodes_file} ${comparison_name}.merged.barcode.matrices.csv

    {
        echo "### merge_barcode_matrices_process for ${comparison_name} ###"
        echo "comparison_name=${comparison_name}"
        echo "treatments.tsv:"
        cat treatments.tsv
        echo "controls.tsv:"
        cat controls.tsv
    } > ${comparison_name}.merge.log.txt
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
    // pure data for downstream
    tuple val(comparison_name),
          path("${comparison_name}.mbarq.meta.csv"),
          emit: mbarq_meta_ch

    // step log
    tuple val(comparison_name),
          path("${comparison_name}.create_meta.log.txt"),
          emit: create_meta_log_ch

    script:
    """
    create_mbarq_meta_matrix.py \
      --treat-ids ${treat_ids.join(' ')} \
      --ctrl-ids  ${ctrl_ids.join(' ')} \
      --batch batch1 \
      --output ${comparison_name}.mbarq.meta.csv

    {
      echo "### create_metadata_file_for_mbarq_process for ${comparison_name} ###"
      echo "created meta: ${comparison_name}.mbarq.meta.csv"
    } > ${comparison_name}.create_meta.log.txt
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
    // filtered matrix only
    tuple val(comparison_name),
          path("${comparison_name}.filtered.merged.barcode.matrices.csv"),
          emit: filtered_matrices_ch

    // step log
    tuple val(comparison_name),
          path("${comparison_name}.filter.log.txt"),
          emit: filter_log_ch

    script:
    """
    {
      echo "### filter_barcodes_in_merged_matrices for ${comparison_name} ###"
      echo "lowly_abundant_barcode_cutoff=${filter_config.lowly_abundant_barcode_cutoff}"
      echo "filter_on_what=${filter_config.filter_on_what}"
      echo "remove_all_0_barcodes=${filter_config.remove_all_0_barcodes}"
    } > ${comparison_name}.filter.log.txt

    filter_merged_barcode_matrix.r \
        ${merged_matrices_path} \
        ${mbarq_meta_path} \
        ${filter_config.filter_on_what} \
        ${filter_config.lowly_abundant_barcode_cutoff} \
        ${filter_config.remove_all_0_barcodes} \
        ${comparison_name}.filtered.merged.barcode.matrices.csv

    echo "filtered matrix: ${comparison_name}.filtered.merged.barcode.matrices.csv" >> ${comparison_name}.filter.log.txt
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

    // step log
    tuple val(comparison_name),
          path("${comparison_name}.mbarq.log.txt"),
          emit: mbarq_log_ch

    script:
    """
    {
      echo "### run_mbarq_process for ${comparison_name} ###"
      echo "norm_method=${mbarq_config.normalization}"
    } > ${comparison_name}.mbarq.log.txt

    mbarq analyze -i ${merged_matrices_path} \
      -s ${mbarq_meta_path} \
      --treatment_column treatment --baseline control \
      --norm_method ${mbarq_config.normalization}

    echo "mbarq finished for comparison ${comparison_name}" >> ${comparison_name}.mbarq.log.txt
    """
}


process plot_volcano_plots {

    tag { comparison_name }
    label 'r_basic'
    publishDir "${params.outdir}/mbarq_analysis/${comparison_name}", mode: 'copy'

    input:
    tuple val(comparison_name),
          path(mbarq_results_path)

    output:
    tuple val(comparison_name),
          path("${comparison_name}.volcano_plot.pdf"),
          emit: volcano_plot_data_ch  

    script:
    """
    plot_volcano.r ${mbarq_results_path} ${comparison_name} ${comparison_name}.volcano_plot.pdf
    """
}


process merge_logs {

    tag { comparison_name }
    // One folder per comparison under mbarq
    publishDir "${params.outdir}/mbarq/${comparison_name}", mode: 'copy'

    input:
    tuple val(comparison_name),
          path(step_logs)

    output:
    tuple val(comparison_name),
          path("${comparison_name}.master.log.txt"),
          emit: master_log_ch

    script:
    """
    cat ${step_logs.join(' ')} > ${comparison_name}.master.log.txt
    """
}


/*
 * WORKFLOW
 */

workflow merge_and_analyze {

  take:
    comparisons_ch
    filter_config
    mbarq_config

  main:
    mr = merge_barcode_matrices_process(comparisons_ch)

    // treat + ctrl for metadata
    ch_joined = mr.treat_ids_ch
                  .join(mr.ctrl_ids_ch)
    // => [comparison_name, treat_ids, ctrl_ids]

    meta_step      = create_metadata_file_for_mbarq_process(ch_joined)
    mbarq_meta_ch  = meta_step.mbarq_meta_ch         // (comparison_name, meta_path)
    create_meta_log_ch = meta_step.create_meta_log_ch

    // merged matrices + meta
    merged_with_meta = mr.merged_matrices_ch
                         .join(mbarq_meta_ch)
    // => [comparison_name, merged_matrices_path, meta_path]

    filtered_step    = filter_barcodes_in_merged_matrices(
        merged_with_meta,
        filter_config
    )
    filtered_mats_ch = filtered_step.filtered_matrices_ch   // (comparison_name, filtered_matrix)
    filter_log_ch    = filtered_step.filter_log_ch

    // join filtered matrix + meta for mbarq
    to_mbarq_data = filtered_mats_ch.join(mbarq_meta_ch)
    // => [comparison_name, filtered_matrix_path, meta_path]

    mbarq_results = run_mbarq_process(
        to_mbarq_data,
        mbarq_config
    )

    mbarq_rra_results      = mbarq_results.mbarq_rra_results
    mbarq_log_ch           = mbarq_results.mbarq_log_ch

    to_volcano = mbarq_rra_results
    volcano_plots = plot_volcano_plots(
        to_volcano
    )

    // collect step logs and build master logs
    ordered_logs = mr.merge_log_ch
                    .join(create_meta_log_ch)   // [id, merge_log, create_meta_log]
                    .join(filter_log_ch)        // [id, merge_log, create_meta_log, filter_log]
                    .join(mbarq_log_ch)         // [id, ..., mbarq_log]
                    .map { id, merge_log, meta_log, filter_log, mbarq_log ->
                        tuple(id, [merge_log, meta_log, filter_log, mbarq_log])
                    }

    master_log_step = merge_logs(ordered_logs)
    master_log_ch   = master_log_step.master_log_ch

  emit:
    mbarq_meta_ch
}
