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

workflow print_summary_table_p {

    take:
    all_counts_list_ch
    sample_qcs_PASSED_CH
    comparisons_status_ch
    filter_ch_comparison_post_filter

    main:
    def total_samples_ch = all_counts_list_ch.map { all_samples_list ->
        all_samples_list.size()
    }

    def samples_after_qc_ch = sample_qcs_PASSED_CH.count()

    def comps_before_ch = comparisons_status_ch.map { x ->
        x.size()
    }

    def comps_after_ch = filter_ch_comparison_post_filter.count()

    print_summary_table(
        total_samples_ch,
        samples_after_qc_ch,
        comps_before_ch,
        comps_after_ch
    )
}

process collate_barseq_qc_results {
    
    label 'python_basic'

    publishDir "${params.outdir}/logs/", mode: 'copy', overwrite: true

    input:
    path metrics_file_paths

    output:
    path("all_samples.barcode_metrics.tsv"), emit: sample_wise_qc_metrics

    script:
    """
    echo -e "sample_id,total_reads,n_barcodes,max_count,mean_count,median_count" > all_samples.barcode_metrics.tsv
    for f in *.csv; do
        tail -n +2 "\$f" >> all_samples.barcode_metrics.tsv
    done
    """

}

process barseq_qc {

    label 'python_basic'

    input:
    tuple val(sample_id), path(counts_path)
    val(minimum_read_sum_for_qc)

    output:
    tuple val(sample_id), path("${sample_id}.barcode_metrics.csv"), path("${sample_id}.passed_qc.txt"), emit: metrics

    """
    barseq_qc.py \
        --sample_id ${sample_id} \
        --input ${counts_path} \
        --output ${sample_id}.barcode_metrics.csv \
        --output_passed ${sample_id}.passed_qc.txt \
        --min_read_sum_for_qc ${minimum_read_sum_for_qc}
    
    """
}

process print_summary_table {

    cache false
    debug true

    input:
    val total_samples
    val samples_after_qc
    val comps_before
    val comps_after

    exec:
    println """
        ============================================
                    Pre-mbarq numbers
        ============================================
        | ${"Metric".padRight(28)} | ${"Number".padLeft(8)} |
        --------------------------------------------
        | ${"BarSeq samples before QC".padRight(28)} | ${String.valueOf(total_samples).padLeft(8)} |
        | ${"BarSeq samples after QC".padRight(28)} | ${String.valueOf(samples_after_qc).padLeft(8)} |
        | ${"Comparisons before filtering".padRight(28)} | ${String.valueOf(comps_before).padLeft(8)} |
        | ${"Comparisons after filtering".padRight(28)} | ${String.valueOf(comps_after).padLeft(8)} |
        ============================================
        """.stripIndent()
}