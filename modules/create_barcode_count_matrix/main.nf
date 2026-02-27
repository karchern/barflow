nextflow.enable.dsl=2

process extract_barcodes_process {
    tag { sample_id }
    publishDir "${params.outdir}/2fast2q", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_path)
    path good_barcodes_csv_path

    output:
    tuple val(sample_id), path("${sample_id}.2fast2q")

    script:
    """
    # put the staged fastq into its own directory (symlink so no big copy)
    # Necessary since 2fast2q expects a directory of FASTQ files, but we want to process one sample at a time for better parallelism and resource management.
    fname=\$(basename "${fastq_path}")
    mkdir -p input
    ln -s "\$PWD/\$fname" input/"\$fname"

    cat ${good_barcodes_csv_path}| cut -d "," -f1 > tmp
    paste tmp tmp > good_barcodes_in_stupid_format.csv

    # run 2fast2q on that single-file folder
    2fast2q -c \\
    --s input \\
    --g good_barcodes_in_stupid_format.csv \\
    --o "${sample_id}.2fast2q" \\
    --ph 1 \\
    --us "${params.upstream_seq}" \\
    --ds "${params.downstream_seq}" \\
    --msu 3 \\
    --msd 3 \\
    --qsu 1 \\
    --qsd 1 \\
    --mo C \\
    --m 2 \\
    --cp ${task.cpus}

    # In order to comply with the output path being emitted, we need to move the resulting counts file to the expected location
    mv ${sample_id}.2fast2q TMP
    find . -name "compiled.csv" | xargs -i mv {} ${sample_id}.2fast2q

    """
}

// workflow create_counts (ch, params) {
//     ch.map { it } | extract_barcodes_process
//     return sample_counts_ch
// }

workflow create_counts {

  take:
  ch
  good_barcodes_ch

  main:
  result = extract_barcodes_process(ch, good_barcodes_ch)

  emit:
  result
}
