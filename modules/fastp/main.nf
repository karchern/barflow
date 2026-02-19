process run_fastp_process {
    tag { sample_id }
    cpus 4
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz") into trimmed_ch
    path "${sample_id}.fastp.json" , emit: fastp_json
    path "${sample_id}.fastp.html"

    script:
    def out_fastq = "${sample_id}.trimmed.fastq.gz"
    def json = "${sample_id}.fastp.json"
    def html = "${sample_id}.fastp.html"
    """
    fastp -i ${reads} -o ${out_fastq} -j ${json} -h ${html} -w ${task.cpus}

    # also print JSON summary to stdout for Nextflow logs
    echo "fastp summary for ${sample_id}:"
    cat ${json}
    """
}

workflow run_fastp (reads_ch, params) {
    reads_ch.into { in_reads }
    in_reads.map { tuple(it[0], it[1]) } | run_fastp_process
    return trimmed_ch
}
