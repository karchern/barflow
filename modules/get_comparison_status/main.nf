// Module exporting the get_comparison_status process so it can be included multiple times
process get_comparison_status {

    publishDir "${params.outdir}/comparisons_status/", mode: 'copy', overwrite: true

    input:
    val comparisons_status_list
    val type

    output:
    path "comparisons_status_${type}.tsv"

    script:
    // Build the table as a single Groovy string
    def lines = comparisons_status_list.collect { name, status, detail ->
        "${name}\t${status}\t${detail ?: ''}"
    }.join('\\n')

    """
    printf '%s\n' "comparison_name\tstatus\tstatus_detail" > comparisons_status_${type}.tsv
    printf '%b\n' "${lines}" >> comparisons_status_${type}.tsv
    """
}
