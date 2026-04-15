import groovy.json.JsonOutput
import groovy.json.JsonSlurper

def buildSampleIndex(List tuples) {
    // tuples are [sample_id, Path]
    def idx = [:]
    tuples.each { sid, p -> idx[sid as String] = p }
    return idx
}

def resolveComparisons(List<List> tuples, Map comparisons) {
    def samplesTsv = java.nio.file.Files.createTempFile('barflow_samples_', '.tsv').toFile()
    def comparisonsJson = java.nio.file.Files.createTempFile('barflow_comparisons_', '.json').toFile()
    def outputJson = java.nio.file.Files.createTempFile('barflow_resolved_', '.json').toFile()

    samplesTsv.deleteOnExit()
    comparisonsJson.deleteOnExit()
    outputJson.deleteOnExit()

    samplesTsv.text = tuples.collect { sid, p -> "${sid}\t${p}" }.join('\n') + '\n'
    comparisonsJson.text = JsonOutput.toJson(comparisons)

    def cmd = [
        'python3',
        "${projectDir}/bin/resolve_comparisons.py",
        '--samples-tsv', samplesTsv.absolutePath,
        '--comparisons-json', comparisonsJson.absolutePath,
        '--output-json', outputJson.absolutePath,
    ].collect { it as String }

    def proc = new ProcessBuilder(cmd)
        .redirectErrorStream(true)
        .start()

    def output = proc.inputStream.text
    def exitCode = proc.waitFor()
    if( exitCode != 0 ) {
        throw new IllegalStateException(
            "Python comparison resolver failed with exit code ${exitCode}. Output:\n${output}"
        )
    }

    return new JsonSlurper().parseText(outputJson.text)
}

def buildComparisonList(List<List> tuples, Map comparisons) {
    def sampleIndex = buildSampleIndex(tuples)
    def resolved = resolveComparisons(tuples, comparisons)

    resolved.collect { cmp ->
        def name = cmp.name as String
        def treatList = (cmp.treat_list as List).collect { sid, p -> [sid, sampleIndex[sid as String]] }
        def ctrlList = (cmp.ctrl_list as List).collect { sid, p -> [sid, sampleIndex[sid as String]] }

        tuple(
            name,
            treatList,
            ctrlList,
            cmp.good_barcodes_file,
            cmp.status,
            cmp.status_detail
        )
    }
}

// Decide and build all_counts_list_ch based on params
def createSampleInputChannelAndDecideIfToRun2Fast2Q(String samplesheet, String twofast2q_folder) {


    if ( samplesheet && !twofast2q_folder ) {

        // def reads_ch = Channel
        //     .fromPath(samplesheet)
        //     .splitCsv(header:false)
        //     .map { row ->
        //         def filename = file(row[0]).name.replaceFirst(/(?i)\.(fastq|fq)(\.gz)?|_sequence\.txt\.gz$/, '')
        //         tuple(filename, row[0])
        //     }

        def reads_ch = Channel
            .fromPath(samplesheet)
            .splitCsv(header:false)
            .map { row ->
                def readsPath = row[0] as String
                def resolvedReadsPath = readsPath.startsWith('/') ? readsPath : "${projectDir}/${readsPath}"
                def reads_file = file(resolvedReadsPath)
                def filename = reads_file.name.replaceFirst(/(?i)(\.(fastq|fq)(\.gz)?|_sequence\.txt\.gz)$/, '')
                tuple(filename, reads_file)
            }        


        return [ reads_ch: reads_ch, run_counts: true ]

    }
    else if( !samplesheet && twofast2q_folder ) {

        def reads_ch = Channel
            .fromPath("${twofast2q_folder}/*2fast2q")
            .map { Path p ->
                def sample_id = p.name.replaceFirst(/(?i)\.2fast2q$/, '')
                tuple(sample_id, p)
            }

    return [ reads_ch: reads_ch, run_counts: false ]
    }
    else {
        log.warn "\u001b[33mExactly one of --samplesheet or --twofast2q_folder must be set (not both, not neither).\u001b[0m"
        System.exit(1)
    }

}

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
    printf '%s\n' "comparison\tstatus\tstatus_detail" > comparisons_status_${type}.tsv
    printf '%b\n' "${lines}" >> comparisons_status_${type}.tsv
    """
}

def stop_after_barcode_extraction_and_warn(all_counts) {
    // materialize channel into list and emit a single view message when the flag is set
    if( params.stop_after_barcode_extraction as boolean ) {
        all_counts.toList()
            .view { tuples ->
                log.warn """
                ================================
                STOP AFTER BARCODE EXTRACTION
                ================================
                Parameter: --stop_after_barcode_extraction true

                Barcode extraction and count matrix generation have completed.
                No downstream merging or mbarq analysis will be performed.
                Disable this behavior by omitting the parameter or setting:
                    --stop_after_barcode_extraction false
                ================================
                """.stripIndent()
            }
        return true
    }
    return false
}