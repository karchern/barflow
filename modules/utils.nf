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

def buildSampleIndex(List tuples) {
    // tuples are [sample_id, Path]
    def idx = [:]
    tuples.each { sid, p -> idx[sid as String] = p }
    return idx
}

// low-level: expand explicit + glob selectors to a Set of IDs
def matchSelectors(List<String> selectors, Set<String> allIds) {
    if( !selectors ) return [] as Set

    def result = [] as Set   // unique
    selectors.each { sel ->
        if( sel.contains('*') || sel.contains('?') ) {
            // treat as glob; translate * and ? to regex
            def regex = '^' + sel
                .replace('.', '\\.')
                .replace('*', '.*')
                .replace('?', '.') + '$'
            def pattern = java.util.regex.Pattern.compile(regex)
            allIds.findAll { sid -> pattern.matcher(sid).matches() }
                  .each { result << it }
        }
        else {
            // explicit ID
            if( allIds.contains(sel) )
                result << sel
        }
    }
    return result
}

// high-level: positive − negative, returns sorted List
def resolveSelectors(List<String> selectors,
                     Set<String> allIds,
                     List<String> negativeSelectors = []) {

    def pos = matchSelectors(selectors ?: [], allIds)
    def neg = matchSelectors(negativeSelectors ?: [], allIds)

    return (pos - neg).toList().sort() as List
}

// helper for buildComparisonList
def asStringList(value) {
    if( value == null )
        return [] as List<String>
    if( value instanceof List )
        return value.collect { it as String }
    // single scalar (String, GString, whatever): wrap in list
    return [ value as String ]
}

def buildComparisonList(List<List> tuples, Map comparisons) {
    // tuples: [ [sid1, path1], [sid2, path2], ... ]
    def sampleIndex = buildSampleIndex(tuples)
    def allIds      = sampleIndex.keySet() as Set<String>

    def requiredKeys = [
        'name', 
        'treatments',
        'controls',
        'good_barcodes_file',
    ]

    def allowedKeys = [
        'name',
        'treatments',
        'controls',
        "good_barcodes_file",
        'treatments_negative_selection',
        'controls_negative_selection'
    ] as Set

    comparisons.comparisons.collect { cmp ->
        // check for unexpected keys
        def cmpKeys      = (cmp as Map).keySet() as Set
        def invalidKeys  = cmpKeys - allowedKeys
        def missingKeys  = requiredKeys - cmpKeys
        if( missingKeys ) {
            throw new IllegalArgumentException(
                "Missing required fields in comparison '${cmp.name}': ${missingKeys.join(', ')}\nRequired fields are: ${requiredKeys.join(', ')}\nYou probably have a typo in your comparisons.json :)"
            )
        }        
        if( invalidKeys ) {
            throw new IllegalArgumentException(
                "Invalid fields in comparison '${cmp.name}': ${invalidKeys.join(', ')}\nAllowed fields are: ${allowedKeys.join(', ')}\nYou probably have a typo in your comparisons.json :)"
            )
        }

        def name        = cmp.name as String

        def treatSel    = asStringList(cmp.treatments)
        def controlSel  = asStringList(cmp.controls)

        def treatNegSel   = asStringList(cmp.treatments_negative_selection)
        def controlNegSel = asStringList(cmp.controls_negative_selection)

        def treatIds   = resolveSelectors(treatSel, allIds, treatNegSel)
        def controlIds = resolveSelectors(controlSel, allIds, controlNegSel)

        def treatNegIds = resolveSelectors(treatNegSel, allIds)
        def controlNegIds = resolveSelectors(controlNegSel, allIds)

        def good_barcodes_file = cmp.good_barcodes_file
        
        // compute missing IDs vs what was requested
        def requestedTreatIds   = treatSel
        def requestedControlIds = controlSel

        def requestedTreatNegIds   = treatNegSel
        def requestedControlNegIds = controlNegSel

        def missingTreat   = (requestedTreatIds - treatIds) as List
        def missingControl = (requestedControlIds - controlIds) as List

        // Remove treatNegSel from missingTreat 
        missingTreat = missingTreat - treatNegSel
        // Remove controlNegSel from missingControl
        missingControl = missingControl - controlNegSel

        def status
        def status_detail

        if( !missingTreat && !missingControl) {
            status        = 'OK'
            status_detail = ''
        }
        // if ALL requested treatments are missing, OR if ALL requested controls are missing, then we consider the comparison as FAILED (because it can't be run at all)
        else if( (requestedTreatIds && !treatIds) || (requestedControlIds && !controlIds) ) {
            status = 'ALL_SAMPLES_MISSING'
            def parts = []
            if( missingTreat )
                parts << "treatments: ${missingTreat.join(', ')}"
            if( missingControl )
                parts << "controls: ${missingControl.join(', ')}"
            status_detail = parts.join(' | ')
        }
        else {
            status = 'SOME_SAMPLES_MISSING'
            def parts = []
            if( missingTreat )
                parts << "treatments: ${missingTreat.join(', ')}"
            if( missingControl )
                parts << "controls: ${missingControl.join(', ')}"
            status_detail = parts.join(' | ')
        }
        tuple(
            name,
            treatIds.collect   { sid -> [sid, sampleIndex[sid]] },
            controlIds.collect { sid -> [sid, sampleIndex[sid]] },
            good_barcodes_file,
            status,
            status_detail
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
                def reads_file = file(row[0])
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