def buildSampleIndex(List tuples) {
    // tuples are [sample_id, Path]
    def idx = [:]
    tuples.each { sid, p -> idx[sid as String] = p }
    return idx
}

def resolveSelectors(List<String> selectors, Set<String> allIds) {
    def result = [] as Set   // unique
    selectors.each { sel ->
        if( sel.contains('*') || sel.contains('?') ) {
            // treat as glob; translate * and ? to regex
            def regex = '^' + sel
                .replace('.', '\\.')
                .replace('*', '.*')
                .replace('?', '.') + '$'
            //def pattern = regex as java.util.regex.Pattern
            def pattern = java.util.regex.Pattern.compile(regex)
            allIds.findAll { sid -> pattern.matcher(sid).matches() }
                  .sort()
                  .each { result << it }
        }
        else {
            // explicit ID
            if( allIds.contains(sel) )
                result << sel
        }
    }
    return result.toList().sort() as List
}

def buildComparisonList(List<List> tuples, Map comparisons) {
    // tuples: [ [sid1, path1], [sid2, path2], ... ]
    def sampleIndex = buildSampleIndex(tuples)
    def allIds      = sampleIndex.keySet() as Set<String>

    //println "Processing comparisons with sample index: ${sampleIndex}"
    //println comparisons.comparisons

    comparisons.comparisons.collect { cmp ->
        def name       = cmp.name as String
            def treatSel   = (cmp.treatments ?: []) as List<String>
            def controlSel = (cmp.controls   ?: []) as List<String>

            def treatIds   = resolveSelectors(treatSel, allIds)
            def controlIds = resolveSelectors(controlSel, allIds)

            
            tuple(
                name,
                treatIds.collect { sid -> [sid, sampleIndex[sid]] },
                controlIds.collect { sid -> [sid, sampleIndex[sid]] }
            )
    }
}
