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

def buildComparisonList(List<List> tuples, Map comparisons) {
    // tuples: [ [sid1, path1], [sid2, path2], ... ]
    def sampleIndex = buildSampleIndex(tuples)
    def allIds      = sampleIndex.keySet() as Set<String>

    def requiredKeys = [
        'name', 
        'treatments',
        'controls'
    ]

    def allowedKeys = [
        'name',
        'treatments',
        'controls',
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

        def treatSel    = (cmp.treatments ?: []) as List<String>
        def controlSel  = (cmp.controls   ?: []) as List<String>

        def treatNegSel   = (cmp.treatments_negative_selection ?: []) as List<String>
        def controlNegSel = (cmp.controls_negative_selection   ?: []) as List<String>

        def treatIds   = resolveSelectors(treatSel, allIds, treatNegSel)
        def controlIds = resolveSelectors(controlSel, allIds, controlNegSel)

        tuple(
            name,
            treatIds.collect   { sid -> [sid, sampleIndex[sid]] },
            controlIds.collect { sid -> [sid, sampleIndex[sid]] }
        )
    }
}

