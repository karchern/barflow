nextflow.enable.dsl=2

process count_barcodes {
    //tag { sample_id }
    label 'fast2q'
    publishDir "${params.outdir}/2fast2q", mode: 'copy'
    errorStrategy 'terminate' // This must absolutely be terminate, as otherwise logging/tracing of samples is fucked

    input:
    tuple val(sample_id), path(fastq_path), path(good_barcodes_csv_path), val(library), val(lib_cfg)
    

    output:
    tuple val(sample_id), path("${sample_id}.2fast2q")

    script:
    """
    # put the staged fastq into its own directory (symlink so no big copy)
    # Necessary since 2fast2q expects a directory of FASTQ files, but we want to process one sample at a time for better parallelism and resource management.
    fname=\$(basename "${fastq_path}")
    mkdir -p input
    ln -s "\$PWD/\$fname" input/"\$fname"

    # Check that ${good_barcodes_csv_path} exists, is not empty, and is a csv file
    if [ ! -s "${good_barcodes_csv_path}" ]; then
        echo "Error: ${good_barcodes_csv_path} does not exist or is empty."
        exit 1
    fi
    if [[ "${good_barcodes_csv_path}" != *.csv ]]; then
        echo "Error: ${good_barcodes_csv_path} is not a CSV file."
        exit 1
    fi
    # Also check that it literally has "," characters in it, and the same number per line across the entire file
    num_commas=\$(head -n 1 "${good_barcodes_csv_path}" | awk -F',' '{print NF-1}')
    if [ \$num_commas -eq 0 ]; then
        echo "Error: ${good_barcodes_csv_path} does not contain any commas, and thus is not a valid CSV file."
        exit 1
    fi
    if ! awk -F',' -v num_commas="\$num_commas" 'NF-1 != num_commas { exit 1 }' "${good_barcodes_csv_path}"; then
        echo "Error: ${good_barcodes_csv_path} does not have the same number of commas on every line, and thus is not a valid CSV file."
        exit 1
    fi



    cat ${good_barcodes_csv_path}| cut -d "," -f1 > tmp
    paste tmp tmp > good_barcodes_in_stupid_format.csv

    echo ${lib_cfg.upstream}
    echo ${lib_cfg.downstream}

    # run 2fast2q on that single-file folder
    2fast2q -c \\
    --s input \\
    --g good_barcodes_in_stupid_format.csv \\
    --o "${sample_id}.2fast2q" \\
    --ph 1 \\
    --us "${lib_cfg.upstream}" \\
    --ds "${lib_cfg.downstream}" \\
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
//     ch.map { it } | count_barcodes
//     return sample_counts_ch
// }

workflow barcode_counter {

    take:
    ch   // emits: [ sample_id, fastq_path, good_barcodes_csv_path, library ]

    main:
    validated = ch.map { sample_id, fastq, barcodes, library ->

        def lib_cfg = params.library_seqs[ library ]

        if( lib_cfg == null ) {
            // fail early with a clear message
            def known = params.library_seqs.keySet().sort().join(', ')
            throw new IllegalArgumentException(
                "Unknown library '${library}'. Known libraries: ${known}"
            )
        }

        // append lib_cfg as an extra field in the tuple
        [ sample_id, fastq, barcodes, library, lib_cfg ]
    }
    
    result = count_barcodes(validated)

    emit:
    result
}
