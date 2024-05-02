
process MULTIQC {
    cpus 1
    memory '3 GB'
    time '1 h'
    container 'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    publishDir params.outdir, mode: 'copy', pattern: "*.html"

    input:
        tuple val(id), path(files)

    output:
        tuple path("${output}_data/multiqc_general_stats.txt"),
              path("${output}.html")

    script:
    output = "${params.clean_id}-${id}_multiqc_report"
    """
    multiqc . \\
        --title "$params.id - $id" \\
        --filename ${output}.html \\
        --module fastp \\
        --module samtools \\
        --module mosdepth \\
        --cl_config "max_table_rows: ${params.multiqc + 1}"
    """
}


