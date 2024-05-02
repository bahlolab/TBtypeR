
process SAMTOOLS_STATS {
    cpus 2
    memory '2 GB'
    time '1 h'
    label 'align'
    tag { sm }

    input:
        tuple val(sm), path(bam), path(bai)

    output:
        tuple val(sm), path(out)

    script:
    out = "${sm}.samtools_stats"
    """
    samtools stats $bam -@ ${task.cpus} > $out
    """
}
