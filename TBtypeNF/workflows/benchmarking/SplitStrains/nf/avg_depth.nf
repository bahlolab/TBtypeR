
process avg_depth {
    cpus 1
    memory '1 GB'
    time '1 h'
    tag { sample }
    container 'bahlolab/tbmixbench:splitstrains'

    input:
        tuple val(sample), path(bam), path(bai)

    output:
        tuple val(sample), path(bam), path(bai), path(depth)

    script:
    depth = sample + '.depth'
    """
    samtools depth $bam |
        awk '{ total += \$3; count++ } END { print int(total/count) }' \\
        > $depth
    """
}