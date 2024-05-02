
process split_strains {
    cpus 1
    memory '1 GB'
    time '1 h'
    tag { sample }
    container 'bahlolab/tbmixbench:splitstrains'

    input:
        tuple val(sample), path(bam), path(bai), val(depth), path(ref)

    output:
        path(out)

    script:
    out = "${sample}.splitStrains.txt"
    """
    splitStrains $bam \\
        -b /SplitStrains/refs/tuberculosis.filtered.gff \\
        -mo gmm \\
        -fe 0.70 \\
        -fes 60 \\
        -g 2 \\
        -f $sample \\
        -s 50000 \\
        -e 4400000 \\
        -r $ref \\
        -o . \\
        -fd $depth |
        tee /dev/stderr |
        grep -v 'INFO:\\|WARNING:' > $out
    """
}