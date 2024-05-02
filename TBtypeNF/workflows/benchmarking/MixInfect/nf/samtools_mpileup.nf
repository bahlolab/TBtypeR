
process samtools_mpileup {
    cpus 2
    memory '2 GB'
    time '2 h'
    container 'quay.io/biocontainers/samtools:1.7--0'
    tag { "$sm" }

    input:
    tuple val(sm), path(bam), path(bai), path(ref), path(ref_files),
          path(regions)

    output:
        tuple val(sm), path(out)

    script:
    out = "${sm}.mpileup.bcf"
    """
    samtools mpileup $bam -l $regions -t AD,DP -I -q 1 -f $ref -g -o $out
    """
}
