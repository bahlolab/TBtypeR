
params.min_qual = 20

process bcftools_call {
    cpus 2
    memory '2 GB'
    time '2 h'
    label 'bcftools'
    tag { "$sm" }

    input:
    tuple val(sm), path(bam), path(bai), path(ref), path(ref_files),
          path(regions)

    output:
        tuple val(sm), path(out)

    script:
    out = "${sm}.vcf.gz"
    """
    bcftools mpileup $bam -B -T $regions -I -q 1 -Q 10 -d 200 -f $ref -a FMT/AD,FMT/DP --threads ${task.cpus} -Ou |
        bcftools call -m -Ou |
        bcftools norm -m-any -Ou |
        bcftools view -i "GT='alt' & QUAL>$params.min_qual & FMT/DP>=10" -Oz -o $out
    """
}

process bcftools_call_old {
    cpus 2
    memory '2 GB'
    time '2 h'
    container 'quay.io/biocontainers/bcftools:1.7--0'
    tag { "$sm" }

    input:
    tuple val(sm), path(mpileup), path(regions)

    output:
    tuple val(sm), path(out)

    script:
    out = "${sm}.vcf.gz"
    """
    bcftools call $mpileup -T $regions --threads ${task.cpus} -c -Ou |
        bcftools norm -m-any -Ou |
        bcftools view -i "GT='alt' & QUAL>20" -Oz -o $out
    """
}
