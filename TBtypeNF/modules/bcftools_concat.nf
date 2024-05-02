
process BCFTOOLS_CONCAT {
    cpus    2
    memory '2 GB'
    time   '1 h'
    label  'bcftools'
    publishDir params.outdir, mode: 'copy', enabled: params.publish

    input:
    tuple path(snp), path(indel)

    output:
    tuple path(out_vcf), path("${out_vcf}.*")

    script:
    out_vcf = "${params.clean_id}.concat.vcf.gz"
    if (indel.name == 'NO_INDELS')
        """
        bcftools view $snp -Oz -o $out_vcf
        bcftools index --threads ${task.cpus} $out_vcf
        """
    else
        """
        bcftools index --threads ${task.cpus} $snp
        bcftools index --threads ${task.cpus} $indel
        bcftools concat $snp $indel --threads ${task.cpus} --allow-overlaps -Oz -o $out_vcf
        bcftools index --threads ${task.cpus} $out_vcf
        """
}