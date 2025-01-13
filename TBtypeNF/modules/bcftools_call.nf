
process BCFTOOLS_CALL {
    cpus   2
    memory '2 GB'
    time   '2 h'
    tag    "$sm" 
    label  'bcftools'

    input:
    tuple val(sm), path(bam), path(bai), path(ref), path(ref_files),
          path(snp_targets), path(indel_regions)

    output:
    path alt_sites,                                emit: 'sites'
    tuple val(sm), path(snp_bcf), path(indel_bcf), emit: 'calls'
    

    script:
    snp_bcf   = "${sm}.snps.bcf"
    alt_sites = "${sm}.snp_alt.gz"
    indel_bcf = indel_regions.name == 'NO_INDELS' ? 'NO_INDELS' : "${sm}.indels.bcf" 
    indel_cmd = indel_regions.name == 'NO_INDELS' ? '' :
    """
    bcftools mpileup $bam --threads ${task.cpus} -q 1 -Q 10 -d 200 -f $ref -R $indel_regions -a FMT/AD -Ou |
        bcftools call -V snps -A -m --prior 1e-2 -Ou |
        bcftools norm -m-any -Ou |
        bcftools annotate -x INFO,^FORMAT/GT,^FORMAT/AD -Ob -o $indel_bcf
    """
    """
    bcftools mpileup $bam --threads ${task.cpus} -q 1 -Q 10 -d 200 -f $ref -a FMT/AD -Ou |
        bcftools call -A -m --prior 1e-2 -C alleles -T $snp_targets -Ou |
        bcftools annotate -x INFO,^FORMAT/GT,^FORMAT/AD -Ob -o $snp_bcf

    bcftools view $snp_bcf --threads ${task.cpus} -i 'GT="AA"' -Ou | 
        bcftools query -f "%CHROM\\t%POS\\n" |
        gzip > $alt_sites
    
    """ + indel_cmd
}
