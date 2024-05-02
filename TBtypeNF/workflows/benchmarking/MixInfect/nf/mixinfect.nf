
params.low_cov = 10
params.min_qual = 20

process mixinfect {
    cpus 1
    memory '2 GB'
    time { task.attempt + ' h'}
    container 'bahlolab/tbmixbench:mixinfect'
    // container 'bahlolab/mixinfect:latest'
    tag {sm}

    input:
    tuple val(sm), path(vcf), path(bam), path(bai)

    output:
    path(out)

    script:
    out = "MixInfect.${sm}.csv"
    // if MixInfect does find enough heterozygous SNPs it will fail - but this is equivalent to not detecting a mixture
    failure_output = "\"Sample name,Mix or Non-mix,Mixed SNPs,Total SNPs,Proportion het/total SNPs,Number of strains in mix,Major strain proportion,SD Major strain proportion,SEM Major strain proportion,CI low Major strain proportion,CI high Major strain proportion,No different genes,Proportion hetSNPs to genes\n${sm},Non-Mix,NA,NA,NA,1,1,NA,NA,NA,NA,NA,NA\""
    """
    MixInfectFB.R $vcf mi $params.min_qual $params.low_cov &&
        sed 's:"::g' mi_output_genes_combined.csv > $out ||
        sed 's:"::g' mi_mixes_separate.csv > $out ||
        echo $failure_output > $out
    """
}