
process BCFTOOLS_MERGE {
    cpus    2
    memory '4 GB'
    time   '1 h'
    label  'bcftools'
    tag    "$id"

    input:
    tuple val(id), path(snps), path(indels), path(ref), path(ref_files), path(sites)

    output:
    tuple val(id), path(snps_out), path(indels_out)

    script:
    snps_out = "snps.merged.${id}.bcf"
    indels_out = indels[0].name == 'NO_INDELS' ? 'NO_INDELS' : "indels.merged.${id}.bcf"
    indel_cmd  = indels[0].name == 'NO_INDELS' ? '' :
    """
    bcftools merge $indels --threads ${task.cpus} --no-index -Ob -o $indels_out
    """
    """
    bcftools merge $snps --threads ${task.cpus} --no-index -Ou |
        bcftools view -T $sites -Ob -o $snps_out
    
    """ + indel_cmd
}
