
process VCF_TO_GDS {
    cpus 1
    memory '2 GB'
    time '1 h'
    publishDir "${params.outdir}/gds", mode: 'copy', enabled: params.publish_gds
    label 'TBtypeR'

    input:
        tuple path(vcf), path(idx)

    output:
        path(out)

    script:
    out = vcf.name.replaceAll('.vcf.gz', '.gds')
    """
    R --slave --vanilla -e "SeqArray::seqVCF2GDS('$vcf', '$out', ignore.chr.prefix = '')"
    """
}