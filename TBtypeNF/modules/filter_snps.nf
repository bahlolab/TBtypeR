
process FILTER_SNPS {
    cpus 2
    memory '2 GB'
    time '1 h'
    label 'bcftools'

    input:
        tuple path(bcf), path(sites)

    output:
        path(out_bcf)

    script:
        out_bcf = "${params.clean_id}.snps.filt.bcf"
        """
        bcftools view $bcf --threads ${task.cpus} -i 'GT="AA"' -Ou |
            bcftools +fill-tags -Ou -- -t AF |
            bcftools view -i 'AF>=$params.min_af' -Ou |
            bcftools annotate -x INFO/AF -Ob -o filt.bcf
        bcftools index --threads ${task.cpus} filt.bcf
            
        bcftools view $bcf --threads ${task.cpus} -T $sites -Ob -o sites.bcf
        bcftools index --threads ${task.cpus} sites.bcf
        
        bcftools concat sites.bcf filt.bcf --threads ${task.cpus} -a -D -Ob -o $out_bcf        
        rm filt.bcf* sites.bcf*
        """
}