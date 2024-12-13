
process MOSDEPTH {
    cpus 2
    memory '2 GB'
    time '1 h'
    container 'quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2'
    tag { sm }

    input:
    tuple val(sm), path(bam), path(bai), path(bed), path(ref), path(ref_idx)

    output:
    tuple val(sm), path("${sm}.mosdepth.*")

    script:
    //TODO - check non-crams unaffected
    if (bed.name == 'NO_REGIONS')
        """
        mosdepth $sm $bam --threads ${task.cpus} --fast-mode --no-per-base --mapq 1 --fasta $ref
        """
    else
        """
        mosdepth $sm $bam --by $bed --threads ${task.cpus} --fast-mode --no-per-base --mapq 1 --fasta $ref
        rm ${sm}.regions.bed.gz*
        """
}