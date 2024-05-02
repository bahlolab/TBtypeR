
process ALIGN_PE {
    cpus 4
    memory '4 GB'
    time '2 h'
    label 'align'
    tag { "$sm" }
    publishDir "${params.outdir}/bam", mode: 'copy', enabled: params.publish_bams

    input:
        tuple val(sm), path(f1), path(f2), path(ref), path(ref_files)

    output:
        tuple val(sm), path("${sm}.bam"), path("${sm}.bam.bai")

    script:
    """
    bwa mem -M -t ${task.cpus} -R '@RG\\tID:${sm}\\tSM:${sm}' $ref $f1 $f2 |
        samtools view -b > unsorted.bam
    samtools sort unsorted.bam -@ ${task.cpus} -n -u |
        samtools fixmate -@ ${task.cpus} -m - fixmate.bam
    samtools sort fixmate.bam -@ ${task.cpus} -u |
        samtools markdup -@ ${task.cpus} -r - ${sm}.bam
    samtools index -@ ${task.cpus} ${sm}.bam
    rm unsorted.bam fixmate.bam
    """
}

process ALIGN_SE {
    cpus 4
    memory '4 GB'
    time '2 h'
    label 'align'
    tag { "$sm" }
    publishDir "${params.outdir}/bam", mode: 'copy', enabled: params.publish_bams

    input:
    tuple val(sm), path(fq), path(ref), path(ref_files)

    output:
    tuple val(sm), path("${sm}.bam"), path("${sm}.bam.bai")

    script:
    """
    bwa mem -M -t ${task.cpus} -R '@RG\\tID:${sm}\\tSM:${sm}' $ref $fq |
        samtools view -b > unsorted.bam
    samtools sort unsorted.bam -@ ${task.cpus} -u |
        samtools markdup -@ ${task.cpus} -r - ${sm}.bam
    samtools index -@ ${task.cpus} ${sm}.bam
    rm unsorted.bam
    """
}