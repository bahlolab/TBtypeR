
process align {
    cpus 4
    memory '4 GB'
    time '2 h'
    tag { sm }
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:8110a70be2bfe7f75a2ea7f2a89cda4cc7732095-0'

    input:
        tuple val(sm), path(f1), path(f2), path(ref), path(ref_files)

    output:
        tuple val(sm), path(output), path("${output}.bai")

    script:
    output= sm + '.bam'
    """
    bwa aln -q 15 -t ${task.cpus} -R '@RG\\tID:mixed\\tSM:mixed\\tLB:None\\tPL:Illumina' $ref $f1 > R1.sai
    bwa aln -q 15 -t ${task.cpus} -R '@RG\\tID:mixed\\tSM:mixed\\tLB:None\\tPL:Illumina' $ref $f2 > R2.sai
    bwa sampe $ref R1.sai R2.sai $f1 $f2 |
        samtools sort -@ ${task.cpus} -n -u |
        samtools fixmate -@ ${task.cpus} -m -r - fixmate.bam
    samtools sort fixmate.bam -@ ${task.cpus} -u |
        samtools markdup -@ ${task.cpus} -r - $output
    samtools index -@ ${task.cpus} $output
    rm R1.sai R2.sai fixmate.bam
    """
}