
process INDEX_REF {
    cpus 1
    memory '1 GB'
    time '1 h'
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:8110a70be2bfe7f75a2ea7f2a89cda4cc7732095-0'

    input:
    path(input)

    output:
    tuple path('ref.fa.gz'), path('ref.fa.gz.*')

    script:
    """
    gzip -t $input && gzip -cd $input | bgzip > ref.fa.gz || 
        bgzip $input > ref.fa.gz
    samtools faidx ref.fa.gz
    bwa index ref.fa.gz
    """
}