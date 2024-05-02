
process FASTP_PE {
    cpus 2
    memory '4 GB'
    time '1 h'
    tag { sm }
    label 'fastp'

    input:
    tuple val(sm), path(fq1), path(fq2)

    output:
    tuple val(sm), path("${sm}.fastp_R1.fq.gz"), path("${sm}.fastp_R2.fq.gz"),
          path("${sm}.fastp.json")

    script:
    """
    fastp -i $fq1 -I $fq2 -o ${sm}.fastp_R1.fq.gz -O ${sm}.fastp_R2.fq.gz \\
        --thread ${task.cpus} --compression 6 --json ${sm}.fastp.json
    sed 's:$fq1:$sm:g' -i ${sm}.fastp.json
    rm fastp.html
    """
}

process FASTP_SE {
    cpus 2
    memory '4 GB'
    time '1 h'
    tag { sm }
    label 'fastp'

    input:
    tuple val(sm), path(fq)

    output:
    tuple val(sm), path("${sm}.fastp.fq.gz"), path("${sm}.fastp.json")

    script:
    """
    fastp -i $fq -o ${sm}.fastp.fq.gz \\
        --thread ${task.cpus} --compression 6 --json ${sm}.fastp.json
    sed 's:$fq:$sm:g' -i ${sm}.fastp.json
    rm fastp.html
    """
}