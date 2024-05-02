
process trimmomatic {
    cpus 2
    memory '2 GB'
    time '1 h'
    tag { sample }
    container 'quay.io/biocontainers/trimmomatic:0.36--6'

    input:
        tuple val(sample), path(fq1), path(fq2)

    output:
        tuple val(sample), path("${pref}_1P.fq.gz"), path("${pref}_2P.fq.gz")

    script:
    pref = sample + '.trimm'
    """
    trimmomatic PE -threads ${task.cpus} $fq1 $fq2 -baseout ${pref}.fq.gz \\
         LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """
}