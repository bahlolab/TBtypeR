
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
    adapter = '/usr/local/share/trimmomatic-0.36-6/adapters/TruSeq3-PE-2.fa'
    """
    trimmomatic PE -phred33 -threads ${task.cpus} $fq1 $fq2 -baseout ${pref}.fq.gz \\
        ILLUMINACLIP:$adapter:3:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:16 MINLEN:40
    """
}