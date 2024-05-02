
process seqkit_sample {
    cpus 2
    memory '1 GB'
    time '1 h'
    container 'quay.io/biocontainers/seqkit:2.2.0--h9ee0642_0'
    tag { "$id:$sam" }

    input:
        tuple val(id), val(sam), val(frac), val(seed), path(fq1), path(fq2)

    output:
        tuple val(id), path(out1), path(out2)

    script:
    out1 = "${id}_${sam}_1.fq.gz"
    out2 = "${id}_${sam}_2.fq.gz"
    """
    zcat $fq1 | seqkit sample -p $frac -s $seed -o $out1 &
    zcat $fq2 | seqkit sample -p $frac -s $seed -o $out2
    wait
    """
}