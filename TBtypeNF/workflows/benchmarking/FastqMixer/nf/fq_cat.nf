
process fq_cat {
    cpus 2
    memory '1 GB'
    time '1 h'
    tag { id }
    publishDir 'output', mode: 'copy'

    input:
        tuple val(id), path(fq1), path(fq2)

    output:
        tuple val(id), path(out1), path(out2)

    script:
    out1 = "${id}_1.fq.gz"
    out2 = "${id}_2.fq.gz"
    if (fq1 instanceof List)
        """
        zcat $fq1 | gzip > $out1 &
        zcat $fq2 | gzip > $out2
        wait
        """
    else
        """
        cp $fq1 $out1 &
        cp $fq2 $out2 &
        """
}