
process FASTLIN {
    cpus 2
    memory { "${ 4 + 4 * task.attempt } GB" }
    time { task.attempt + ' h' }
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    tag { sm }
    container 'bahlolab/tbmixbench:fastlin'

    input:
        path(barcodes)
        tuple val(sm), path(f1), path(f2)

    output:
        tuple val(sm), path("${out}.gz")

    script:
    out = "${sm}_fastlin.tsv"
    """
    fastlin -d . -c 1 -b $barcodes -o $out
    gzip $out
    """
}
