
process fastlin {
    cpus 2
    memory { "${ 4 + 4 * task.attempt } GB" }
    time { task.attempt + ' h' }
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    tag { sm }
    container 'bahlolab/tbmixbench:fastlin'

    input:
        tuple val(sm), path(f1), path(f2)
        path(barcode)

    output:
        tuple val(sm), path("$out")

    script:
    out = "${sm}_fastlin.tsv"
    """
    fastlin -d . -b ${barcode.name != 'DEFAULT' ? "$barcode" : '/fastlin/MTBC_barcodes.tsv'} -o $out
    sed 's:#sample:sample:g' -i $out
    """
}
