
process quanttb {
    cpus 2
    memory { "${ 4 + 4 * task.attempt } GB" }
    time { task.attempt + ' h' }
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    tag { sm }
    container 'bahlolab/tbmixbench:quanttb'

    input:
        tuple val(sm), path(f1), path(f2)
        path barcode

    output:
        tuple val(sm), path("$out")

    script:
    out = "${sm}.out.txt"
    """
    export PYTHON_EGG_CACHE="\$PWD/python_egg_cache"
    quanttb quant -f $f1 $f2 -o $out ${barcode.name != 'DEFAULT' ? "-db $barcode" : ''}
    mv output/$out .
    rm python_egg_cache -r
    """
}
