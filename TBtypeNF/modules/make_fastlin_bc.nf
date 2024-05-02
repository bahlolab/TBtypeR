
process MAKE_FASTLIN_BC {
    cpus 1
    memory '1 GB'
    time '1 h'
    label 'TBtypeR'

    input:
    tuple path(ref), path(ref_idx)
    path(panel)

    output:
    path(output)

    script:
    output = 'TBtypeR.fastlin_bc.tsv'
    """
    R --slave --vanilla -e "
        TBtypeR:::make_fastlin_bc(
            ref_fasta = '$ref',
            panel     = '$panel',
            output    = '$output')
        "
    """
}