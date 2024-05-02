
process make_fastlin_bc {
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
    output = 'TBtypeR_Fastlin_barcode.tsv'
    """
    R --slave --vanilla -e "
        TBtypeR:::make_fastlin_bc(
            ref_fasta = '$ref',
            panel     = '$panel',
            output    = '$output',
            mode      = 'fastlin')
        "
    """
}