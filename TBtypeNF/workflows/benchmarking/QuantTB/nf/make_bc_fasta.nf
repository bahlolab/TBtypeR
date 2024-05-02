
process make_bc_fasta {
    cpus 1
    memory '1 GB'
    time '1 h'
    label 'TBtypeR'

    input:
    tuple path(ref), path(ref_idx)
    path(panel)

    output:
    path("*.fna")

    script:
    """
    R --slave --vanilla -e "
        TBtypeR:::make_quanttb_fasta(
            ref_fasta = '$ref',
            panel     = '$panel',
            output    = '.')
        "
    """
}