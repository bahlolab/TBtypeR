
process make_tbprofiler_bc {
    cpus 1
    memory '1 GB'
    time '1 h'
    label 'TBtypeR'

    input:
    path(panel)

    output:
    path(output)

    script:
    output = 'TBtypeR_TBProfiler_barcode.bed'
    """
    R --slave --vanilla -e "
        TBtypeR:::make_tbprofiler_bc(
            panel     = '$panel',
            output    = '$output')
        "
    """
}