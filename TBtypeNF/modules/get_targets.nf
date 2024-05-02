
process GET_TARGETS {
    cpus 1
    memory '1 GB'
    time '1 h'
    label 'TBtypeR'

    input:
    tuple path(ref), path(ref_idx)
    path(regions)
    path(panel)
    path(dr_panel)

    output:
    path 'snp_targets.tsv.gz', emit: 'snp'
    path 'panel_sites.tsv.gz', emit: 'panel'
    path 'indel_regions.bed.gz', emit: 'indel', optional: true

    script:
    regions_  = regions.name  == 'NO_REGIONS'  ? 'NULL' : "'$regions'"
    dr_panel_ = dr_panel.name == 'NO_DR_PANEL' ? 'NULL' : "'$dr_panel'"
    """
    R --slave --vanilla -e "
        TBtypeR::create_calling_targets(
            ref_fasta       = '$ref',
            panel           = '$panel',
            dr_panel        = $dr_panel_,
            calling_regions = $regions_,
            output_dir      = '.')
        readr::write_tsv(
            dplyr::select(TBtypeR::read_panel('$panel'), chrom, pos),
            'panel_sites.tsv.gz',
            col_names = F)
        "
    """
}