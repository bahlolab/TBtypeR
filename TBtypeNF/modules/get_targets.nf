
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
        
        panel_sites <- TBtypeR::read_panel('$panel')
        try(panel_sites <- dplyr::bind_rows(panel_sites, TBtypeR::read_panel($dr_panel_, phylo = FALSE)))
        panel_sites <- dplyr::filter(panel_sites, nchar(ref) == 1, nchar(alt) == 1)
        panel_sites <- dplyr::select(panel_sites, chrom, pos)
        panel_sites <- dplyr::arrange(panel_sites, chrom, pos)
        readr::write_tsv(
            dplyr::distinct(panel_sites),
            'panel_sites.tsv.gz',
            col_names = F)
        "
    """
}