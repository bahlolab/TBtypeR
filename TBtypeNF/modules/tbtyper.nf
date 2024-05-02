
process TBTYPER {
    cpus 10
    memory '30 GB'
    time '2 h'
    publishDir params.outdir, mode: 'copy', enabled: params.publish
    label 'TBtypeR'

    input:
        tuple path(gds), path(panel), path(dr_panel), val(args_json), path(manifest)

    output:
        tuple path("${params.clean_id}.tbt_res.rds"),
              path("${params.clean_id}.tbt_calls.csv")
        path "${params.clean_id}.html", optional: true

    script:
    """
    TBType.R $gds $panel \\
        ${args_json ? "--args-json $args_json" : ''} \\
        ${dr_panel.name != 'NO_DR_PANEL' ? "--dr-panel $dr_panel" : ''} \\
        ${params.no_report ? "--no-report" : ''} \\
        --sample-meta $manifest \\
        --output $params.clean_id \\
        --threads ${task.cpus} \\
        --max-mix $params.max_mix
    """
}