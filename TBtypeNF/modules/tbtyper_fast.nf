
process TBTYPER {
    cpus 10
    memory '30 GB'
    time '2 h'
    publishDir params.outdir, mode: 'copy', enabled: params.publish
    label 'TBtypeR'
    tag { sample }

    input:
        path(panel)
        val(args_json)
        tuple val(sample), path(fastlin_out)

    output:
        tuple path("${sample}.tbt_res.rds"),
              path("${sample}.tbt_calls.csv")

    script:
    """
    TBTypeFast.R $fastlin_out $panel \\
        ${args_json ? "--args-json $args_json" : ''} \\
        --sample-id $sample \\
        --threads ${task.cpus} \\
        --max-mix $params.max_mix \\
        --min-mix-prop $params.min_mix_prop
    """
}