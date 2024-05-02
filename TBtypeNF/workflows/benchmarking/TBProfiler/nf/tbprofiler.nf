
process tbprofiler {
    cpus 4
    memory { "8 GB" }
    time { task.attempt + ' h' }
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    tag { sm }
    container 'bahlolab/tbmixbench:tbprofiler'
    // Note - the following container options only work for Singularity
    containerOptions {
        barcode.name == 'DEFAULT' ? 
            null :
            "-B \$PWD/$barcode:/opt/conda/share/tbprofiler/tbdb.barcode.bed"
    }

    input:
        tuple val(sm), path(f1), path(f2)
        path(barcode)

    output:
        path("$out")

    script:
    out = "${sm}.tbprofiler.json"
    """
    tb-profiler profile -1 $f1 -2 $f2 \\
       --no_delly \\
       --threads $task.cpus \\
       --ram $task.memory.giga
    rm bam vcf -r
    jq '. | { sample:"$sm",  sublin: .sublin, lineage: [.lineage[].lin],  frac: [.lineage[].frac ]} ' results/tbprofiler.results.json \\
        > $out
    """
}
