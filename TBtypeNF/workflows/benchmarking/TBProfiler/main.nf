#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import static Helpers.*
import groovy.json.JsonSlurper

params.manifest = 'manifest.tsv'
params.outdir = 'output'

include { tbprofiler         } from './nf/tbprofiler'
include { make_tbprofiler_bc } from './nf/make_tbprofiler_bc'


workflow TBProfiler {
    main:

    if (params.tbprofiler_use_tbt_bc) {
        panel_ch = Channel.value(path(params.panel))
        barcode_ch = make_tbprofiler_bc(panel_ch)
    } else {
        barcode_ch = Channel.value(path("$projectDir/resources/DEFAULT"))
    }

    fastqs = Channel.fromList(
        read_tsv(path(params.manifest),
                 required:['sample', 'fastq1', 'fastq2'])) |
        map { [it.sample, path(it.fastq1), path(it.fastq2)] }

    results = tbprofiler(fastqs, barcode_ch) |
        map {
            def data = new JsonSlurper().parseText(it.text)
            if (data.sublin == "") {
                [data.sample, "ERROR: ${data.lineage.join(',')}; ${data.frac.join(',')} ", "NA"].join('\t')
            } else {
                def lins = data.sublin.split(';') as ArrayList
                def fracs = data.frac[
                    lins.collect {  lin ->
                        data.lineage.withIndex().find { l, i -> l == lin}[1]
                    }
                ]
                fracs = fracs.collect { it / fracs.sum() }
                lins.withIndex().collect { lin, i -> 
                    [data.sample, lin, fracs[i]].join('\t')
                }.join('\n')
            }
        } |
        collectFile(name: 'TBProfiler.results.tsv',
                    newLine:true,
                    sort:true,
                    seed: 'sample\tstrain\tproportion',
                    storeDir: params.outdir)

    emit: results
}

workflow {
    TBProfiler()
}
