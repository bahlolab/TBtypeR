#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import static Helpers.*

println "\n------------ nf-fq-mix ------------"
workflow.onComplete { println "Done." }

params.manifest = 'manifest.tsv'
params.outdir = 'output'

manifest = read_tsv(path(params.manifest),
                    required: ['id', 'sample', 'frac', 'seed', 'fq1', 'fq2'])
    .collect { it + [fq1: path(it.fq1), fq2:path(it.fq2)] }
    .collect { [it.id, it.sample, it.frac, it.seed, it.fq1, it.fq2 ]  }

include { seqkit_sample } from './nf/seqkit_sample'
include { fq_cat } from './nf/fq_cat'
file('output').mkdirs()

workflow {
    sampled =
        Channel.fromList(manifest) |
        seqkit_sample |
        groupTuple(by:0, sort:true) |
        fq_cat |
        map {
            it.collect {
                it instanceof Path ? file('./output/' + it.name.toString()).toString() : it
            }.join('\t')
        } |
        collectFile(name: 'output_manifest.tsv',
                    seed: 'id\tfastq1\tfastq2',
                    storeDir: params.outdir,
                    sort: true,
                    newLine: true)
}
