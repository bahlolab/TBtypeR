#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import static Helpers.*

params.manifest = 'manifest.tsv'
params.outdir = 'output'

include { quanttb } from './nf/quanttb'
include { make_bc_fasta } from './nf/make_bc_fasta'
include { make_bc_snps } from './nf/make_bc_snps'
include { make_snp_db } from './nf/make_snp_db'

workflow QuantTB {
    take: opts

    main:

    ref_ch = opts.ref ?: index_ref(path(params.ref_fasta))

    if (params.quanttb_use_tbt_bc) {
        panel_ch = Channel.value(path(params.panel))
        barcode_ch = make_bc_fasta(ref_ch, panel_ch) \
            | flatMap \
            | make_bc_snps \
            | collect \
            | make_snp_db 
    } else {
        barcode_ch = Channel.value(path("$projectDir/resources/DEFAULT"))
    }

    fastqs = Channel.fromList(
        read_tsv(path(params.manifest),
                 required:['sample', 'fastq1', 'fastq2'])) |
        map { [it.sample, path(it.fastq1), path(it.fastq2)] }

    results = quanttb(fastqs, barcode_ch) |
        splitCsv(header:true) |
        map { sm, res -> [sm, res + [sample:sm]  ] } |
        join (fastqs.map { [it[0]] }, remainder: true) |
        map { sm, res -> // No result, report all as NAs
            res != null ? res : [
                sample: sm,
                refname: 'NA',
                totscore: 'NA',
                relabundance: 'NA',
                depth: 'NA'
            ]
        } |
        map { it.values().join('\t') } |
        collectFile(name: 'quanttb.results.tsv',
                    newLine:true,
                    sort:true,
                    seed: 'sample\trefname\ttotscore\trelabundance\tdepth',
                    storeDir: params.outdir)

    emit: results
}

workflow {
    QuantTB()
}
