#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import static Helpers.*

params.manifest = 'manifest.tsv'
params.outdir = 'output'

include { fastlin         } from './nf/fastlin'
include { make_fastlin_bc } from './nf/make_fastlin_bc'

workflow Fastlin {
    take: opts

    main:

    ref_ch = opts.ref ?: index_ref(path(params.ref_fasta))

    if (params.fastlin_use_tbt_bc) {
        panel_ch = Channel.value(path(params.panel))
        barcode_ch = make_fastlin_bc(ref_ch, panel_ch)
    } else {
        barcode_ch = Channel.value(path("$projectDir/resources/DEFAULT"))
    }

    fastqs = Channel.fromList(
        read_tsv(path(params.manifest),
                 required:['sample', 'fastq1', 'fastq2'])) |
        map { [it.sample, path(it.fastq1), path(it.fastq2)] }

    results = fastlin(fastqs, barcode_ch) |
        splitCsv(header:true, sep: '\t') |
        map { sm, res -> res + [sample:sm] } |
        map { it.values().join('\t') } |
        collectFile(name: 'fastlin.results.tsv',
                    newLine:true,
                    sort:true,
                    seed: 'sample\tdata_type\tk_cov\tmixture\tlineages\tlog_barcodes\tlog_errors',
                    storeDir: params.outdir)

    emit: results
}

workflow {
    Fastlin()
}
