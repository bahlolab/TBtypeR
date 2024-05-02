#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import static Helpers.*

params.manifest = 'manifest.tsv'
params.ref_fasta = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2/GCA_000195955.2_ASM19595v2_genomic.fna.gz'
params.outdir = 'output'

include { index_ref     } from './nf/index_ref'
include { trimmomatic   } from './nf/trimmomatic'
include { align         } from './nf/align'
include { avg_depth     } from './nf/avg_depth'
include { split_strains } from './nf/split_strains'

workflow SplitStrains {
    take: opts

    main:

    ref = opts.ref ?: index_ref(path(params.ref_fasta))

    fastqs = Channel.fromList(
        read_tsv(path(params.manifest),
                 required:['sample', 'fastq1', 'fastq2'])) |
        map { [it.sample, path(it.fastq1), path(it.fastq2)] }

    results =
        fastqs |
        trimmomatic |
        combine(ref) |
        align |
        avg_depth |
        map { it[0..2] + [(0.7*(it[3].toFile().text as int)) as int]} |
        combine(ref.map {it[0]} ) |
        split_strains |
        collectFile(name: 'SplitStrains.results.tsv',
                    storeDir: params.outdir,
                    keepHeader:true)

    emit: results
}

workflow {
    SplitStrains([:])
}
