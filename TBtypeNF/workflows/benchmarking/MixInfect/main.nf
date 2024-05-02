nextflow.enable.dsl=2

import static Helpers.path
import static Helpers.read_tsv

params.manifest = 'manifest.tsv'
params.ref_fasta = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2/GCA_000195955.2_ASM19595v2_genomic.fna.gz'
params.mi_calling_regions = "${workflow.projectDir}/misc/calling_regions.bed.gz"

include { trimmomatic   } from './nf/trimmomatic'
include { align         } from './nf/align'
include { bcftools_call } from './nf/bcftools_call'
include { mixinfect     } from './nf/mixinfect'

workflow MixInfect {
    take: opts

    main:

    ref_ch = opts.ref ?: index_ref(path(params.ref_fasta))

    fastqs = Channel.fromList(
        read_tsv(path(params.manifest),
                 required:['sample', 'fastq1', 'fastq2'])) |
        map { [it.sample, path(it.fastq1), path(it.fastq2)] }

    bams = fastqs |
        trimmomatic |
        combine(ref_ch) |
        align

    results = bams |
        combine(ref_ch) |
        combine([path(params.mi_calling_regions)]) |
        bcftools_call |
        combine(bams, by:0) |
        mixinfect |
        collectFile(name: 'MixInfect.results.tsv',
                    storeDir: params.outdir,
                    keepHeader:true)

    emit: results
}

workflow {
    MixInfect([:])
}