#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import static Helpers.path

/* ---------- TBtypeNF parameters ---------- */
params.manifest        = null
params.id              = 'TBtypeNF-run'
params.ref_fasta       = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2/GCA_000195955.2_ASM19595v2_genomic.fna.gz'
params.calling_regions = "$projectDir/resources/h37rv_core_genome.bed.gz" // Optional. Set NULL to skip SNP calling outside of panel/barcode
params.panel           = "$projectDir/resources/tbt_panel.tsv.bz2" // Optional. Set to null for whole genome.
params.min_af          = 0.001
params.max_mix         = 3
params.min_mix_prop    = 0.005
params.max_merge       = 50
params.multiqc         = 500
params.args_json       = null // null or '{\\"min_median_depth\\":5}'
params.clean_id        = params.id.replaceAll('[\\s-]+', '-')
params.publish         = true
params.publish_bams    = false
params.publish_gds     = true
params.no_report       = false // disable TBtypeR report generation
params.outdir          = 'output'
params.fast            = false

/* ---------- Experimental parameters ---------- */
params.dr_panel        = null // null or "$projectDir/resources/who_mut_cat_2021.tsv.bz2"

/* --------- Benchmarking parameters --------- */
params.benchmark = false
params.methods = ['TBtypeR', 'MixInfect', 'SplitStrains', 'quanttb', 'TBProfiler', 'fastlin']
params.mi_calling_regions = "$projectDir/workflows/benchmarking/MixInfect/misc/calling_regions.bed.gz"
params.fastlin_use_tbt_bc = false
params.tbprofiler_use_tbt_bc = false
params.quanttb_use_tbt_bc = false


include { TBtypeNF      } from './workflows/TBtype'
include { FastTBtypeNF  } from './workflows/FastTBtype'
include { Benchmark    } from './workflows/Benchmark'

workflow {
    if (params.benchmark) {
        Benchmark()
    } else {
        if (params.fast) {
            log.info("""\
                FastTBtypeNF
                ===================================
                id              : ${params.id}
                manifest        : ${params.manifest}
                ref_fasta       : ${params.ref_fasta}
                panel           : ${params.panel}
                """
                .stripIndent())
            FastTBtypeNF([:])
        } else {
            log.info("""\
                TBtypeNF
                ===================================
                id              : ${params.id}
                manifest        : ${params.manifest}
                ref_fasta       : ${params.ref_fasta}
                panel           : ${params.panel}
                dr_panel        : ${params.dr_panel}
                calling_regions : ${params.calling_regions}
                """
                .stripIndent())
            TBtypeNF([:])
        }
    }
}

workflow.onComplete { 
    log.info( 
        workflow.success ? 
        "Pipeline completed sucessfully." : 
        "Pipeline failed.")
}
