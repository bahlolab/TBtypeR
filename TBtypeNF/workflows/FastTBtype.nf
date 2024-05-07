import static Helpers.path
// FUNCTIONS
include { read_manifest  } from '../functions/functions'
// SUBWORKFLOWS
include { DOWNLOAD } from '../subworkflows/download'
// PROCESSES
include { INDEX_REF       } from '../modules/index_ref'
include { MAKE_FASTLIN_BC } from '../modules/make_fastlin_bc'
include { FASTLIN         } from '../modules/fastlin'
include { TBTYPER         } from '../modules/tbtyper_fast'

workflow FastTBtypeNF {
    take:
    opts

    main:
    manifest =  read_manifest(path(params.manifest))
    types = manifest.collect { it.type }.unique()
    if (types.contains('bam')) { error ("Bams not supported in fast mode") }
    n_sample = manifest.findAll{ it.type != 'error' }.collect{it.sample}.unique().size()
    
    manifest_ch = Channel.fromPath(params.manifest)

    ref_ch = opts.ref ?: INDEX_REF(path(params.ref_fasta))
    panel_ch = Channel.value(path(params.panel))

    inputs =  DOWNLOAD(manifest) \
        | branch { 
            paired: it.type == 'paired'
            single: it.type == 'single' 
          }
  
    bc_ch = MAKE_FASTLIN_BC(ref_ch, panel_ch)

    fq_ch = inputs.paired \
        | map { [it.sample, it.fastq1, it.fastq2 ] } \
        | mix(inputs.single.map {[it.sample, it.fastq1, path("$projectDir/resources/DEFAULT")] })

    FASTLIN(bc_ch, fq_ch)

    results = TBTYPER(panel_ch, params.args_json ?: false, FASTLIN.out) \
        | map { it[1] } \
        | collectFile(
            name: 'FastTBtypeR.results.csv',
            storeDir: params.outdir,
            keepHeader:true)

    emit:
    results
}

