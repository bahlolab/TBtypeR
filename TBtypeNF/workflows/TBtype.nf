import static Helpers.path
// FUNCTIONS
include { read_manifest  } from '../functions/functions'
// SUBWORKFLOWS
include { DOWNLOAD } from '../subworkflows/download'
include { MERGE    } from '../subworkflows/merge'
// PROCESSES
include { INDEX_REF       } from '../modules/index_ref'
include { GET_TARGETS     } from '../modules/get_targets'
include { FASTP_PE        } from '../modules/fastp'
include { FASTP_SE        } from '../modules/fastp'
include { ALIGN_PE        } from '../modules/align'
include { ALIGN_SE        } from '../modules/align'
include { SAMTOOLS_STATS  } from '../modules/samtools_stats'
include { BCFTOOLS_CALL   } from '../modules/bcftools_call'
include { COMBINE_SITES   } from '../modules/combine_sites'
include { BCFTOOLS_CONCAT } from '../modules/bcftools_concat'
include { FILTER_SNPS     } from '../modules/filter_snps'
include { VCF_TO_GDS      } from '../modules/vcf2gds'
include { TBTYPER         } from '../modules/tbtyper'
include { MOSDEPTH        } from '../modules/mosdepth'
include { MULTIQC         } from '../modules/multiqc'


workflow TBtypeNF {
    take:
    opts

    main:
    manifest =  read_manifest(path(params.manifest))
    types = manifest.collect { it.type }.unique()
    n_sample = manifest.findAll{ it.type != 'error' }.collect{it.sample}.unique().size()
    
    manifest_ch = Channel.fromPath(params.manifest)

    ref_ch = opts.ref ?: INDEX_REF(path(params.ref_fasta))
    regions_ch = Channel.value(
        params.calling_regions ? path(params.calling_regions) : path("$projectDir/resources/NO_REGIONS")
    )
    panel_ch = Channel.value(path(params.panel))
    dr_panel_ch = Channel.value(
        params.dr_panel ? path(params.dr_panel) : path("$projectDir/resources/NO_DR_PANEL")
    )

    targets = GET_TARGETS(ref_ch,
                          regions_ch,
                          panel_ch,
                          dr_panel_ch) 
    targets_indel = targets.indel.ifEmpty(path("$projectDir/resources/NO_INDELS"))

    inputs =  DOWNLOAD(manifest) \
        | branch { 
            paired: it.type == 'paired'
            single: it.type == 'single' 
            bam:    it.type == 'bam' 
          }

    fastp_json = Channel.fromList([])

    if (types.contains('bam')) {
        bams = inputs.bam \
            | map { [it.sample, it.bam, it.bai ] } \
    } else {
        bams = Channel.empty()
    }

    if (types.contains('single')) {
        bams = inputs.single \
            | map { [it.sample, it.fastq1 ] } \
            | FASTP_SE \
            | map { it.take(2) } \
            | combine(ref_ch) \
            | ALIGN_SE
        fastp_json = fastp_json.mix(FASTP_SE.out.map {it[[0,3]]})
    }

    if (types.contains('paired')) {
        bams = inputs.paired \
            | map { [it.sample, it.fastq1, it.fastq2 ] } \
            | FASTP_PE \
            | map { it.take(3) } \
            | combine(ref_ch) \
            | ALIGN_PE \
            | mix(bams)
        fastp_json = fastp_json.mix(FASTP_PE.out.map {it[[0,3]]})
    }

    calls = bams \
        | combine(ref_ch) \
        | combine(targets.snp) \
        | combine(targets_indel) \
        | combine(targets.panel) \
        | BCFTOOLS_CALL

    sites = calls.sites \
        | mix(targets.panel) \
        | toSortedList { a, b -> a.name <=> b.name } \
        | COMBINE_SITES


    merged = MERGE(calls.calls, ref_ch, sites, n_sample)
    
    merged.snp \
        | combine(targets.panel) \
        | first \
        | FILTER_SNPS \
        | combine(merged.indel) \
        | first \
        | BCFTOOLS_CONCAT \
        | VCF_TO_GDS \
        | combine(panel_ch) \
        | combine(dr_panel_ch) \
        | combine([params.args_json]) \
        | combine(manifest_ch) \
        | first() \
        | TBTYPER

    if (params.multiqc) {
        bams \
            | combine(regions_ch) \
            | MOSDEPTH \
            | combine(SAMTOOLS_STATS(bams), by:0) \
            | join(fastp_json, by:0, remainder:true) \
            | toSortedList { a, b -> a[0] <=> b[0] } \
            | flatMap { it.collect { it.drop(1).flatten() }} \
            | collate(params.multiqc) \
            | toList() \
            | flatMap { 
                it.withIndex().collect { item, idx -> 
                    [idx + 1, item.flatten().findAll { it != null } ] 
                } } \
            | MULTIQC \
            | map { it[0] } \
            | collectFile(name: "${params.clean_id}.multiqc_general_stats.tsv",
                        storeDir: 'output', keepHeader: true)
    }

    emit:
    calls = TBTYPER.out[0].map { it[1] }
    bams = bams
}