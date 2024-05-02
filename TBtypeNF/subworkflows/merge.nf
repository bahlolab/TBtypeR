
include { BCFTOOLS_MERGE as PRE_MERGE; BCFTOOLS_MERGE } from '../modules/bcftools_merge'

workflow MERGE {
    take:
    vcfs // set, vcfs
    ref_ch // ref, ref_sites
    sites
    n_sample

    main:
    if (n_sample > params.max_merge) {
        vcfs = vcfs \
            | toSortedList { a, b -> a[0] <=> b[0] } \
            | flatMap { it.collect { it.drop(1) } } \
            | collate(params.max_merge) \
            | toList() \
            | flatMap { it.withIndex().collect { item, idx -> [idx + 1] + item.transpose() } } \
            | map { id, snps, indels ->
                 indels[0].name == 'NO_INDELS' ? [id, snps, indels[[0]]] : [id, snps, indels] 
              } \
            | combine(ref_ch) \
            | combine(sites) \
            | PRE_MERGE
    }

    merged = vcfs \
        | toSortedList { a, b -> a[0] <=> b[0] } \
        | map { ['all'] + it.collect { it.drop(1) }.transpose() } \
        | map { id, snps, indels ->
             indels[0].name == 'NO_INDELS' ? [id, snps, indels[[0]]] : [id, snps, indels] 
          } \
        | combine(ref_ch) \
        | combine(sites) \
        | BCFTOOLS_MERGE

    emit:
        snp   = merged.map { it[1] }
        indel = merged.map { it[2] }

}