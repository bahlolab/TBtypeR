
import static Helpers.read_tsv
import static Helpers.path

def read_manifest(Path file) {

    ArrayList manifest  = read_tsv(file, required: ['sample'])
        .collectMany {
            if (it.fastq1 != null & it.fastq2 != null) {
                [
                    [type: 'paired', sample: it.sample, read: 'R1', file: it.fastq1],
                    [type: 'paired', sample: it.sample, read: 'R2', file: it.fastq2]
                ]
            } else if (it.fastq1) {
                [
                    [type: 'single', sample: it.sample, read: 'R1', file: it.fastq1]
                ]
            } else if (it.bam) {
                [
                    [type: 'bam', sample: it.sample, read: 'bam', file: it.bam],
                    [type: 'bam', sample: it.sample, read: 'bai', file: "${it.bam}.bai" ]
                ]
            } else if (it.cram) {
                [
                    [type: 'bam', sample: it.sample, read: 'bam', file: it.cram],
                    [type: 'bam', sample: it.sample, read: 'bai', file: "${it.cram}.crai" ]
                ]
            } else {
                [
                    [type: 'error', sample: it.sample]
                ]
            }
        }
    return manifest
}

def standardise(source, data) {

    // sample, strain, proportion
    if (source == 'TBtypeR' | source == 'FastTBtypeR') {
        result =
            data |
                splitCsv(header:true) |
                map { [method: source,
                       sample: it.sample_id,
                       strain: it.mix_phylotype ?: 'NA',
                       proportion: it.mix_prop ?: 'NA'] }
    } else if (source == 'MixInfect') {
        result =
            data |
                splitCsv(header: true) |
                flatMap {
                    res = [method: source, sample: it['Sample name']]
                    num_mix = it['Number of strains in mix'] as int
                    if (num_mix == 1) {
                        [res + [strain: 'NA', proportion: 1]]
                    } else {
                        maj_prop = it['Major strain proportion'] as double
                        (1..num_mix).collect {
                            res + [strain: 'NA', proportion: it == 1 ? maj_prop : it == 2 ? 1 - maj_prop : 'NA']
                        }
                    }
                }
    } else if (source == 'quanttb') {
        result =
            data |
                splitCsv(sep: '\t', header:true) |
                map { [method: source,
                       sample: it.sample,
                       strain: it.refname,
                       proportion: it.relabundance ] }
    } else if (source == 'SplitStrains') {
        result = data |
            splitCsv(sep: '\t', header:true) |
            flatMap {
                it.proportions.split(' ').collect { pr ->
                    [method: source,
                     sample: it.file.replaceAll('\\.bam$', ''),
                     strain: 'NA',
                     proportion: pr ] }
            }
    } else if (source == 'fastlin') {
        result = data |
            splitCsv(sep: '\t', header:true) |
            flatMap {
                if (it.lineages == "") {
                    // no kmer coverage - couldn't assign lineage
                    return [[
                        method: source,
                        sample: it.sample,
                        strain: 'NA',
                        proportion: 'NA'
                    ]]
                }
                lins = it.lineages.split(',').collect { (it =~ /(.+) \((.+)\)/)[0][1] }
                counts = it.lineages.split(',').collect { (it =~ /(.+) \((.+)\)/)[0][2] as int}
                props = counts.collect { it / counts.sum() }
                lins.withIndex().collect { lin, i ->
                    [method: source,
                     sample: it.sample,
                     strain: lin,
                     proportion: props[i] ]
                }
            }
    } else if (source == 'TBProfiler') {
        result = data |
            splitCsv(sep: '\t', header:true) |
            map { [
                method: source,
                sample: it.sample,
                strain: it.strain -~ /^lineage/,
                proportion: it.proportion 
            ] }
    }
    return result
}

