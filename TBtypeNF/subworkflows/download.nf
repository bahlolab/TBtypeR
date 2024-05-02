
import static Helpers.path

include { WGET } from '../modules/wget'

workflow DOWNLOAD {
    take: manifest

    main:
    any_download = manifest.any { it.file ==~ /^https?:.+/ }

    outputs = Channel.fromList(manifest) \
        | filter { it.type != 'error '} \
        | filter { !(it.file ==~ /^https?:.+/) } \
        | map { it + [file: path(it.file)] }

    if (any_download) {
        outputs = Channel.fromList(manifest) \
            | filter { it.type != 'error '} \
            | filter { it.file ==~ /^https?:.+/ } \
            | map { [[sample: it.sample, type: it.type, read: it.read], it.file] } \
            | WGET \
            | map { meta, file -> meta + [file: file] } \
            | mix(outputs)
    } 

    output_types = outputs \
        | map { [ it.type, it.read, it.sample, it.file] } \
        | branch {
            single: it[0] == 'single'
            paired: it[0] == 'paired'
            bam:    it[0] == 'bam'
        }
    
    single = output_types.single \
        | map { tp, rd, sm, fl -> [sample: sm, type: 'single', fastq1: fl] } 
    
    paired = output_types.paired \
        | filter { it[1] == 'R1'} \
        | map { it[2..3] } \
        | join(
            output_types.paired \
                | filter { it[1] == 'R2'} \
                | map { it[2..3] },
            by:0, remainder: false 
            ) \
        | map { sm, f1, f2 ->  [sample: sm, type: 'paired', fastq1: f1, fastq2: f2] }

    bam = output_types.bam \
        | filter { it[1] == 'bam'} \
        | map { it[2..3] } \
        | join(
            output_types.bam \
                | filter { it[1] == 'bai'} \
                | map { it[2..3] },
            by:0, remainder: false 
            ) \
        | map { sm, bam, bai ->  [sample: sm, type: 'bam', bam: bam, bai: bai] }

    outputs = single.mix(paired).mix(bam)

    emit: outputs
}
