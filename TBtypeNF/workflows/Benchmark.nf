import static Helpers.path

// functions
include { standardise  } from '../functions/functions'
// processes
include { INDEX_REF    } from '../modules/index_ref'
// subworkflows
include { TBtypeNF     } from './TBtype'
include { FastTBtypeNF } from './FastTBtype'
include { SplitStrains } from './benchmarking/SplitStrains/main'
include { QuantTB      } from './benchmarking/QuantTB/main'
include { MixInfect    } from './benchmarking/MixInfect/main'
include { Fastlin      } from './benchmarking/Fastlin/main'
include { TBProfiler   } from './benchmarking/TBProfiler/main'

workflow Benchmark {

    results = Channel.of()
    ref = INDEX_REF(path(params.ref_fasta))
    
    if (params.methods.contains('TBtypeR')) {
        TBtypeNF(ref: ref)
        results = results.mix(standardise('TBtypeR', TBtypeNF.out.calls))
        TBtypeNF.out.calls |
            map { it.text } |
            collectFile(name: 'TBtypeR.results.csv',
                        storeDir: params.outdir,
                        keepHeader:true)
    }
    if (params.methods.contains('FastTBtypeR')) {
        results = results.mix(standardise('FastTBtypeR', FastTBtypeNF(ref: ref)))
    }
    if (params.methods.contains('MixInfect')) {
        MixInfect(ref: ref)
        results = results.mix(standardise('MixInfect', MixInfect.out))
    }
    if (params.methods.contains('SplitStrains')) {
        SplitStrains(ref: ref)
        results = results.mix(standardise('SplitStrains', SplitStrains.out))
    }
    if (params.methods.contains('quanttb')) {
        results = results.mix(standardise('quanttb', QuantTB(ref: ref)))
    }
    if (params.methods.contains('fastlin')) {
        results = results.mix(standardise('fastlin', Fastlin(ref: ref)))
    }
    if (params.methods.contains('TBProfiler')) {
        results = results.mix(standardise('TBProfiler', TBProfiler()))
    }

    results_file =
        results |
        map { it.values().join('\t') } |
        collectFile(name: 'combined_results.tsv',
                    newLine:true, sort:true,
                    seed: 'method\tsample\tstrain\tproportion',
                    storeDir: params.outdir)

    results_file |
        map { "results saved to: $it" } |
        view
    
}