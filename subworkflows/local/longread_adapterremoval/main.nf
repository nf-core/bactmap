//
// Process long raw reads with porechop or porechop_abi
//

include { PORECHOP_PORECHOP } from '../../../modules/nf-core/porechop/porechop/main'
include { PORECHOP_ABI      } from '../../../modules/nf-core/porechop/abi/main'

workflow LONGREAD_ADAPTERREMOVAL {
    take:
    ch_reads
    ch_custom_adapters

    main:
    ch_multiqc_files = channel.empty()

    if (params.longread_adapterremoval_tool == 'porechop_abi') {
        PORECHOP_ABI(ch_reads, ch_custom_adapters)
        ch_processed_reads = PORECHOP_ABI.out.reads.map { meta, chopped_reads -> [meta + [single_end: true], chopped_reads] }
        ch_multiqc_files = ch_multiqc_files.mix(PORECHOP_ABI.out.log)
    }
    else if (params.longread_adapterremoval_tool == 'porechop') {
        PORECHOP_PORECHOP(ch_reads)
        ch_processed_reads = PORECHOP_PORECHOP.out.reads.map { meta, chopped_reads -> [meta + [single_end: true], chopped_reads] }
        ch_multiqc_files = ch_multiqc_files.mix(PORECHOP_PORECHOP.out.log)
    }
    else {
        ch_processed_reads = ch_reads
    }

    emit:
    reads    = ch_processed_reads // channel: [ val(meta), [ reads ] ]
    mqc      = ch_multiqc_files
}
