/*
 * Sub-sampling subworkflow
 */
params.mash_sketch_options    = [:]
params.rasusa_options    = [:]

include { MASH_SKETCH } from '../modules/nf-core/software/mash/sketch/main' addParams( mash_sketch: params.mash_sketch_options )
include { RASUSA } from '../modules/nf-core/software/rasusa/main'           addParams( options: params.rasusa_options)
include { find_genome_size } from '../modules/local/functions'



workflow SUB_SAMPLING {
    take:
    reads     // channel: INPUT_CHECK or FASTP
    
    main:
    
    //If genome size is not defined
    if (params.genome_size) {
        reads_and_genome_size = reads.combine([params.genome_size])
    } else {
        MASH_SKETCH (
            reads
        )
        genome_size = MASH_SKETCH.out.stats.map { meta, file -> [meta, find_genome_size(file.text)]}
        reads_and_genome_size = reads.combine(genome_size, by: 0)
        
    }
    RASUSA(reads_and_genome_size, params.subsampling_depth_cutoff)

    emit:
    reads     = RASUSA.out.reads      // channel: [ reads ]
    version   = params.genome_size ? RASUSA.out.version : MASH_SKETCH.out.version.mix(RASUSA.out.version.first().ifEmpty(null)) //    path: *.version.txt
}
