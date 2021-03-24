/*
 * Phylogenies subworkflow
 */

params.fasttree_options    = [:]
//params.iqtree_options    = [:]

include { MASH_SKETCH } from '../mash_sketch' addParams( mash_sketch: params.mash_sketch_options )
include { RASUSA } from '../rasusa' addParams( options: params.rasusa_options)
include { find_genome_size } from '../functions'



workflow SUB_SAMPLING {
    take:
    reads     // channel: INPUT_CHECK or FASTP
    
    main:
    
    //If genome size is not defined
    if (!params.genome_size) {
         MASH_SKETCH (
            reads
        )

        ch_genome_size = MASH_SKETCH.out.stats.map { meta, file -> find_genome_size(meta.id, file.text)}
    
        // INPUT_CHECK.out.sample_info.view()
    
        // ch_genome_size.view()
        //TO-DO: May need to find neater way 
        joined = reads
                      .map{it -> tuple(it[0].id,it[0],it[1])}
                      .join(ch_genome_size)
                      .map{it -> tuple(it[1],it[2],it[3])}
        // joined.view()
    
    } else {
        //Create a value channel for genome size
        ch_genome_size = Channel.from(params.genome_size)
        
        joined = MASH_SKETCH.out.stats.combine(ch_genome_size)
    }
        //Sub sampling give a depth cutoff and genome size
        RASUSA(joined)

    emit:
    reads     = RASUSA.out.reads      // channel: [ reads ]
    version   = MASH_SKETCH.out.version.mix(RASUSA.out.version.first().ifEmpty(null)) //    path: *.version.txt
}
