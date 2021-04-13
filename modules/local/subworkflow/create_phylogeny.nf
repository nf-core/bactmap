/*
 * Phylogenies subworkflow
 */
params.rapidnj_options  = [:]
params.fasttree_options = [:]
params.iqtree_options   = [:]

include { RAPIDNJ }  from '../../nf-core/software/rapidnj/main'  addParams( options: params.rapidnj_options )
include { FASTTREE } from '../../nf-core/software/fasttree/main' addParams( options: params.fasttree_options )
include { IQTREE  }  from '../../nf-core/software/iqtree/main'   addParams( options: params.iqtree_options )


workflow CREATE_PHYLOGENY {
    take:
    fasta                 // channel: aligned pseudogenomes or filtered version
    constant_sites_string // val: string of constant sites A,C,G,T 
    
    main:
    /*
    * MODULE Fasttree
    */
    if (params.rapidnj_options.build){
        RAPIDNJ(fasta)
    }
    /*
     * MODULE Fasttree
     */
    if (params.fasttree_options.build){
        FASTTREE(fasta)
    }
    /*
     * MODULE Iqtree
     */
    if (params.iqtree_options.build){
        IQTREE(fasta, constant_sites_string)
    }

    emit:
    rapidnj_tree      = RAPIDNJ.out.phylogeny  // channel: [ phylogeny ]
    fasttree_tree     = FASTTREE.out.phylogeny // channel: [ phylogeny ]
    iqtree_tree       = IQTREE.out.phylogeny   // channel: [ phylogeny ]
    rapidnj_version   = RAPIDNJ.out.version    // path: *.version.txt
    fasttree_version  = FASTTREE.out.version   // path: *.version.txt
    iqtree_version    = IQTREE.out.version     // path: *.version.txt

}
