//
// Phylogenies subworkflow
//
params.rapidnj_options  = [:]
params.fasttree_options = [:]
params.iqtree_options   = [:]
params.raxmlng_options  = [:]

include { RAPIDNJ  } from '../modules/nf-core/software/rapidnj/main'  addParams( options: params.rapidnj_options  )
include { FASTTREE } from '../modules/nf-core/software/fasttree/main' addParams( options: params.fasttree_options )
include { IQTREE   } from '../modules/nf-core/software/iqtree/main'   addParams( options: params.iqtree_options   )
include { RAXMLNG  } from '../modules/nf-core/software/raxmlng/main'  addParams( options: params.raxmlng_options  )

workflow CREATE_PHYLOGENY {
    take:
    fasta                 // channel: aligned pseudogenomes or filtered version
    constant_sites_string // val: string of constant sites A,C,G,T

    main:
    //
    // MODULE: rapidnj
    //
    rapidnj_tree = Channel.empty()
    rapidnj_version = null
    if (params.rapidnj) {
        RAPIDNJ(fasta)
        rapidnj_tree = RAPIDNJ.out.phylogeny
        rapidnj_version = RAPIDNJ.out.version
    }

    //
    // MODULE: fasttree
    //
    fasttree_tree = Channel.empty()
    fasttree_version = null
    if (params.fasttree) {
        FASTTREE(fasta)
        fasttree_tree = FASTTREE.out.phylogeny
        fasttree_version = FASTTREE.out.version
    }

    //
    // MODULE: iqtree
    //
    iqtree_tree = Channel.empty()
    iqtree_version = null
    if (params.iqtree) {
        IQTREE(fasta, constant_sites_string)
        iqtree_tree = IQTREE.out.phylogeny
        iqtree_version = IQTREE.out.version
    }

    //
    // MODULE: raxmlng
    //
    raxmlng_tree = Channel.empty()
    raxmlng_version = null
    if (params.raxmlng) {
        RAXMLNG(fasta)
        raxmlng_tree = RAXMLNG.out.phylogeny
        raxmlng_version = RAXMLNG.out.version
    }

    emit:
    rapidnj_tree      = rapidnj_tree     // channel: [ phylogeny ]
    fasttree_tree     = fasttree_tree    // channel: [ phylogeny ]
    iqtree_tree       = iqtree_tree      // channel: [ phylogeny ]
    raxmlng_tree      = raxmlng_tree     // channel: [ phylogeny ]
    rapidnj_version   = rapidnj_version  // path: *.version.txt
    fasttree_version  = fasttree_version // path: *.version.txt
    iqtree_version    = iqtree_version   // path: *.version.txt
    raxmlng_version   = raxmlng_version  // path: *.version.txt
}
