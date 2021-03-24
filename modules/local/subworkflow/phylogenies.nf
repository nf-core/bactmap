/*
 * Phylogenies subworkflow
 */

params.fasttree_options    = [:]
//params.iqtree_options    = [:]

include { FASTTREE } from '../fasttree' addParams( options: params.fasttree_options )
//include { IQTREE  } from  '../iqtree' addParams( options: params.iqtree_options )


workflow PHYLOGENIES {
    take:
    fasta     // channel: aligned pseudogenomes or filtered version
    
    main:
    /*
     * MODULE Fasttree
     */
    FASTTREE ( fasta )
    /*
     * MODULE Iqtree
     */
    //IQTREE ( fasta  )

    emit:
    fasttree_tree     = FASTTREE.out.phylogeny      // channel: [ tree ]
    //iqtree_tree       = IQTREE.out.phylogeny    // channel: [ tree ]
    fasttree_version  = FASTTREE.out.version //    path: *.version.txt
    //iqtree_version    = IQTREE.out.version //    path: *.version.txt

}
