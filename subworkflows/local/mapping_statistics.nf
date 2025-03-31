//
// Summarize mappping statistics with SEQTK
//

include { SEQTK_COMP } from '../../modules/nf-core/seqtk/comp/main'
include { SEQTK_PARSE } from '../../modules/local/seqtk_parse/main'

workflow MAPPING_STATISTICS {
   take:
   fasta // [[meta], [fasta]]

   main:
   ch_versions = Channel.empty()

   SEQTK_COMP( fasta )
   ch_versions = ch_versions.mix( SEQTK_COMP.out.versions )

   ch_seqtk_stats  =  SEQTK_COMP.out.seqtk_stats
                       .map( it -> it[1])
                       .collect()

   SEQTK_PARSE( ch_seqtk_stats )
   ch_versions = ch_versions.mix( SEQTK_PARSE.out.versions )

   emit:
   ch_mapping_stats = SEQTK_PARSE.out.tsv // Channel: [ mapping_summary.tsv ]
   ch_versions      // Channel: [ version.yml ]

}

