// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { MINIMAP2_INDEX      } from '../../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN      } from '../../../modules/nf-core/minimap2/align/main'

workflow MINIMAP2_ALIGNMENT {

    take:

    ch_ref 
    ch_fasta

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    MINIMAP2_INDEX ( ch_ref )
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions.first())

    MINIMAP2_ALIGN ( MINIMAP2_INDEX.out,  )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    align      = MINIMAP2_ALIGN.out.bam           // channel: [ val(meta), [ bam ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

