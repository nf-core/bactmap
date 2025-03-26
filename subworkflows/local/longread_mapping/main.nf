// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules
// minimap2, bam_sort_stats_samtools, clair3, bctools sort, index, filter, stats
include { MINIMAP2_ALIGNMENT      } from '../minimap2_alignment/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../../nf-core/bam_sort_stats_samtools/main'
include { CLAIR3                  } from '../../../modules/nf-core/clair3/main'
include { BCFTOOLS_SORT           } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_INDEX          } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_FILTER         } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_STATS          } from '../../../modules/nf-core/bcftools/stats/main'

workflow LONGREAD_MAPPING {

    take:
    ch_ref      // channel: [meta, ref]
    ch_ref_idx  // channel: [meta, ref index]
    ch_fasta    // channel: [meta2, fasta/fastq]

    main:

    ch_versions = Channel.empty()

    MINIMAP2_ALIGNMENT( ch_ref, ch_fasta )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGNMENT.out.versions.first())

    BAM_SORT_STATS_SAMTOOLS ( MINIMAP2_ALIGNMENT.out.minimap_align,  ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions.first())
    
    ch_clair3_input = MINIMAP2_ALIGNMENT.out.minimap_align.join(MINIMAP2_ALIGNMENT.out.minimap_index).map {
        meta, bam, bai -> 
        [meta, bam, bai, file(params.clair3_model), params.clair3_platform]
    }
    CLAIR3 (ch_clair3_input, ch_ref, ch_ref_idx)
    ch_versions = ch_versions.mix(CLAIR3.out.versions.first())

    BCFTOOLS_SORT ( CLAIR3.out.vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())

    BCFTOOLS_INDEX ( BCFTOOLS_SORT.out.vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    ch_bcftool_filter_input = BCFTOOLS_SORT.out.vcf.join(BCFTOOLS_INDEX.out.tbi)
    BCFTOOLS_FILTER ( ch_bcftool_filter_input )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())
    
    ch_bcftool_stat_input = BCFTOOLS_FILTER.out.vcf.join(BCFTOOLS_FILTER.out.tbi)

    // create channelS with additional input for BCFTOOLS_STAT

    ch_regions = params.bcftools_regions ? Channel.fromPath( params.bcftools_regions, checkIfExists: true ).map(path -> return([description:"regions"], path)) : Channel.of([description:"regions"], []) 
    ch_targets = params.bcftools_targets ? Channel.fromPath( params.bcftools_targets, checkIfExists: true ).map(path -> return([description:"targets"], path)) : Channel.of([description:"targets"], [])
    ch_samples = params.bcftools_samples ? Channel.fromPath( params.bcftools_samples, checkIfExists: true ).map(path -> return([description:"samples"], path)) : Channel.of([description:"samples"], [])
    ch_exons   = params.bcftools_exons ? Channel.fromPath( params.bcftools_exons, checkIfExists: true ).map(path -> return([description:"exons"], path)) : Channel.of([description:"exons"], [])

    BCFTOOLS_STATS ( ch_vcf_index, ch_regions, ch_targets, ch_samples, ch_exons, ch_fasta)
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())
    

    emit:
    stats = BCFTOOLS_STATS.out.stats // channel: [meta, stats]
    vcf   = BCFTOOLS_FILTER.out.vcf  // channel: [meta, vcf]
    tbi   = BCFTOOLS_FILTER.out.tbi  // channel; [meta, tbi]

    versions = ch_versions           // channel: [ versions.yml ]
}

