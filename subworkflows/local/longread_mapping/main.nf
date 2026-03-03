//
// Perform long read mapping and variant calling
//

include { MINIMAP2_ALIGNMENT      } from '../minimap2_alignment/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../../nf-core/bam_sort_stats_samtools/main'
include { CLAIR3                  } from '../../../modules/local/clair3/main'
include { BCFTOOLS_SORT           } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_INDEX          } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_VIEW           } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_NORM           } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_STATS          } from '../../../modules/nf-core/bcftools/stats/main'
include { CONSENSUS_BCFTOOLS      } from '../consensus_bcftools/main'
include { SEQTK_COMP              } from '../../../modules/nf-core/seqtk/comp/main.nf'

workflow LONGREAD_MAPPING {

    take:
    ch_fasta // channel: [meta, ref]
    ch_faidx // channel: [meta, ref index]
    reads    // channel: [meta2, fasta/fastq]

    main:
    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    MINIMAP2_ALIGNMENT( ch_fasta, reads )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGNMENT.out.versions)

    BAM_SORT_STATS_SAMTOOLS ( MINIMAP2_ALIGNMENT.out.minimap_align,  ch_fasta )
    ch_multiqc_files = ch_multiqc_files.mix( BAM_SORT_STATS_SAMTOOLS.out.stats )
    ch_versions      = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    ch_clair3_input = BAM_SORT_STATS_SAMTOOLS.out.bam
        .join(BAM_SORT_STATS_SAMTOOLS.out.bai)
            .multiMap{
                    meta, bam, bai ->
                    bam: [ meta, bam, bai, params.clair3_model, params.clair3_platform ]
            }

    CLAIR3 (ch_clair3_input.bam, ch_fasta, ch_faidx)
    ch_versions = ch_versions.mix(CLAIR3.out.versions.first())

    BCFTOOLS_SORT ( CLAIR3.out.vcf )

    BCFTOOLS_INDEX ( BCFTOOLS_SORT.out.vcf )

    ch_bcftool_view_input = BCFTOOLS_SORT.out.vcf.join(BCFTOOLS_INDEX.out.tbi)
    BCFTOOLS_VIEW ( ch_bcftool_view_input, [], [], [] )

    ch_bcftool_norm_input = BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.tbi)
    BCFTOOLS_NORM ( ch_bcftool_norm_input, ch_fasta )

    ch_bcftool_stats_input = BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.tbi)

    BCFTOOLS_STATS ( ch_bcftool_stats_input, [ [:], [] ], [ [:], [] ], [ [:], [] ], [ [:], [] ], [ [:], [] ] )
    ch_multiqc_files = ch_multiqc_files.mix( BCFTOOLS_STATS.out.stats )

    CONSENSUS_BCFTOOLS ( BAM_SORT_STATS_SAMTOOLS.out.bam, BCFTOOLS_NORM.out.vcf, BCFTOOLS_NORM.out.tbi, ch_fasta )
    ch_versions = ch_versions.mix( CONSENSUS_BCFTOOLS.out.versions )

    SEQTK_COMP( CONSENSUS_BCFTOOLS.out.consensus )

    emit:
    bam         = BAM_SORT_STATS_SAMTOOLS.out.bam  // channel: [ val(meta), [ bam ] ]
    bai         = BAM_SORT_STATS_SAMTOOLS.out.bai  // channel: [ val(meta), [ bai ] ]
    vcf         = BCFTOOLS_NORM.out.vcf          // channel: [meta, vcf]
    csi         = BCFTOOLS_NORM.out.csi          // channel: [ val(meta), path(csi) ]
    tbi         = BCFTOOLS_NORM.out.tbi          // channel; [meta, tbi]
    stats       = BCFTOOLS_STATS.out.stats         // channel: [meta, stats]
    consensus   = CONSENSUS_BCFTOOLS.out.consensus // channel: [ val(meta), path(consensus) ]
    seqtk_stats = SEQTK_COMP.out.seqtk_stats       // channel: [meta, stats]
    versions    = ch_versions                      // channel: [ versions.yml ]
    mqc         = ch_multiqc_files                 // channel: [ val(meta), [ multiqc files ] ]
}

