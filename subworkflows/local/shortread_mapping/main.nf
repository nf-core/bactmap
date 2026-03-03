//
// Perform short read mapping and variant calling
//

include { FASTQ_ALIGN_BWAMEM2                         } from '../fastq_align_bwamem2/main'
include { FASTQ_ALIGN_BOWTIE2                         } from '../../nf-core/fastq_align_bowtie2/main'
include { BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS } from '../../local/bam_variant_calling_sort_freebayes_bcftools/main'
include { BCFTOOLS_FILTER                             } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_NORM                               } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_STATS                              } from '../../../modules/nf-core/bcftools/stats/main'
include { CONSENSUS_BCFTOOLS                          } from '../consensus_bcftools/main'
include { SEQTK_COMP                                  } from '../../../modules/nf-core/seqtk/comp/main.nf'

workflow SHORTREAD_MAPPING {

    take:
    reads    // channel: [ val(meta), [ reads ] ]
    ch_fasta // channel: [meta, ref]
    ch_index // channel: [meta, ref index]
    ch_faidx // channel: [meta, ref fai]

    main:
    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    if (params.shortread_mapping_tool == 'bowtie2') {
        FASTQ_ALIGN_BOWTIE2 (
            reads,
            ch_index,
            false,
            false,
            ch_fasta
        )
        ch_bam           = FASTQ_ALIGN_BOWTIE2.out.bam
        ch_bam_index     = FASTQ_ALIGN_BOWTIE2.out.bai
        ch_multiqc_files = ch_multiqc_files.mix( FASTQ_ALIGN_BOWTIE2.out.stats )
        ch_versions      = ch_versions.mix( FASTQ_ALIGN_BOWTIE2.out.versions )
    } else {
        FASTQ_ALIGN_BWAMEM2 (
            reads,
            ch_index,
            ch_fasta,
            false
        )
        ch_bam           = FASTQ_ALIGN_BWAMEM2.out.bam
        ch_bam_index     = FASTQ_ALIGN_BWAMEM2.out.bai
        ch_multiqc_files = ch_multiqc_files.mix( FASTQ_ALIGN_BWAMEM2.out.stats )
        ch_versions      = ch_versions.mix( FASTQ_ALIGN_BWAMEM2.out.versions )
    }

    // Prepare inputs for FreeBayes
    ch_freebayes_fasta = ch_fasta // channel: [ val(meta), path(reference), path(fai)]
        .join( ch_faidx )

    freebayes_input = ch_bam  // channel: [ val(meta), path(bam) ]
        .join( ch_bam_index ) // channel: [ val(meta), path(bam), path(bam_index)]
            .map{
                meta, bam, bai -> [ meta, bam, bai, [], [], [] ]
            }

    BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS (freebayes_input,
                        ch_freebayes_fasta.first(),
                        [ [:], [] ],
                        [ [:], [] ],
                        [ [:], [] ]
    )
    ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS.out.versions)

    ch_bcftool_filter_input = BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS.out.vcf
        .join(BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS.out.tbi)

    BCFTOOLS_FILTER ( ch_bcftool_filter_input )

    ch_bcftool_norm_input = BCFTOOLS_FILTER.out.vcf.join(BCFTOOLS_FILTER.out.tbi)
    BCFTOOLS_NORM ( ch_bcftool_norm_input, ch_fasta )

    ch_bcftool_stats_input = BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.tbi)

    BCFTOOLS_STATS ( ch_bcftool_stats_input, [ [:], [] ], [ [:], [] ], [ [:], [] ], [ [:], [] ], [ [:], [] ] )
    ch_multiqc_files = ch_multiqc_files.mix( BCFTOOLS_STATS.out.stats )

    CONSENSUS_BCFTOOLS ( ch_bam, BCFTOOLS_NORM.out.vcf, BCFTOOLS_NORM.out.tbi, ch_fasta )
    ch_versions = ch_versions.mix( CONSENSUS_BCFTOOLS.out.versions )

    SEQTK_COMP( CONSENSUS_BCFTOOLS.out.consensus )

    emit:
    bam         = ch_bam                           // channel: [ val(meta), [ bam ] ]
    bai         = ch_index                         // channel: [ val(meta), [ bai ] ]
    vcf         = BCFTOOLS_NORM.out.vcf            // channel: [meta, vcf]
    csi         = BCFTOOLS_NORM.out.csi            // channel: [ val(meta), path(csi) ]
    tbi         = BCFTOOLS_NORM.out.tbi            // channel; [meta, tbi]
    stats       = BCFTOOLS_STATS.out.stats         // channel: [meta, stats]
    consensus   = CONSENSUS_BCFTOOLS.out.consensus // channel: [ val(meta), path(consensus) ]
    seqtk_stats = SEQTK_COMP.out.seqtk_stats       // channel: [meta, stats]
    versions    = ch_versions                      // channel: [ versions.yml ]
    mqc         = ch_multiqc_files                 // channel: [ val(meta), [ multiqc files ] ]
}
