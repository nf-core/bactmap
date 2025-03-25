//
// Perform read mapping and variant calling
//

include { BWAMEM2_MEM                                 } from '../../../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_FAIDX                              } from '../../../../modules/nf-core/samtools/faidx/main'

include { FASTQ_ALIGN_BOWTIE2                         } from '../../nf-core/fastq_align_bowtie2/main'
include { BAM_SORT_STATS_SAMTOOLS                     } from '../../nf-core/bam_sort_stats_samtools/main'
include { BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS } from '../../nf-core/bam_variant_calling_sort_freebayes/main'

workflow SHORTREAD_MAPPING {

    take:
    reads    // channel: [ val(meta), [ reads ] ]
    ch_fasta // channel: /path/to/reference.fasta
    ch_index // channel: /path/to/index/
    
    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    if (params.shortread_mapping_tool == 'bowtie2') {
        ch_unsorted_bam       = FASTQ_ALIGN_BOWTIE2 ( reads, bowtie_index, false, false, ch_fasta ).bam
        ch_unsorted_bam_index = FASTQ_ALIGN_BOWTIE2.out.bai
        ch_multiqc_files      = ch_multiqc_files.mix( FASTQ_ALIGN_BOWTIE2.out.stats )
        ch_versions           = ch_versions.mix( FASTQ_ALIGN_BOWTIE2.out.versions )
    } else {
        ch_unsorted_bam = BWAMEM2_MEM ( reads, ch_index, ch_fasta, false ).bam
        ch_versions     = ch_versions.mix( BWAMEM2_MEM.out.versions )
    }
    
    if (params.shortread_mapping_tool == 'bowtie2') {
        ch_processed_bam       = ch_unsorted_bam
        ch_processed_bam_index = ch_unsorted_bam_index
    } else {
        ch_processed_bam       = BAM_SORT_STATS_SAMTOOLS ( ch_unsorted_bam, ch_fasta ).bam
        ch_processed_bam_index = BAM_SORT_STATS_SAMTOOLS.out.bai
        ch_multiqc_files       = ch_multiqc_files.mix( BAM_SORT_STATS_SAMTOOLS.out.stats )
        ch_versions            = ch_versions.mix( BAM_SORT_STATS_SAMTOOLS.out.versions )
    }
    
    // Prepare inputs for FreeBayes
    
    ch_fasta                            // channel: [ val(meta), path(reference), path(fai)]
        .join(SAMTOOLS_FAIDX.fai)
        .set { ch_freebayes_fasta }

    freebayes_input = ch_processed_bam // channel: [ val(meta), path(bam) ]
        .join(ch_processed_bam_index)  // channel: [ val(meta), path(bam), path(bam_index)]
            .multiMap{
                meta, bam, bai ->
                   reads: [ meta, bam, bai, [], [], [] ]
            }

    BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS (freebayes_input.reads,
                        ch_freebayes_fasta,
                        [],
                        [],
                        []
    )
    ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS.out.versions)

    emit:
    bam      = ch_processed_bam                                    // channel: [ val(meta), [ bam ] ]
    bai      = ch_processed_bam_index                              // channel: [ val(meta), [ bai ] ]
    vcf      = BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS.out.vcf // channel: [ val(meta), path(vcf) ]
    csi      = BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS.out.csi // channel: [ val(meta), path(csi) ]
    tbi      = BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS.out.tbi // channel: [ val(meta), path(tbi) ]
    versions = ch_versions                                         // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}