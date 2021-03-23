/*
 * Variant calling and downstream processing for BCFTools
 */

params.bcftools_mpileup_options    = [:]

include { BCFTOOLS_MPILEUP } from '../../nf-core/software/bcftools/mpileup/main' addParams( options: params.bcftools_mpileup_options )

workflow VARIANTS_BCFTOOLS {
    take:
    bam       // channel: [ val(meta), [ bam ] ]
    fasta     // channel: /path/to/genome.fasta
    
    main:
    /*
     * Call variants
     */
    BCFTOOLS_MPILEUP ( bam, fasta )

    emit:
    vcf              = BCFTOOLS_MPILEUP.out.vcf     // channel: [ val(meta), [ vcf ] ]
    tbi              = BCFTOOLS_MPILEUP.out.tbi     // channel: [ val(meta), [ tbi ] ]
    stats            = BCFTOOLS_MPILEUP.out.stats   // channel: [ val(meta), [ txt ] ]
    bcftools_version = BCFTOOLS_MPILEUP.out.version //    path: *.version.txt

}
