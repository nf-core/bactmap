/*
 * Variant calling and downstream processing for BCFTools
 */

params.bcftools_mpileup_options    = [:]
params.bcftools_filter_options    = [:]

include { BCFTOOLS_MPILEUP } from '../modules/nf-core/software/bcftools/mpileup/main' addParams( options: params.bcftools_mpileup_options )
include { BCFTOOLS_FILTER  } from  '../modules/nf-core/software/bcftools/filter/main' addParams( options: params.bcftools_filter_options )


workflow VARIANTS_BCFTOOLS {
    take:
    bam       // channel: [ val(meta), [ bam ] ]
    fasta     // channel: /path/to/genome.fasta
    
    main:
    /*
     * MODULE Call variants
     */
    BCFTOOLS_MPILEUP ( bam, fasta )
    /*
     * MODULE Filter variants
     */
    BCFTOOLS_FILTER ( BCFTOOLS_MPILEUP.out.vcf  )

    emit:
    filtered_vcf     = BCFTOOLS_FILTER.out.vcf      // channel: [ val(meta), [ vcf ] ]
    vcf              = BCFTOOLS_MPILEUP.out.vcf     // channel: [ val(meta), [ vcf ] ]
    tbi              = BCFTOOLS_MPILEUP.out.tbi     // channel: [ val(meta), [ tbi ] ]
    stats            = BCFTOOLS_MPILEUP.out.stats   // channel: [ val(meta), [ txt ] ]
    bcftools_version = BCFTOOLS_MPILEUP.out.version //    path: *.version.txt

}
