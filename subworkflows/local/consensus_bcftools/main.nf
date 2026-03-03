//
// Consensus calling with BCFTools
//

include { BCFTOOLS_QUERY      } from '../../../modules/nf-core/bcftools/query/main'
include { BEDTOOLS_GENOMECOV  } from '../../../modules/local/bedtools/genomecov/main'
include { BEDTOOLS_SUBTRACT   } from '../../../modules/nf-core/bedtools/subtract/main'
include { BCFTOOLS_CONSENSUS  } from '../../../modules/local/bcftools/consensus/main'
include { RENAME_FASTA_HEADER } from '../../../modules/local/fasta_rename/main'
include { CONCATENATE_FASTA   } from '../../../modules/local/concatenate_fasta/main'

workflow CONSENSUS_BCFTOOLS {

    take:
    ch_bam          // channel: [ val(meta), [ bam ] ]
    ch_vcf          // channel: [ val(meta), [ vcf ] ]
    ch_tbi          // channel: [ val(meta), [ tbi ] ]
    ch_fasta        // channel: /path/to/genome.fasta

    main:
    ch_versions = channel.empty()

    //
    // Filter variants by allele frequency, zip and index
    //
    ch_query = ch_vcf.join(ch_tbi)
        .map{
            meta, vcf, tbi -> [ meta, vcf, tbi ]
        }

    BCFTOOLS_QUERY  (
        ch_query,
        [],
        [],
        []
    )

    ch_genomecov = ch_bam
        .multiMap {
            meta, bam ->
            genomecov: [meta, bam, params.genomecov_scale ]
        }

    BEDTOOLS_GENOMECOV (
        ch_genomecov.genomecov
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

    //
    // Make the bed mask
    //
    ch_bedtools_subtract = BEDTOOLS_GENOMECOV.out.genomecov.join(BCFTOOLS_QUERY.out.output)

    BEDTOOLS_SUBTRACT (
        ch_bedtools_subtract
        )

    ch_consensus = ch_vcf
        .join( ch_tbi)
        .join( BEDTOOLS_SUBTRACT.out.bed ).map{
            meta, vcf, tbi, bed -> [ meta, vcf, tbi, bed ]
        }

    BCFTOOLS_CONSENSUS (
        ch_consensus,
        ch_fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions.first())

    //
    // Rename consensus header adding sample name
    //
    RENAME_FASTA_HEADER (
        BCFTOOLS_CONSENSUS.out.fasta
    )
    ch_versions = ch_versions.mix(RENAME_FASTA_HEADER.out.versions.first())

    CONCATENATE_FASTA (
        RENAME_FASTA_HEADER.out.fasta
    )
    ch_versions = ch_versions.mix(RENAME_FASTA_HEADER.out.versions.first())


    emit:
    consensus = CONCATENATE_FASTA.out.fasta // channel: [ val(meta), [ fasta ] ]
    versions  = ch_versions                 // channel: [ versions.yml ]
}
