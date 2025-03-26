//
// Consensus calling with BCFTools
//

include { BCFTOOLS_QUERY            } from '../../modules/nf-core/bcftools/filter/main'
include { BEDTOOLS_GENOMECOV        } from '../../modules/local/bedtools/genomecov/main'
include { BEDTOOLS_SUBTRACT         } from '../../modules/nf-core/bedtools/subtract/main'
include { BCFTOOLS_CONSENSUS        } from '../../modules/nf-core/bcftools/consensus/main'
include { RENAME_FASTA_HEADER       } from '../../modules/local/fasta_rename/main'

workflow CONSENSUS_BCFTOOLS {
    take:
    bam          // channel: [ val(meta), [ bam ] ]
    vcf          // channel: [ val(meta), [ vcf ] ]
    tbi          // channel: [ val(meta), [ tbi ] ]
    fasta        // channel: /path/to/genome.fasta

    main:

    ch_versions = Channel.empty()

    //
    // Filter variants by allele frequency, zip and index
    //
    ch_query = vcf.join(tbi)

    BCFTOOLS_QUERY  (
        ch_query,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions.first())

    ch_genomecov = bam.map { meta, bam ->
        [meta, bam, [] ]
        }


    BEDTOOLS_GENOMECOV (
        ch_genomecov
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

    //
    // Make the bed mask
    //
    BEDTOOLS_SUBTRACT (
        BEDTOOLS_GENOMECOV.out.genomecov.join(BBCFTOOLS_QUERY.out.output)
        )
    ch_versions = ch_versions.mix(BEDTOOLS_SUBTRACT.out.versions.first())

    //
    // Build consensus
    //
    ch_consensus = vcf.join(tbi, fasta, BEDTOOLS_SUBTRACT.out.bed)

    BCFTOOLS_CONSENSUS (
        ch_consensus
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions.first())

    //
    // Rename consensus header adding sample name
    //
    RENAME_FASTA_HEADER (
        BCFTOOLS_CONSENSUS.out.fasta
    )
    ch_versions = ch_versions.mix(RENAME_FASTA_HEADER.out.versions.first())


    emit:
    consensus        = RENAME_FASTA_HEADER.out.fasta     // channel: [ val(meta), [ fasta ] ]
    versions         = ch_versions                       // channel: [ versions.yml ]
}
