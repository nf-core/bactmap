/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_bactmap_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

// include { FASTQSCANPARSE as FASTQSCANPARSE_TRIM      } from '../modules/local/fastqscanparse/main'
// include { FASTQSCANPARSE as FASTQSCANPARSE_SUBSAMPLE } from '../modules/local/fastqscanparse/main'
// include { SEQTK_COMP                                 } from '../modules/local/seqtk_comp/main'
// include { SEQTK_PARSE                                } from '../modules/local/seqtk_parse/main'
// include { ALIGNPSEUDOGENOMES                         } from '../modules/local/alignpseudogenomes/main'

// include { SHORTREAD_PREPROCESSING                    } from '../subworkflows/local/shortread_preprocessing'
// include { LONGREAD_PREPROCESSING                     } from '../subworkflows/local/longread_preprocessing'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// NF-CORE MODULES/PLUGINS
//
// include { BWAMEM2_INDEX                               } from '../modules/nf-core/bwamem2/index/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
// include { CAT_FASTQ as MERGE_RUNS                     } from '../modules/nf-core/cat/fastq/main'
// include { FASTQSCAN as FASTQSCAN_TRIM                 } from '../modules/nf-core/modules/fastqscan/main'
// include { RASUSA                                      } from '../modules/nf-core/rasusa/main'
// include { FASTQSCAN as FASTQSCAN_SUBSAMPLE            } from '../modules/nf-core/modules/fastqscan/main'
// include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_SHORT      } from '../modules/nf-core/minimap2/align/main'
// include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_LONG       } from '../modules/nf-core/minimap2/align/main'
// include { SAMTOOLS_INDEX                              } from '../modules/nf-core/samtools/index/main'
// include { SAMTOOLS_FAIDX                              } from '../modules/nf-core/samtools/faidx/main'
// include { BCFTOOLS_FILTER                             } from '../modules/nf-core/bcftools/filter/main'
// include { BCFTOOLS_CONSENSUS                          } from '../modules/nf-core/bcftools/consensus/main'
// include { SNPSITES                                    } from '../modules/nf-core/snpsites/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'

//
// NF-CORE SUBWORKFLOWS
//

// include { BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS } from '../subworkflows/nf-core/bam_variant_calling_sort_freebayes_bcftools/main'
// include { BAM_STATS_SAMTOOLS                          } from '../subworkflows/nf-core/bam_variant_calling_sort_freebayes_bcftools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BACTMAP {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    //
    // MODULE: Run FastQC
    //
    FASTQC(
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'bactmap_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
