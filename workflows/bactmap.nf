/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_bactmap_pipeline'

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.reference, params.multiqc_config,
                           params.shortread_qc_adapterlist, params.multiqc_logo, 
                           params.multiqc_methods_description ]
                            
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if ( params.input ) {
    ch_input = file(params.input, checkIfExists: true)
} else {
    error("Input samplesheet not specified")
}

if (params.reference) { ch_reference =  Channel.fromPath(params.reference) } else { exit 1, 'Reference sequence FASTA not specified!' }

// Modify reference channel to include meta data
ch_reference_meta = ch_reference.map{ it -> [[id:it[0].baseName], it] }.collect()

// Get genome size from reference
//TO DO!!
genome_size = 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { FASTQSCANPARSE as FASTQSCANPARSE_TRIM      } from '../modules/local/fastqscanparse/main'
include { FASTQSCANPARSE as FASTQSCANPARSE_SUBSAMPLE } from '../modules/local/fastqscanparse/main'
include { SEQTK_COMP                                 } from '../modules/local/seqtk_comp/main'
include { SEQTK_PARSE                                } from '../modules/local/seqtk_parse/main'
include { ALIGNPSEUDOGENOMES                         } from '../modules/local/alignpseudogenomes/main'

include { SHORTREAD_PREPROCESSING                    } from '../subworkflows/local/shortread_preprocessing'
include { LONGREAD_PREPROCESSING                     } from '../subworkflows/local/longread_preprocessing'
include { SHORTREAD_MAPPING                          } from '../subworkflows/local/shortread_mapping/main'
include { LONGREAD_MAPPING                           } from '../subworkflows/local/longread_mapping/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// NF-CORE MODULES/PLUGINS
//
include { BOWTIE2_BUILD                          } from '../modules/nf-core/bowtie2/build/main'
include { BWAMEM2_INDEX                          } from '../modules/nf-core/bwamem2/index/main'
include { FASTQC                                 } from '../modules/nf-core/fastqc/main'
include { FALCO                                  } from '../modules/nf-core/falco/main'
include { CAT_FASTQ as MERGE_RUNS                } from '../modules/nf-core/cat/fastq/main'
include { FASTQSCAN as FASTQSCAN_TRIM            } from '../modules/nf-core/modules/fastqscan/main'
include { RASUSA                                 } from '../modules/nf-core/rasusa/main'
include { FASTQSCAN as FASTQSCAN_SUBSAMPLE       } from '../modules/nf-core/modules/fastqscan/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_SHORT } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_LONG  } from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_INDEX                         } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX                         } from '../modules/nf-core/samtools/faidx/main'
include { BCFTOOLS_FILTER                        } from '../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_CONSENSUS                     } from '../modules/nf-core/bcftools/consensus/main'
include { SNPSITES                               } from '../modules/nf-core/snpsites/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'

//
// NF-CORE SUBWORKFLOWS
//

include { BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS } from '../subworkflows/nf-core/bam_variant_calling_sort_freebayes_bcftools/main'
include { BAM_STATS_SAMTOOLS } from '../subworkflows/nf-core/bam_variant_calling_sort_freebayes_bcftools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BACTMAP {

    take:
    samplesheet  // channel: samplesheet read in from --input
    ch_reference // channel: path(reference.fasta)
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    // Validate input files and create separate channels for FASTQ, FASTA, and Nanopore data
    ch_input = samplesheet
        .map { meta, run_accession, instrument_platform, fastq_1, fastq_2, fasta ->
            meta.run_accession = run_accession
            meta.instrument_platform = instrument_platform

            // Define single_end based on the conditions
     if ( !fastq_1 ) {
            error("ERROR: Please check input samplesheet: entry `fastq_1` doesn't exist!")
     meta.single_end = !fastq_2
     if (meta.single_end && meta.instrument_platform == 'OXFORD_NANOPORE') {
          error("Error: Please check input samplesheet: for Oxford Nanopore reads entry `fastq_2` should be empty!")
     
        }
        .branch { meta, fastq_1, fastq_2 ->
           nanopore : meta.instrument_platform == 'OXFORD_NANOPORE'
                return [ meta + [type: "long"], [fastq_1, fastq_2]
           fastq : meta.instrument_platform != 'OXFORD_NANOPORE'
               return [ meta + [ type: "short" ], fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ] ]
           }
        }
    ch_input_for_fastqc = ch_input.fastq.mix( ch_input.fastq )
    
    /*
        Reference indexing
    */
    if (params.shortread_mapping_tool == 'bowtie2') {
        ch_index = BOWTIE_BUILD ( ch_reference ).index
        ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions.first())
    } else {
        ch_index = BWAMEM2_INDEX ( ch_reference ).index
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())
    }
    
    /*
        MODULE: Run FastQC
    */

    if ( !params.skip_preprocessing_qc ) {
        if ( params.preprocessing_qc_tool == 'falco' ) {
            FALCO ( ch_input_for_fastqc )
            ch_versions = ch_versions.mix(FALCO.out.versions.first())
        } else {
            FASTQC ( ch_input_for_fastqc )
            ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        }
    }

    /*
        SUBWORKFLOW: PERFORM PREPROCESSING
    */

    if ( params.perform_shortread_qc ) {
        ch_shortreads_preprocessed = SHORTREAD_PREPROCESSING ( ch_input.fastq, adapterlist ).reads
        ch_versions = ch_versions.mix( SHORTREAD_PREPROCESSING.out.versions )
    } else {
        ch_shortreads_preprocessed = ch_input.fastq
    }

    if ( params.perform_longread_qc ) {
        ch_longreads_preprocessed = LONGREAD_PREPROCESSING ( ch_input.nanopore ).reads
                                        .map { it -> [ it[0], [it[1]] ] }
        ch_versions = ch_versions.mix( LONGREAD_PREPROCESSING.out.versions )
    } else {
        ch_longreads_preprocessed = ch_input.nanopore
    }

    /*
        Run indexing
    */
    if ( params.perform_runmerging ) {

        ch_reads_for_cat_branch = ch_shortreads_preprocessed
            .mix( ch_longreads_preprocessed )
            .map {
                meta, reads ->
                    def meta_new = meta - meta.subMap('run_accession')
                    [ meta_new, reads ]
            }
            .groupTuple()
            .map {
                meta, reads ->
                    [ meta, reads.flatten() ]
            }
            .branch {
                meta, reads ->
                // we can't concatenate files if there is not a second run, we branch
                // here to separate them out, and mix back in after for efficiency
                cat: ( meta.single_end && reads.size() > 1 ) || ( !meta.single_end && reads.size() > 2 )
                skip: true
            }

        ch_reads_runmerged = MERGE_RUNS ( ch_reads_for_cat_branch.cat ).reads
            .mix( ch_reads_for_cat_branch.skip )
            .map {
                meta, reads ->
                [ meta, [ reads ].flatten() ]
            }
            .mix( ch_input.fasta_short, ch_input.fasta_long)

        ch_versions = ch_versions.mix(MERGE_RUNS.out.versions)

    } else {
        ch_reads_runmerged = ch_shortreads_preprocessed
            .mix( ch_longreads_preprocessed, ch_input.fasta_short, ch_input.fasta_long )
    }
    
    /*
        MODULE: Run fastq-scan
    */
    FASTQSCAN_TRIM (
        ch_reads_runmerged
    )
    ch_fastqscantrim_fastqscanparse = FASTQSCAN_RAW.out.json
    ch_fastqscantrim_readstats      = FASTQSCAN_RAW.out.json
    ch_versions                     = ch_versions.mix(FASTQSCAN_TRIM.out.versions.first())
    
    /*
        MODULE: Run fastqscanparse
    */
    FASTQSCANPARSE_TRIM (
        ch_fastqscantrim_fastqscanparse.collect{it[1]}.ifEmpty([])
    )
    ch_versions = ch_versions.mix(FASTQSCANPARSE_TRIM.out.versions.first())
    
    /*
        MODULE: Perform subsampling
    */
     
    if ( params.perform_subsampling ) {
        ch_reads_subsampled = RASUSA( ch_reads_runmerged, genome_size, subsampling_depth_cutoff ).reads
        ch_versions = ch_versions.mix( RASUSA.out.versions )
    } else {
        ch_reads_subsampled = ch_reads_runmerged
    }
    
    /*
        MODULE: Run fastq-scan
    */
    FASTQSCAN_SUBSAMPLE (
        ch_reads_subsampled
    )
    ch_fastqscansubsample_fastqscanparse = FASTQSCAN_SUBSAMPLE.out.json
    ch_fastqscansubsample_readstats      = FASTQSCAN_SUBSAMPLE.out.json
    ch_versions                          = ch_versions.mix(FASTQSCAN_SUBSAMPLE.out.versions.first())
    
    /*
        MODULE: Run fastqscanparse
    */
    FASTQSCANPARSE_SUBSAMPLE (
        ch_fastqscansubsample_fastqscanparse.collect{it[1]}.ifEmpty([])
    )
    ch_versions = ch_versions.mix(FASTQSCANPARSE_SUBSAMPLE.out.versions.first())
    
    /*
        MODULE: Map reads
    */


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'bactmap_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
