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

include { GET_GENOME_SIZE                            } from '../modules/local/get_genome_size/main'
include { FASTQSCANPARSE as FASTQSCANPARSE_RAW       } from '../modules/local/fastq_scan_parse/main'
include { FASTQSCANPARSE as FASTQSCANPARSE_PROCESSED } from '../modules/local/fastq_scan_parse/main'
include { READ_STATS                                 } from '../modules/local/read_stats/main'
include { READSTATS_PARSE                            } from '../modules/local/read_stats_parse/main'
include { SEQTK_PARSE                                } from '../modules/local/seqtk_parse'
include { ALIGNPSEUDOGENOMES                         } from '../modules/local/alignpseudogenomes/main'

include { SHORTREAD_PREPROCESSING                    } from '../subworkflows/local/shortread_preprocessing/main'
include { LONGREAD_PREPROCESSING                     } from '../subworkflows/local/longread_preprocessing/main'
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
include { SAMTOOLS_FAIDX                         } from '../modules/nf-core/samtools/faidx/main'
include { GUNZIP                                 } from '../modules/nf-core/gunzip/main'
include { FASTQSCAN as FASTQSCAN_RAW             } from '../modules/nf-core/fastqscan/main'
include { FASTQC                                 } from '../modules/nf-core/fastqc/main'
include { FALCO                                  } from '../modules/nf-core/falco/main'
include { FASTQSCAN as FASTQSCAN_PROCESSED       } from '../modules/nf-core/fastqscan/main'
include { CAT_FASTQ as MERGE_RUNS                } from '../modules/nf-core/cat/fastq/main'
include { RASUSA                                 } from '../modules/nf-core/rasusa/main'
include { SNPSITES                               } from '../modules/nf-core/snpsites/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BACTMAP {

    take:
    samplesheet // channel: samplesheet read in from --input
    fasta       // channel: path(reference.fasta)
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    adapterlist = params.shortread_qc_adapterlist ? file(params.shortread_qc_adapterlist) : []
    custom_adapters = params.longread_qc_adapterlist ? file(params.longread_qc_adapterlist, checkIfExists: true) : []

    if ( params.shortread_qc_adapterlist ) {
        if ( params.shortread_qc_tool == 'adapterremoval' && !(adapterlist.extension == 'txt') ) error "[nf-core/bactmap] ERROR: AdapterRemoval2 adapter list requires a `.txt` format and extension. Check input: --shortread_qc_adapterlist ${params.shortread_qc_adapterlist}"
        if ( params.shortread_qc_tool == 'fastp' && !adapterlist.extension.matches(".*(fa|fasta|fna|fas)") ) error "[nf-core/bactmap] ERROR: fastp adapter list requires a `.fasta` format and extension (or fa, fas, fna). Check input: --shortread_qc_adapterlist ${params.shortread_qc_adapterlist}"
    }

    // Validate input files and create separate channels for FASTQ, FASTA, and Nanopore data
    ch_input = samplesheet
        .map { meta, fastq_1, fastq_2 ->

            // Define single_end based on the conditions
            if ( !fastq_1 ) {
                error("ERROR: Please check input samplesheet: entry `fastq_1` doesn't exist!")
            }
            meta.single_end = !fastq_2

            if (meta.instrument_platform == 'OXFORD_NANOPORE' && !meta.single_end) {
                error("Error: Please check input samplesheet: for Oxford Nanopore reads entry `fastq_2` should be empty!")
            }
            return [ meta, fastq_1, fastq_2 ]
        }
        .branch { meta, fastq_1, fastq_2 ->
            nanopore : meta.instrument_platform == 'OXFORD_NANOPORE'
                return [ meta + [type: "long"], [fastq_1]]
            fastq : meta.instrument_platform != 'OXFORD_NANOPORE'
                return [ meta + [ type: "short" ], fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ] ]
        }

    ch_input_for_fastqc = ch_input.nanopore.mix( ch_input.fastq )

    /*
        MODULE: Index reference file with Samtools faidx
    */

    ch_fasta = channel.of([])
    // Uncompress FASTA if needed
    if (fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP ([ [:], file(fasta, checkIfExists: true) ]).gunzip.map { tuple -> tuple[1] }
    } else {
        ch_fasta = channel.value(file(fasta, checkIfExists: true))
    }

    SAMTOOLS_FAIDX(ch_fasta.map { item -> [ [:], item, [] ] }, true )
    sizes = SAMTOOLS_FAIDX.out.sizes

    // Construct a reference fasta channel with empty meta for downstream subworkflows
    // Note: SAMTOOLS_FAIDX.out.fa only emits for .fa/.fasta extensions, so we use ch_fasta directly
    ch_fasta_meta = ch_fasta.map { item -> [ [:], item ] }

    /*
        MODULE: Get genome size
    */

    genome_size = GET_GENOME_SIZE(sizes).ch_genome_size

    /*
        Reference indexing
    */

    if (params.shortread_mapping_tool == 'bowtie2') {
        ch_index = BOWTIE2_BUILD ( ch_fasta.map { item -> [ [:], item ] } ).index
    } else {
        ch_index = BWAMEM2_INDEX ( ch_fasta.map { item -> [ [:], item ] }, 'bwamem2' ).index
    }

    /*
        MODULE: Run fastq-scan
    */
    FASTQSCAN_RAW (
        ch_input_for_fastqc
    )

    ch_fastqscanraw_fastqscanparse = FASTQSCAN_RAW.out.json
        .map { it[1] }
        .collect()

    ch_fastqscanraw_readstats = FASTQSCAN_RAW.out.json

    /*
        MODULE: Run fastqscanparse
    */
    FASTQSCANPARSE_RAW (
        ch_fastqscanraw_fastqscanparse
    )

    /*
        MODULE: Run FastQC
    */

    if ( !params.skip_preprocessing_qc ) {
        if ( params.preprocessing_qc_tool == 'falco' ) {
            FALCO ( ch_input_for_fastqc )
        } else {
            FASTQC ( ch_input_for_fastqc )
        }
    }

    /*
        SUBWORKFLOW: PERFORM PREPROCESSING
    */

    if (params.perform_shortread_qc) {
        SHORTREAD_PREPROCESSING(ch_input.fastq, adapterlist)
        ch_shortreads_preprocessed = SHORTREAD_PREPROCESSING.out.reads
    }
    else {
        ch_shortreads_preprocessed = ch_input.fastq
    }

    if ( params.perform_longread_qc ) {
        ch_longreads_preprocessed = LONGREAD_PREPROCESSING ( ch_input.nanopore, custom_adapters ).reads
                                        .map { it -> [ it[0], [it[1]] ] }
    } else {
        ch_longreads_preprocessed = ch_input.nanopore
    }

    ch_reads_for_fastqscan = ch_shortreads_preprocessed
            .mix( ch_longreads_preprocessed )

    /*
        MODULE: Run fastq-scan
    */
    FASTQSCAN_PROCESSED (
        ch_reads_for_fastqscan
    )

    ch_fastqscanprocessed_fastqscanparse = FASTQSCAN_PROCESSED.out.json
        .map { it[1] }
        .collect()

    /*
        MODULE: Run fastqscanparse
    */
    FASTQSCANPARSE_PROCESSED (
        ch_fastqscanprocessed_fastqscanparse
    )

    /*
        MODULE: Calculate read stats
    */
    ch_fastqscanraw_readstats                           // tuple val(meta), path(json)
        .join( FASTQSCAN_PROCESSED.out.json )           // tuple val(meta), path(json)
        .set { ch_readstats }                           // tuple val(meta), path(json), path(json)

    READ_STATS (
        ch_readstats
    )

    ch_readstats_readstatsparse = READ_STATS.out.csv
        .map { it[1] }
        .collect()

    /*
        MODULE: Summarise read stats outputs
    */
    READSTATS_PARSE (
        ch_readstats_readstatsparse
    )
    /*
        Run merging
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

    } else {
        ch_reads_runmerged = ch_shortreads_preprocessed
            .mix( ch_longreads_preprocessed )
    }

    /*
        MODULE: Perform subsampling
    */
    if ( params.perform_subsampling ) {
        ch_input_for_rasusa = ch_reads_runmerged
            .combine( genome_size )
        ch_reads_subsampled = RASUSA( ch_input_for_rasusa, params.subsampling_depth_cutoff ).reads
    } else {
        ch_reads_subsampled = ch_reads_runmerged
    }

    // Create separate channels for FASTQ and Nanopore data
    ch_mapping_input = ch_reads_subsampled
        .branch { meta, reads ->
                nanopore : meta.instrument_platform == 'OXFORD_NANOPORE'
                fastq : meta.instrument_platform != 'OXFORD_NANOPORE'
            }

    /*
        MODULE: Map short-reads
    */
    SHORTREAD_MAPPING (
        ch_mapping_input.fastq,
        ch_fasta_meta,
        ch_index,
        SAMTOOLS_FAIDX.out.fai
    )

    /*
        MODULE: Map long-reads
    */
    LONGREAD_MAPPING (
        ch_fasta_meta,
        SAMTOOLS_FAIDX.out.fai,
        ch_mapping_input.nanopore
    )

    /*
        MODULE: Summarise seqtk outputs
    */
    ch_seqtk_seqtkparse = SHORTREAD_MAPPING.out.seqtk_stats
        .mix( LONGREAD_MAPPING.out.seqtk_stats )

    SEQTK_PARSE (
        ch_seqtk_seqtkparse.map { tsv -> tsv[1] }.collect()
    )

    /*
        MODULE: Align pseudogenomes
    */
    ch_align_pseudogenomes = SHORTREAD_MAPPING.out.consensus
        .mix( LONGREAD_MAPPING.out.consensus )

    ALIGNPSEUDOGENOMES (
        ch_align_pseudogenomes.map { consensus -> consensus[1] }.collect(),
        ch_fasta_meta
    )

    ALIGNPSEUDOGENOMES.out.aligned_pseudogenomes
        .map { fasta_file ->
            def count = fasta_file.text.count('>')
            [count, fasta_file]
        }
        .branch { count, _fasta_file ->
            ALIGNMENT_NUM_PASS: count >= 4
            ALIGNMENT_NUM_FAIL: true
        }
        .set { aligned_pseudogenomes_branch }

    // Don't proceed further if too few genomes
    aligned_pseudogenomes_branch.ALIGNMENT_NUM_FAIL.view { count, _fasta ->
        "Insufficient (${count}) genomes after filtering to continue. Check results/pseudogenomes/low_quality_pseudogenomes.tsv for details"
    }

    aligned_pseudogenomes_branch.ALIGNMENT_NUM_PASS
        .map { _count, fasta_file -> fasta_file }
        .set { aligned_pseudogenomes }

    SNPSITES(
        aligned_pseudogenomes
    )

    /*
        Collate and save software versions
    */
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'bactmap_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    /*
        MODULE: MultiQC
    */
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def summary_params = paramsSummaryMap( parameters_schema: "nextflow_schema.json" )
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    if (!params.skip_preprocessing_qc) {
        if (params.preprocessing_qc_tool == 'falco') {
            // only mix in files actually used by MultiQC
            ch_multiqc_files = ch_multiqc_files.mix(
                FALCO.out.txt.map { _meta, reports -> reports }.flatten().filter { path -> path.name.endsWith('_data.txt') }.ifEmpty([])
            )
        }
        else {
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] }.ifEmpty([]))
        }
    }

    if (params.perform_shortread_qc) {
        ch_multiqc_files = ch_multiqc_files.mix( SHORTREAD_PREPROCESSING.out.mqc.collect{it[1]}.ifEmpty([]) )
    }

    if (params.perform_longread_qc) {
        ch_multiqc_files = ch_multiqc_files.mix( LONGREAD_PREPROCESSING.out.mqc.collect{it[1]}.ifEmpty([]) )
    }

    ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_MAPPING.out.mqc.collect{it[1]}.ifEmpty([]))

    ch_multiqc_files = ch_multiqc_files.mix(LONGREAD_MAPPING.out.mqc.collect{it[1]}.ifEmpty([]))

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'bactmap'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> report } // channel: /path/to/multiqc_report.html
    versions       = ch_versions                                        // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
