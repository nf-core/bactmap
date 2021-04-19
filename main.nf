#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/bactmap
========================================================================================
 nf-core/bactmap Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/bactmap
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    // TODO nf-core: Update typical command used to run pipeline
    def command = "nextflow run nf-core/bactmap --input samplesheet.csv --reference ref.fasta -profile docker"
    log.info Schema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check that conda channels are set-up correctly
if (params.enable_conda) {
    Checks.check_conda_channels(log)
}

// Check AWS batch settings
Checks.aws_batch(workflow, params)

// Check the hostnames against configured profiles
Checks.hostname(workflow, params, log)

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
checkPathParamList = [ params.input, params.reference ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.reference) { ch_reference = file(params.reference) } else { exit 1, 'Reference fasta file not specified!' }

////////////////////////////////////////////////////
/* --       IMPORT MODULES / SUBWORKFLOWS      -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

// Local: Modules
include { GET_SOFTWARE_VERSIONS } from './modules/local/get_software_versions'     addParams( options: [publish_files : ['csv':'']] )
include { VCF2PSEUDOGENOME      } from './modules/local/vcf2pseudogenome'          addParams( options: modules['vcf2pseudogenome'])
include { ALIGNPSEUDOGENOMES } from './modules/local/alignpseudogenomes'           addParams( options: modules['alignpseudogenomes'])

// nf-core: Modules
include { GUBBINS } from './modules/nf-core/software/gubbins/main'                 addParams( options: modules['gubbins'])
include { SNPSITES } from './modules/nf-core/software/snpsites/main'               addParams( options: modules['snpsites'])

// Local: Sub-workflows
include { INPUT_CHECK       } from './modules/local/subworkflow/input_check'       addParams( options: [:] )
include { BAM_SORT_SAMTOOLS } from './modules/local/subworkflow/bam_sort_samtools' addParams( samtools_sort_options: modules['samtools_sort'], samtools_index_options : modules['samtools_index'], bam_stats_options: modules['bam_stats'])

def bcftools_filter_options = modules['bcftools_filter']
bcftools_filter_options.args = params.override_vcf_filter ? params.override_vcf_filter : bcftools_filter_options.args
println("BCF OPTIONS" + bcftools_filter_options)
include { VARIANTS_BCFTOOLS } from './modules/local/subworkflow/variants_bcftools' addParams( bcftools_mpileup_options: modules['bcftools_mpileup'], bcftools_filter_options: bcftools_filter_options)
include { SUB_SAMPLING } from './modules/local/subworkflow/sub_sampling'           addParams( mash_sketch_options: modules['mash_sketch'], rasusa_options: modules['rasusa'])
include { CREATE_PHYLOGENY } from './modules/local/subworkflow/create_phylogeny'   addParams( rapidnj_options: modules['rapidnj'], fasttree_options: modules['fasttree'], iqtree_options: modules['iqtree'], raxmlng_options: modules['raxmlng'])

include { find_genome_size } from './modules/local/functions.nf'

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////
def fastp_options   = modules['fastp']
if (fastp_options.adapter_fasta){
    fastp_options.args +=" --adapter_fasta adapter.fasta"
    ch_adapter_fasta = [ modules['fastp']['adapter_fasta'] ]
} else {
    ch_adapter_fasta = []
}
include { FASTP     } from './modules/nf-core/software/fastp/main'     addParams( options: fastp_options )
include { BWA_INDEX } from './modules/nf-core/software/bwa/index/main' addParams( options: modules['bwa_index'] )
include { BWA_MEM   } from './modules/nf-core/software/bwa/mem/main'   addParams( options: modules['bwa_mem'] )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report = []

workflow {

    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK ( 
        ch_input
    )

    /*
     * MODULE: Run bwa index
     */
    BWA_INDEX (
        ch_reference
    )
    /*
     * MODULE: Run fastp
     */
    if (params.trim){
        FASTP (
            INPUT_CHECK.out.sample_info,
            ch_adapter_fasta
        )
        ch_software_versions = ch_software_versions.mix(FASTP.out.version.first().ifEmpty(null))
        ch_reads = FASTP.out.reads
    } else {
        ch_reads = INPUT_CHECK.out.sample_info
    }

    if (params.depth_cutoff) {
        SUB_SAMPLING(ch_reads)
        ch_reads = SUB_SAMPLING.out.reads
    }
    

    BWA_MEM (
            ch_reads,
            BWA_INDEX.out.index
        )
        ch_software_versions = ch_software_versions.mix(BWA_MEM.out.version.first().ifEmpty(null))
    /*
     * SUBWORKFLOW: Sort bam files
     */
    BAM_SORT_SAMTOOLS ( 
        BWA_MEM.out.bam
    )
    ch_software_versions = ch_software_versions.mix(BAM_SORT_SAMTOOLS.out.samtools_version.first().ifEmpty(null))
    /*
     * SUBWORKFLOW: Call variants
     */
    VARIANTS_BCFTOOLS ( 
        BAM_SORT_SAMTOOLS.out.bam,
        ch_reference
    )
    ch_software_versions = ch_software_versions.mix(VARIANTS_BCFTOOLS.out.bcftools_version.first().ifEmpty(null))

    /*
     * MODULE: Make pseudogenome from VCF
     */
    VCF2PSEUDOGENOME (
        VARIANTS_BCFTOOLS.out.filtered_vcf,
        ch_reference
    )
    ch_software_versions = ch_software_versions.mix(VCF2PSEUDOGENOME.out.version.first().ifEmpty(null))

    /*
     * MODULE: make pseudogenome alignment
     */
    ALIGNPSEUDOGENOMES (
        VCF2PSEUDOGENOME.out.pseudogenome.map { pseudogenome -> pseudogenome[1] }.collect(),
        ch_reference
    )
    ch_software_versions = ch_software_versions.mix(ALIGNPSEUDOGENOMES.out.version.ifEmpty(null))
    ALIGNPSEUDOGENOMES.out.aligned_pseudogenomes
        .branch { 
            aligned_pseudogenomes ->
            ALIGNMENT_NUM_PASS: aligned_pseudogenomes[0].toInteger() >= 4
            ALIGNMENT_NUM_FAIL: aligned_pseudogenomes[0].toInteger() < 4
        }
        .set { aligned_pseudogenomes_branch }
    
    // Don't proceeed further if two few genonmes
    aligned_pseudogenomes_branch.ALIGNMENT_NUM_FAIL.view { "Insufficient (${it[0]}) genomes after filtering to continue. Check results/pseudogenomes/low_quality_pseudogenomes.tsv for details"}

    aligned_pseudogenomes_branch.ALIGNMENT_NUM_PASS
        .map{ it[1] }
        .set { aligned_pseudogenomes }

    SNPSITES(
        aligned_pseudogenomes
    )
    ch_software_versions = ch_software_versions.mix(SNPSITES.out.version.ifEmpty(null))
    
    /*
     * MODULE: remove recombination
     */
    if (params.remove_recombination){
        GUBBINS (
            aligned_pseudogenomes
        )
        ch_software_versions = ch_software_versions.mix(GUBBINS.out.version.ifEmpty(null))
        /*
        * SUBWORKFLOW: Create phylogenies
        */
        CREATE_PHYLOGENY (
            GUBBINS.out.fasta,
            SNPSITES.out.constant_sites_string
        )
    } else {
        /*
        * SUBWORKFLOW: Create phylogenies
        */
        CREATE_PHYLOGENY (
            SNPSITES.out.fasta,
            SNPSITES.out.constant_sites_string
        )
    }
    ch_software_versions = ch_software_versions.mix(CREATE_PHYLOGENY.out.rapidnj_version.ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(CREATE_PHYLOGENY.out.fasttree_version.ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(CREATE_PHYLOGENY.out.iqtree_version.ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(CREATE_PHYLOGENY.out.raxmlng_version.ifEmpty(null))

    /*
     * MODULE: Pipeline reporting
     */
    GET_SOFTWARE_VERSIONS ( 
        ch_software_versions.map { it }.collect()
    )

}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
