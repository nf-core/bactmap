/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowBactmap.initialise(params, log)

checkPathParamList = [ params.input, params.reference, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.reference) { ch_reference = file(params.reference) } else { exit 1, 'Reference fasta file not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

// Local: Modules
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['csv':'']] )
include { VCF2PSEUDOGENOME      } from '../modules/local/vcf2pseudogenome'      addParams( options: modules['vcf2pseudogenome'])
include { ALIGNPSEUDOGENOMES    } from '../modules/local/alignpseudogenomes'    addParams( options: modules['alignpseudogenomes'])
include { MULTIQC               } from '../modules/local/multiqc'               addParams( options: multiqc_options )

// Local: Sub-workflows
include { INPUT_CHECK       } from '../subworkflows/input_check'       addParams( options: [:] )
include { BAM_SORT_SAMTOOLS } from '../subworkflows/bam_sort_samtools' addParams( samtools_sort_options: modules['samtools_sort'], samtools_index_options : modules['samtools_index'], bam_stats_options: modules['bam_stats'])

include { VARIANTS_BCFTOOLS } from '../subworkflows/variants_bcftools' addParams( bcftools_mpileup_options: modules['bcftools_mpileup'], bcftools_filter_options: modules['bcftools_filter'])
include { SUB_SAMPLING      } from '../subworkflows/sub_sampling'      addParams( mash_sketch_options: modules['mash_sketch'], rasusa_options: modules['rasusa'])

include { CREATE_PHYLOGENY  } from '../subworkflows/create_phylogeny'   addParams( rapidnj_options: modules['rapidnj'], fasttree_options: modules['fasttree'], iqtree_options: modules['iqtree'], raxmlng_options: modules['raxmlng'] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def fastp_options = modules['fastp']
if (params.trim && params.adapter_file){
    fastp_options.args +=" --adapter_fasta adapter.fasta"
    ch_adapter_file = [ params.adapter_file ]
} else {
    ch_adapter_file = []
}
include { FASTP     } from '../modules/nf-core/software/fastp/main'     addParams( options: fastp_options )
include { BWA_INDEX } from '../modules/nf-core/software/bwa/index/main' addParams( options: modules['bwa_index'] )
include { BWA_MEM   } from '../modules/nf-core/software/bwa/mem/main'   addParams( options: modules['bwa_mem'] )
include { GUBBINS   } from '../modules/nf-core/software/gubbins/main'   addParams( options: modules['gubbins'] )
include { SNPSITES  } from '../modules/nf-core/software/snpsites/main'  addParams( options: modules['snpsites'] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow BACTMAP {
    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )

    //
    // MODULE: Run bwa index
    //
    BWA_INDEX (
        ch_reference
    )

    //
    // MODULE: Run fastp
    //
    if (params.trim && params.adapter_file){
        FASTP (
            INPUT_CHECK.out.sample_info,
            ch_adapter_file
        )
        ch_software_versions = ch_software_versions.mix(FASTP.out.version.first().ifEmpty(null))
        ch_reads = FASTP.out.reads
    } else {
        ch_reads = INPUT_CHECK.out.sample_info
    }

    if (!params.subsampling_off) {
        SUB_SAMPLING(
            ch_reads
        )
        ch_reads = SUB_SAMPLING.out.reads
    }

    //
    // MODULE: Map reads
    //
    BWA_MEM (
        ch_reads,
        BWA_INDEX.out.index
    )
    ch_software_versions = ch_software_versions.mix(BWA_MEM.out.version.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Sort bam files
    //
    BAM_SORT_SAMTOOLS (
        BWA_MEM.out.bam
    )
    ch_software_versions = ch_software_versions.mix(BAM_SORT_SAMTOOLS.out.samtools_version.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Call variants
    //
    VARIANTS_BCFTOOLS (
        BAM_SORT_SAMTOOLS.out.bam,
        ch_reference
    )
    ch_software_versions = ch_software_versions.mix(VARIANTS_BCFTOOLS.out.bcftools_version.first().ifEmpty(null))

    //
    // MODULE: Make pseudogenome from VCF
    //
    VCF2PSEUDOGENOME (
        VARIANTS_BCFTOOLS.out.filtered_vcf,
        ch_reference
    )

    //
    // MODULE: make pseudogenome alignment
    //
    ALIGNPSEUDOGENOMES (
        VCF2PSEUDOGENOME.out.pseudogenome.map { pseudogenome -> pseudogenome[1] }.collect(),
        ch_reference
    )
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

    //
    // MODULE: remove recombination
    //
    if (params.remove_recombination){
        GUBBINS (
            aligned_pseudogenomes
        )
        ch_software_versions = ch_software_versions.mix(GUBBINS.out.version.ifEmpty(null))

        //
        // SUBWORKFLOW: Create phylogenies
        //
        CREATE_PHYLOGENY (
            GUBBINS.out.fasta,
            SNPSITES.out.constant_sites_string
        )
    } else {
        //
        // SUBWORKFLOW: Create phylogenies
        //
        CREATE_PHYLOGENY (
            SNPSITES.out.fasta,
            SNPSITES.out.constant_sites_string
        )
    }
    ch_software_versions = ch_software_versions.mix(CREATE_PHYLOGENY.out.rapidnj_version.ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(CREATE_PHYLOGENY.out.fasttree_version.ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(CREATE_PHYLOGENY.out.iqtree_version.ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(CREATE_PHYLOGENY.out.raxmlng_version.ifEmpty(null))

    //
    // MODULE: Pipeline reporting
    //
    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowBactmap.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            GET_SOFTWARE_VERSIONS.out.yaml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            FASTP.out.json.collect{it[1]}.ifEmpty([]),
            BAM_SORT_SAMTOOLS.out.stats.collect{it[1]}.ifEmpty([]),
            VARIANTS_BCFTOOLS.out.stats.collect{it[1]}.ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
