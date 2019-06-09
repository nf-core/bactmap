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
// DSL 2
nextflow.preview.dsl=2

def helpMessage() {
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'

     nf-core/bactmap v${workflow.manifest.version}
    =======================================================

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/bactmap --reads '*_R{1,2}.fastq.gz' -profile standard,docker

    Mandatory arguments:
      --reference_sequence  The path to a fasta file the reference sequence to which the reads will be mapped
      --outdir          Path to output dir
      -profile              Configuration profile to use. Can use multiple (comma separated)
                            Available: standard, docker, test

    Alternate mandatory arguments:
    Either indir and fastq_pattern must be specified if using local short reads or accession_number_file if
      specifying samples in the SRA for which the fastqs will be downloaded

      --indir                 Path to input dir. This must be used in conjunction with fastq_pattern
      --fastq_pattern             The regular expression that will match fastq files wrapped in quotes e.g '*{R,_}{1,2}.fastq.gz'
      --accession_number_file     Path to a text file containing a list of accession numbers (1 per line)
      
    Optional arguments:
      --adapter_sequences         The path to a fasta file containing adapter sequences to trim from reads. If not specified
                                  those in the workflow basedir will be used
      --depth_cutoff              The estimated depth to downsample each sample to. If not specified no downsampling will occur
      --filtering_conditions      A string containing the expression used to filter the vcf files
                                  If this is not specified the default is :
                                  %QUAL<25 || FORMAT/DP<10 || MAX(FORMAT/ADF)<5 || MAX(FORMAT/ADR)<5 || MAX(FORMAT/AD)/SUM(FORMAT/DP)<0.9 || MQ<30 || MQ0F>0.1
      --tree                      Whether to create a maximum likelihood tree
      --remove_recombination      Whether to remove recombination from the combined alignment using gubbins before
                                  producing the ML tree
      --non_GATC_bases_threshold  The maximum fraction of non-GATC bases in a pseudogenome to accept. Those not passing will be
                                  written into a file $outdir/pseudogenomes/low_quality_pseudogenomes.tsv. Default is 0.5

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

// Show help emssage


/*
 * SET UP VARIABLES
 */
params.help = false
params.indir = false
params.accession_number_file = false
params.fastq_pattern = false
params.adapter_file = false
params.outdir = false
params.reference = false
params.depth_cutoff = false
params.tree = false
params.remove_recombination = false
params.filtering_conditions = false
params.non_GATC_bases_threshold = false

if (params.help){
    helpMessage()
    exit 0
}
include 'lib/helper_functions'

check_optional_parameters(params, ['indir', 'accession_number_file', 'read_paths'])

//  check a pattern has been specified
if (params.indir){
  fastq_pattern = check_mandatory_parameter(params, 'fastq_pattern')
}

// set up path to adapter sequences 
if ( params.adapter_file ) {
  adapter_file = params.adapter_file
} else {
  adapter_file = "${workflow.projectDir}/adapters.fas"
}

// set up output directory
outdir = check_mandatory_parameter(params, 'outdir') - ~/\/$/

// set up input from reference sequnce
reference_sequence = file(check_mandatory_parameter(params, 'reference'))

// assign depth cutoff
if ( params.depth_cutoff ) {
  depth_cutoff = params.depth_cutoff
} else {
  depth_cutoff = false
}

// assign filtering parameters
if ( params.filtering_conditions ) {
  filtering_conditions = params.filtering_conditions
} else {
  filtering_conditions = '%QUAL<25 || FORMAT/DP<10 || MAX(FORMAT/ADF)<5 || MAX(FORMAT/ADR)<5 || MAX(FORMAT/AD)/SUM(FORMAT/DP)<0.9 || MQ<30 || MQ0F>0.1'
}

// assign non GATC bases threshold
if ( params.non_GATC_bases_threshold ) {
  non_GATC_bases_threshold = params.non_GATC_bases_threshold
} else {
  non_GATC_bases_threshold = 0.5
}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
  // Check workDir/outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}


// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/bactmap v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/bactmap'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
if(params.indir) summary['Reads'] = params.indir + "/" + params.fastq_pattern
if(params.accession_number_file) summary['File with accession numbers'] = params.accession_number_file 
summary['Reference sequence']    = params.reference_sequence
summary['Threshold for non GATC bases'] = non_GATC_bases_threshold
summary['Read coverage depth threshold'] = depth_cutoff
summary['VCF filtering conditions'] = filtering_conditions
summary['Remove recombination'] = params.remove_recombination
summary['Produce ML tree'] = params.tree
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-bactmap-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/bactmap Workflow Summary'
    section_href: 'https://github.com/nf-core/bactmap'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
INCLUDE processes from module file
*/
include 'lib/bactmap_processes' params(
    outdir: outdir,
    depth_cutoff: depth_cutoff,
    filtering_conditions: filtering_conditions,
    non_GATC_bases_threshold: non_GATC_bases_threshold,
    remove_recombination: params.remove_recombination,
    cpus_for_tree_building: params.cpus_for_tree_building,
    memory_for_tree_building: params.memory_for_tree_building
  )

/*
 * STEP 1 - index reference sequence for subsequent mapping
 */

prepared_reference_files = prepare_reference(reference_sequence)

/*
 * STEP 2 - setup fastq reads from either local diretory or by fetching from the ENA
 */

if (params.accession_number_file){
  // Fetch samples from ENA
  accession_number_file = params.accession_number_file - ~/\/$/
  accession_numbers = Channel
      .fromPath(accession_number_file)
      .splitText()
      .map{ x -> x.trim()}
  
  raw_fastqs = fetch_from_ena(accession_numbers)

} else if (params.read_paths) {
  raw_fastqs = Channel.from( params.read_paths )
    .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
} else if (params.indir) {
  indir = params.indir - ~/\/$/
  fastqs = indir + '/' + fastq_pattern
  raw_fastqs = Channel
    .fromFilePairs( fastqs )
    .ifEmpty { error "Cannot find any reads matching: ${fastqs}" }
}

/*
 * STEP 3 - Assess read length prior to trimming. Will set  MIN LEN for trimmomatic 1/3 of this value
 */

min_read_length = determine_min_read_length(raw_fastqs)
min_read_length_and_raw_fastqs = min_read_length.join(raw_fastqs)

/*
 * STEP 4 - Trim reads
 */

trimmed_fastqs = trim_reads(min_read_length_and_raw_fastqs, adapter_file)

/*
 * STEP 5 - Estimate genome size based on kmers with mash
 */

mash_output = genome_size_estimation(trimmed_fastqs)
// channel to output genome size from mash output
genome_size_estimation_for_downsampling =  mash_output.map { sample_id, file -> find_genome_size(sample_id, file.text) } 

// /*
//  * STEP 6 - Estimate total number of reads
// */

seqtk_fqchk_output = count_number_of_reads(trimmed_fastqs)
read_counts = seqtk_fqchk_output.map { sample_id, file -> find_total_number_of_reads(sample_id, file.text) }

fastqs_and_genome_size_and_read_count = trimmed_fastqs.join(genome_size_estimation_for_downsampling).join(read_counts).map{ tuple -> [tuple[0], tuple[1], tuple[2], tuple[3], tuple[4]]}

/*
===========================================
  run snp pipeline to create filtered vcfs
===========================================
*/

/*
 * STEP 7 - Map reads with bwa
*/

sorted_bam_files = map_reads(reference_sequence, prepared_reference_files, fastqs_and_genome_size_and_read_count)

// /*
//  * STEP 8 - Call variants
// */

bcf_files = call_variants(reference_sequence, sorted_bam_files)

/*
 * STEP 9 - Filter variants
*/
filtered_bcf_files = filter_variants(bcf_files)


// /*
//  * STEP 10 - Create pseudogenome from filtered bcf file
// */
pseudogenomes = create_pseudogenome(reference_sequence, filtered_bcf_files)

// /*
//  * STEP 11 - Combine all sample pseudogenomes into an fasta alignment
// */
combined_pseudogenomes = pseudogenomes.collect { it }

(aligned_pseudogenome, low_quality_pseudogenomes) = create_pseudogenome_alignment(reference_sequence, combined_pseudogenomes)

/*
 * STEP 12 - Remove non-informative positions from alignment
*/
variant_only_aligned_pseudogenome = create_variant_only_alignment(aligned_pseudogenome)
/*
 * STEP 13 - Build a ML tree using IQTree
*/

if (params.tree) {
  build_tree(variant_only_aligned_pseudogenome)
}

workflow.onComplete {
  complete_message(params, workflow, version)
}

workflow.onError {
  error_message(workflow)
}
