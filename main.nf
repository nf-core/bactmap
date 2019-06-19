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
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/bactmap --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --reference_sequence          The path to a fasta file the reference sequence to which the reads will be mapped
      --outdir                      Path to output dir
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.

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

/*
 * SET UP CONFIGURATION VARIABLES
 */


// Show help emssage
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
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}



// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
if(params.indir) summary['Reads'] = params.indir + "/" + params.fastq_pattern
if(params.accession_number_file) summary['File with accession numbers'] = params.accession_number_file 
summary['Reference sequence']    = params.reference
summary['Threshold for non GATC bases'] = non_GATC_bases_threshold
summary['Read coverage depth threshold'] = depth_cutoff
summary['VCF filtering conditions'] = filtering_conditions
summary['Remove recombination'] = params.remove_recombination
summary['Produce ML tree'] = params.tree
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(params.email) {
  summary['E-mail Address']  = params.email
  summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

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
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".csv") > 0) filename
        else null
    }

    output:
    file 'software_versions_mqc.yaml'
    file 'software_versions.csv'

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    bwa 2>&1 | grep Version > v_bwa.txt || true
    samtools 2>&1 | grep Version > v_samtools.txt || true
    bcftools 2>&1 | grep Version > v_bcftools.txt || true
    trimmomatic -version > v_trimmomatic.txt
    mash -version > v_mash.txt
    python3 -c "import pysam; print(pysam.__version__)" > v_pysam.txt
    seqtk 2>&1| grep Version > v_seqtk.txt || true
    iqtree -version  | grep version  > v_iqtree.txt
    snp-sites -V > v_snp-sites.txt
    gubbins 2>&1 | grep Version > v_gubbins.txt || true
    echo \$(fasttree 2>&1 | grep version) > v_fasttree.txt

    scrape_software_versions.py &> software_versions_mqc.yaml
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
 * STEP 0 - write out software versions
 */
get_software_versions()


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



/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/bactmap] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/bactmap] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.maxMultiqcEmailFileSize)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList){
                log.warn "[nf-core/bactmap] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/bactmap] Could not attach MultiQC report to summary email"
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/bactmap] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/bactmap] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCountFmt > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCountFmt} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCountFmt} ${c_reset}"
    }

    if(workflow.success){
        log.info "${c_purple}[nf-core/bactmap]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[nf-core/bactmap]${c_red} Pipeline completed with errors${c_reset}"
    }

}


def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/bactmap v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if(params.hostnames){
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if(hostname.contains(hname) && !workflow.profile.contains(prof)){
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
