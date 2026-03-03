//
// Subworkflow with functionality specific to the nf-core/bactmap pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/bactmap ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/taxprofiler/blob/master/CITATIONS.md
"""
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    channel
        .fromList(samplesheetToList(params.input, "assets/schema_input.json"))
        .set { ch_samplesheet }


    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def text_seq_qc = [
        "Sequencing quality control with",
        params.preprocessing_qc_tool == "falco" ? "Falco (de Sena Brandine and Smith 2021)." : "FastQC (Andrews 2010)."
    ].join(' ').trim()

    def text_shortread_qc = [
        "Short read preprocessing was performed with:",
        params.shortread_qc_tool == "adapterremoval" ? "AdapterRemoval (Schubert et al. 2016)." : "",
        params.shortread_qc_tool == "fastp" ? "fastp (Chen et al. 2018)." : "",
    ].join(' ').trim()

    def text_longread_qc = [
        "Long read preprocessing was performed with:",
        params.longread_adapterremoval_tool == "porechop_abi" ? "Porechop_ABI (Bonenfant et al. 2023)," : "",
        params.longread_adapterremoval_tool == "porechop" ? "Porechop (Wick et al. 2017)," : "",
        params.longread_filter_tool == "filtlong" ? "Filtlong (Wick 2021)." : "",
        params.longread_filter_tool == "nanoq" ? "Nanoq (Steinig and Coin 2022)." : "",
    ].join(' ').trim()

    def text_subsampling = [
        "Read subsampling was done with Rasusa (Hall et al. 2019)."
    ].join(' ').trim()

    def shortread_mapping = [
        "Short read mapping and variant calling was performed with",
        params.shortread_mapping_tool == "bowtie2" ? "Bowtie2 (Langmead and Salzberg 2012)," : "BWA-MEM2 (Li 2021),",
        "FreeBayes (Garrison and Marth 2012)",
        "and BCFtools (Li 2011).",
    ].join(' ').trim()

    def longread_mapping = [
        "Long read mapping and variant calling was performed with",
        "Minimap2 (Li 2018),",
        "Clair3 (Zhang et al. 2020)",
        "and BCFtools (Li 2011).",
    ].join(' ').trim()

    def consensus = [
        "Consensus FASTA sequences were created with",
        "BCFtools (Li 2011)",
        "and bedtools (Quinlan and Hall 2010).",
    ].join(' ').trim()

    def snps = [
        "Variant and constant sites were extracted from the alignment using SNP-sites (Page et al. 2016).",
    ].join(' ').trim()

    def statistics = [
        "Read, mapping and variant statistics were generated with",
        "fastq-scan (Petit 2022),",
        "Samtools (Li et al. 2009),",
        "BCFtools (Li 2011)",
        "and seqtk (Li 2012).",
    ].join(' ').trim()

    def citation_text = [
        "Tools used in the workflow included:",
        text_seq_qc,
        params.perform_shortread_qc ? text_shortread_qc : "",
        params.perform_longread_qc ? text_longread_qc : "",
        params.perform_subsampling ? text_subsampling : "",
        shortread_mapping,
        longread_mapping,
        consensus,
        snps,
        statistics,
        "Pipeline results statistics were summarised with MultiQC (Ewels et al. 2016)."
        ].join(' ').trim().replaceAll("[,|.] +\\.", ".")

    return citation_text
}

def toolBibliographyText() {
    def text_seq_qc = [
        params.preprocessing_qc_tool == "falco"  ? "<li>de Sena Brandine, G., & Smith, A. D. (2021). Falco: high-speed FastQC emulation for quality control of sequencing data. F1000Research, 8(1874), 1874.  <a href=\"https://doi.org/10.12688/f1000research.21142.2\">10.12688/f1000research.21142.2</li>" : "",
        params.preprocessing_qc_tool == "fastqc" ? "<li>Andrews S. (2010) FastQC: A Quality Control Tool for High Throughput Sequence Data, URL: <a href=\"https://www.bioinformatics.babraham.ac.uk/projects/fastqc/\">https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</a></li>" : "",
    ].join(' ').trim()

    def text_shortread_qc = [
        params.shortread_qc_tool == "adapterremoval" ? "<li>Schubert, M., Lindgreen, S., & Orlando, L. (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. BMC Research Notes, 9, 88. <a href=\"https://doi.org/10.1186/s13104-016-1900-2\">10.1186/s13104-016-1900-2</a></li>" : "",
    ].join(' ').trim()

    def text_longread_qc = [
        params.longread_adapterremoval_tool == "porechop_abi" ? "<li>Bonenfant, Q., Noé, L., & Touzet, H. (2023). Porechop_ABI: discovering unknown adapters in Oxford Nanopore Technology sequencing reads for downstream trimming. Bioinformatics Advances, 3(1):vbac085. <a href=\"https://10.1093/bioadv/vbac085\">10.1093/bioadv/vbac085</a></li>" : "",
        params.longread_adapterremoval_tool == "porechop" ? "<li>Wick, R. R., Judd, L. M., Gorrie, C. L., & Holt, K. E. (2017). Completing bacterial genome assemblies with multiplex MinION sequencing. Microbial Genomics, 3(10), e000132. <a href=\"https://doi.org/10.1099/mgen.0.000132\">10.1099/mgen.0.000132</a></li>" : "",
        params.longread_filter_tool == "filtlong" ? "<li>Wick R. (2021) Filtlong, URL:  <a href=\"https://github.com/rrwick/Filtlong\">https://github.com/rrwick/Filtlong</a></li>" : "",
        params.longread_filter_tool == "nanoq" ? "<li>Steinig, E., & Coin, L. (2022). Nanoq: ultra-fast quality control for nanopore reads. Journal of Open Source Software, 7(69). <a href=\"https://doi.org/10.21105/joss.02991\">10.21105/joss.02991</a></li>" : ""
    ].join(' ').trim()

    def text_subsampling = [
        "<li>Hall, M. B. (2019). Rasusa: Randomly subsample sequencing reads to a specified coverage. <a href=\"https:/doi.org/10.5281/zenodo.3731394\">10.5281/zenodo.3731394</a></li>"
    ].join(' ').trim()

    def text_shortreadmapping = [
        params.shortread_mapping_tool == "bowtie2" ? "<li>Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357–359. <a href=\"https://doi.org/10.1038/nmeth.1923\">10.1038/nmeth.1923</a></li>" : "<li>Vasimuddin, M., Misra, S., Li H., & Aluru S. (2019). Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE International Parallel and Distributed Processing Symposium (IPDPS). 2019, pp. 314-324. <a href=\"https://doi.org/10.1109/IPDPS.2019.00041\">10.1109/IPDPS.2019.00041</a></li>",
        "<li>Garrison, E., & Marth, G. T. (2012). Haplotype-based variant detection from short-read sequencing. arXiv, 1207.3907 <a href=\"https://doi.org/10.48550/arXiv.1207.3907\">10.48550/arXiv.1207.3907</a></li>",
        "<li>Li, H. (2011). A statistical framework for SNP calling, mutation discovery, and population genetical parameter estimation from sequencing data. Bioinformatics, 27(21), 2987–2993. <a href=\"https://doi.org/10.1093/bioinformatics/btr509\">10.1093/bioinformatics/btr509</a></li>",
    ].join(' ').trim()

    def text_longreadmapping = [
        "<li>Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics , 34(18), 3094–3100. <a href=\"https://doi.org/10.1093/bioinformatics/bty191\">10.1093/bioinformatics/bty191</a></li>",
        "<li>Zheng Z., Li S., Su J., Leung A. W., Lam T. W., & Luo R. (2022). Clair3: fast and accurate long-read variant calling with deep neural networks. Nature Computer Science, 2(12):797-803. <a href=\"https://doi.org/10.1038/s43588-022-00387-x\">10.1038/s43588-022-00387-x</a></li>",
    ].join(' ').trim()

    def text_consensus = [
        "<li>Quinlan A. R., & Hall I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6):841-2. <a href=\"https:/doi.org/10.1093/bioinformatics/btq033\">10.1093/bioinformatics/btq033</a></li>"
    ].join(' ').trim()

    def text_snps = [
        "<li>Page A. J., Taylor B., Delaney A. J., Soares J., Seemann T., Keane J. A., Harris S. R. (2016). SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments. Microbial Genomics, 29;2(4):e000056. <a href=\"https:/doi.org/10.1099/mgen.0.000056\">10.1099/mgen.0.000056</a></li>"
    ].join(' ').trim()

    def text_statistics = [
        "<li>Petit, R. (2022). fastq-scan: A tools for reading a FASTQ from STDIN and outputting summary statistics. URL: <a href=\"https://github.com/rpetit3/fastq-scan\">https://github.com/rpetit3/fastq-scan</a></li>",
        "<li>Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., ... & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. <a href=\"https://doi.org/10.1093/bioinformatics/btp352\">10.1093/bioinformatics/btp352</a></li>",
        "<li>Li H. (2012). seqtk: Toolkit for processing sequences in FASTA/Q formats. URL: <a href=\"https://github.com/lh3/seqtk\">https://github.com/lh3/seqtk</a></li>"
    ].join(' ').trim()

    def reference_text = [
        text_seq_qc,
        text_shortread_qc,
        params.perform_longread_qc ? text_longread_qc : "",
        params.perform_subsampling ? text_subsampling : "",
        text_shortreadmapping,
        text_longreadmapping,
        text_consensus,
        text_snps,
        text_statistics,
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
    ].join(' ').trim().replaceAll("[,|.] +\\.", ".")

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    // meta["tool_citations"] = ""
    // meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

