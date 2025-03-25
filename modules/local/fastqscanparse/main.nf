process FASTQSCANPARSE {
    label 'process_low'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-9adca5a7d3b24119897cfc20386da6c7fa47bdab:77c1885b47edc369aceb4cccf161a549bdac3d4b-0' :
        'quay.io/biocontainers/mulled-v2-9adca5a7d3b24119897cfc20386da6c7fa47bdab:77c1885b47edc369aceb4cccf161a549bdac3d4b-0' }"
    
    input:
    path json

    output:
    path("*.tsv"),         emit: tsv
    path  "versions.yml" , emit: versions
    
    script: // This script is bundled with the pipeline in avantonder/bovisanalyzer/bin/
    def prefix = task.ext.prefix ?: ''
    def parser_version = '1.0'
    """
    fastqscan_parser.py
    mv fastq-scan_summary.tsv ${prefix}_fastq-scan_summary.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqscan_parser.py: ${parser_version}
    END_VERSIONS 
    """
}