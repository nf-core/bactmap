process SEQTK_PARSE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-9adca5a7d3b24119897cfc20386da6c7fa47bdab:77c1885b47edc369aceb4cccf161a549bdac3d4b-0' :
        'quay.io/biocontainers/mulled-v2-9adca5a7d3b24119897cfc20386da6c7fa47bdab:77c1885b47edc369aceb4cccf161a549bdac3d4b-0' }"

    input:
    path tsv

    output:
    path "mapping_summary.tsv", emit: tsv
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    script: // This script is bundled with the pipeline in avantonder/mapbac/bin/
    def parser_version = '1.0'
    """
    seqtk_parser.py
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk_parser.py: ${parser_version}
    END_VERSIONS 
    """
}