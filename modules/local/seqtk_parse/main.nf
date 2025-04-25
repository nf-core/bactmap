process SEQTK_PARSE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92491bcdc9430c8a9feed2e23772ab75d9ff0671853f4eaaab30befdacf33d54/data' :
        'community.wave.seqera.io/library/pip_pandas:8dba217693f72600' }"

    input:
    path tsv

    output:
    path "mapping_summary.tsv", emit: tsv
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline in avantonder/bactmap/bin/
    def parser_version = '1.0'
    """
    seqtk_parser.py
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk_parser.py: ${parser_version}
    END_VERSIONS 
    """
}