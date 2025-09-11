process READ_STATS {
    label 'process_low'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c0/c03d8fbb44376692e10bb9e56be2577fc35446c193637735de9fed182e6b58df/data' :
        'community.wave.seqera.io/library/pandas:2.3.1--139e2fa6c1f18206' }"

    input:
    tuple val(meta), path(json), path(json)

    output:
    tuple val(meta), path("*.read_stats.csv"), emit: csv
    path  "versions.yml",                      emit: versions

    script: // This script is bundled with the pipeline in avantonder/bovisanalyzer/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    def parser_version = '1.0'
    """
    read_stats.py
    mv read_stats.csv ${prefix}.read_stats.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        read_stats.py: ${parser_version}
    END_VERSIONS
    """
}
