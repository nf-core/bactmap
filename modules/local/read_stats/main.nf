process READ_STATS {
    label 'process_low'
    tag "$meta.id"
   
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92491bcdc9430c8a9feed2e23772ab75d9ff0671853f4eaaab30befdacf33d54/data' :
        'community.wave.seqera.io/library/pip_pandas:8dba217693f72600' }"
    
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