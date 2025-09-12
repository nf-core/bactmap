process READSTATS_PARSE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c0/c03d8fbb44376692e10bb9e56be2577fc35446c193637735de9fed182e6b58df/data' :
        'community.wave.seqera.io/library/pandas:2.3.1--139e2fa6c1f18206' }"

    input:
    path csv

    output:
    path "read_stats_summary.tsv", emit: tsv
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline in avantonder/bovisanalyzer/bin/
    def parser_version = '1.0'
    """
    read_stats_parser.py
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        read_stats_parser.py: ${parser_version}
    END_VERSIONS
    """
}
