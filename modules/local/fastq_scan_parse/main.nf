process FASTQSCANPARSE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c0/c03d8fbb44376692e10bb9e56be2577fc35446c193637735de9fed182e6b58df/data' :
        'community.wave.seqera.io/library/pandas:2.3.1--139e2fa6c1f18206' }"

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
