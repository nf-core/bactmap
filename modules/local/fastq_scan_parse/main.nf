process FASTQSCANPARSE {
    label 'process_low'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92491bcdc9430c8a9feed2e23772ab75d9ff0671853f4eaaab30befdacf33d54/data' :
        'community.wave.seqera.io/library/pip_pandas:8dba217693f72600' }"
    
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