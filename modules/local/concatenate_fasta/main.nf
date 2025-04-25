process CONCATENATE_FASTA {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'quay.io/biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def aligner_version = '1.0'
    """
    multi2single_sequence.py -r ${fasta} -o ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multi2single_sequence.py: ${aligner_version}
    END_VERSIONS
    """
}
