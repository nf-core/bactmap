process GET_GENOME_SIZE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(size_file)

    output:
    env 'genome_size', emit: ch_genome_size
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def cmd = '{sum+=$2} END {print sum}'
    """
    genome_size=`cat ${size_file} | awk '${cmd}'`

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
