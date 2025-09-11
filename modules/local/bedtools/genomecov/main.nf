process BEDTOOLS_GENOMECOV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/6397750e9730a3fbcc5b4c43f14bd141c64c723fd7dad80e47921a68a7c3cd21/data':
        'community.wave.seqera.io/library/bedtools_coreutils:a623c13f66d5262b' }"

    // simpler input, having removed the sorting option
    input:
    tuple val(meta), path(intervals), val(scale)

    output:
    tuple val(meta), path("*.bed"), emit: genomecov
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args  ?: ''
    def args_list = args.tokenize()
    args += (scale > 0 && scale != 1) ? " -scale $scale" : ""
    if (!args_list.contains('-bg') && (scale > 0 && scale != 1)) {
        args += " -bg"
    }
    // Sorts output file by chromosome and position using additional options for performance and consistency
    // See https://www.biostars.org/p/66927/ for further details
    def buffer   = task.memory ? "--buffer-size=${task.memory.toGiga().intdiv(2)}G" : ''

    def prefix = task.ext.prefix ?: "${meta.id}"
    def cmd = '$4 <'
    // hard-coded for bcftools_consensus subworkflow. mainly the problem was the awk pipe at the end, we decided to go with a local version of genomecov.
    // this is a simpler version than the nf-core version, it expects a bam as interval file, and it does not expect sorting either. we did add a threshold parameter in the nextflow.config in case the user wants to specify the threshold for low coverage
    """
    bedtools \\
        genomecov \\
        -ibam $intervals \\
        $args \\
        | awk '${cmd}'$params.genomecov_threshold \\
        > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
