// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RAXMLNG {
    tag "$alignment"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::raxml-ng:1.0.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/raxml-ng:1.0.2--h7447c1b_0"
    } else {
        container "quay.io/biocontainers/raxml-ng:1.0.2--h7447c1b_0"
    }

    input:

    path alignment

    output:
    path "*.log", emit: raxml_log
    path "*.rba", emit: binary_variant_alignment
    path "*.rba", emit: phylogeny
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    raxml-ng --parse \\
        --msa $alignment \\
        --model GTR+G \\



    
    raxml-ng --all
    
    
    raxml-ng \\
        --msa $variant_alignment \\
        --model GTR+G \\
        --threads $task.cpus
        --bs-trees 1000

    echo \$(raxml-ng --version 2>&1) | sed 's/^.*RAxML-NG v. //; s/released.*\$//' > ${software}.version.txt
    """
}
