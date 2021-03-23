// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process IQTREE {
    tag "$variant_alignment"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::iqtree=2.0.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/iqtree:2.0.3--h176a8bc_1"
    } else {
        container "quay.io/biocontainers/iqtree:2.0.3--h176a8bc_1"
    }

    input:

    path variant_alignment
    path constant_sites

    output:
    path "*.treefile", emit: phylogeny
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    iqtree \\
        -fconst $(cat $constant_sites) \\ // I assume this'll work
        -s $variant_alignment \\
        -nt AUTO \\
        -ntmax $task.cpus \\
        -mem $task.memory \\
        -bb 1000 \\
        -m MFP
    
    // Not sure this is correct!
    
    echo (iqtree -version 2>&1) | sed 's/^.IQ-TREE multicore version //' > ${software}.version.txt
    """
}
