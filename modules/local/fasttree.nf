// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FASTTREE {
    tag "$variant_alignment"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::fasttree=2.1.10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fasttree:2.1.10--h516909a_4"
    } else {
        container "quay.io/biocontainers/fasttree:2.1.10--h516909a_4"
    }

    input:
    
    path variant_alignment

    output:
    path "*.tre", emit: phylogeny
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    //Does the biocontainer come with both fasttree and fasttreeMP?
    
    FastTreeMP \\
        -gtr \\
        -gamma \\
        -fastest \\
        -nt $variant_alignment \\
        2>&1 fasttree_phylogeny.tre
    
    //This needs to be checked
    
    echo \$(iqtree 2>&1) | sed 's/^.Usage for FastTree version //' > ${software}.version.txt
    """
}
