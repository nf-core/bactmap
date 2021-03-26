// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SNPSITES {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::snp-sites=2.5.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/snp-sites:2.5.1--hed695b0_0"
    } else {
        container "quay.io/biocontainers/snp-sites:2.5.1--hed695b0_0"
    }

    input:

    path alignment

    output:
    
    path "aligned_pseudogenomes.variants_only.fas", emit: variant_only_alignment
    path "*.sites.txt", emit: constant_sites
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    # output alignment with variant sites only
    snp-sites \\
        -o aligned_pseudogenomes.variants_only.fas \\
        $alignment
    # output constant sites
    snp-sites \\
        -C \\
        -o constant.sites.txt \\
        $alignment
    
    echo \$(snp-sites -V 2>&1) | sed 's/snp-sites //' > ${software}.version.txt
    """
}
