// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RAPIDNJ {
    tag "$variant_alignment"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::rapidnj=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/rapidnj:2.3.2--hc9558a2_0"
    } else {
        container "quay.io/biocontainers/rapidnj:2.3.2--hc9558a2_0"
    }

    input:

    path variant_alignment

    output:
    path "*.sth", emit: stockholm_alignment
    path "*.tre", emit: phylogeny
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    python3.6 \\
        -c 'from Bio import SeqIO; SeqIO.convert("$variant_alignment", "fasta", "aligned_pseudogenome.variants_only.sth", "stockholm")'
    
    rapidnj \\
        aligned_pseudogenome.variants_only.sth \\
        -i sth \\
        -t d \\
        -b 1000 \\
        -n \\
        -c 1 \\
        -x rapidnj_phylogeny.tre
    
    //Doesn't appear to be a way of getting the version number
    
    echo 2.3.2 > ${software}.version.txt
    """
}
