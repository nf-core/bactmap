// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process VCF2PSEUDOGENOME {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? 'bioconda::pysam=0.16.0.1 conda-forge::biopython=1.78' :  null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0"
    }

    input:
    tuple val(meta), path(vcf)
    path reference

    output:
    tuple val(meta), path("${meta.id}.fas"), emit: pseudogenome

    script: // This script is bundled with the pipeline, in nf-core/bactmap/bin/
    def software = getSoftwareName(task.process)
    """
    vcf2pseudogenome.py  -r ${reference} -b ${vcf} -o ${meta.id}.fas
    """
}
