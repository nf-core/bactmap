process ALIGNPSEUDOGENOMES {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'quay.io/biocontainers/biopython:1.78' }"
    
    input:
    path pseudogenomes
    path reference

    output:
    tuple env(NUM_ALIGNMENT_GENOMES), path("aligned_pseudogenomes.fas"), emit: aligned_pseudogenomes
    path "low_quality_pseudogenomes.tsv",                                emit: low_quality_metrics
    path  "versions.yml",                                                emit: versions

    script: // This script is bundled with the pipeline, in avantonder/mapbac/bin/
    def aligner_version = '1.0'
    """
    touch low_quality_pseudogenomes.tsv
    touch aligned_pseudogenomes.fas
    for pseudogenome in ${pseudogenomes}
    do
        fraction_non_GATC_bases=\$(calculate_fraction_of_non_GATC_bases.py -f \$pseudogenome | tr -d '\\n')
        if awk 'BEGIN { exit !(\$fraction_non_GATC_bases < ${params.non_GATC_threshold}) }'; then
            cat \$pseudogenome >> aligned_pseudogenomes.fas
        else
            echo "\$pseudogenome\t\$fraction_non_GATC_bases" >> low_quality_pseudogenomes.tsv
        fi
    done
    reference2single_sequence.py -r ${reference} -o final_reference.fas
    cat final_reference.fas >> aligned_pseudogenomes.fas

    NUM_ALIGNMENT_GENOMES=\$(grep -c ">" aligned_pseudogenomes.fas)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reference2single_sequence.py: ${aligner_version}
    END_VERSIONS
    """
}