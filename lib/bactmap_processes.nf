// index reference sequence
process prepare_reference {
  tag {'prepare reference'}

  input:
  file ('reference.fasta')

  output:
  file 'reference.*' 

  """
  bwa index reference.fasta
  """
}

// get fastq files from the read archive at EBI
process fetch_from_ena {
    tag { sample_id }

    publishDir "${params.outdir}/fastqs",
        mode: 'copy',
        saveAs: { file -> file.split('\\/')[-1] }

    input:
    val sample_id

    output:
    set sample_id, file("${sample_id}/*.fastq.gz") into raw_fastqs

    """
    enaDataGet -f fastq -as /aspera.ini ${sample_id}
    """
}

// determine min read length to keep post trimming
process determine_min_read_length {
  tag { sample_id }

  input:
  set sample_id, file(file_pair)

  output:
  set sample_id, stdout

  script:
  """
  seqtk sample -s123 ${file_pair[0]} 1000 | printf "%.0f\n" \$(awk 'NR%4==2{sum+=length(\$0)}END{print sum/(NR/4)/3}')
  """
}

// Trim reads
process trim_reads {
  tag { sample_id }


  publishDir "${params.outdir}/trimmed_fastqs",
  mode: 'copy'

  input:
  set sample_id, min_read_length, file(reads)
  file('adapter_file.fas')

  output:
  set sample_id, file("${sample_id}_1.trimmed.fastq.gz"), file("${sample_id}_2.trimmed.fastq.gz")
  """
  trimmomatic PE -phred33 ${reads} ${sample_id}_1.trimmed.fastq.gz ${sample_id}_1_unpaired.trimmed.fastq.gz ${sample_id}_2.trimmed.fastq.gz ${sample_id}_2_unpaired.trimmed.fastq.gz ILLUMINACLIP:adapter_file.fas:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:${min_read_length}
    """
}

 // Genome Size Estimation
process genome_size_estimation {
  tag { sample_id }

  input:
  set sample_id, file(for_reads), file(rev_reads)

  output:
  set sample_id, file('mash_stats.out')

  """
  mash sketch -o /tmp/sketch_${sample_id}  -k 32 -m 3 -r ${for_reads}  2> mash_stats.out
  """
}

def find_genome_size(sample_id, mash_output) {
  m = mash_output =~ /Estimated genome size: (.+)/
  genome_size = Float.parseFloat(m[0][1]).toInteger()
  return [sample_id, genome_size]
}

// Estimate total number of reads
process count_number_of_reads {
  tag { sample_id }

  input:
  set sample_id, file(for_reads), file(rev_reads)

  output:
  set sample_id, file('seqtk_fqchk.out')

  """
  seqtk fqchk -q 25 ${for_reads} > seqtk_fqchk.out
  """
}

def find_total_number_of_reads(sample_id, seqtk_fqchk_ouput){
  m = seqtk_fqchk_ouput =~ /ALL\s+(\d+)\s/
  total_reads = m[0][1].toInteger() * 2 // the *2 is an estimate since number of reads >q25 in R2 may not be the same
  return [sample_id, total_reads]
}

// map reads
process map_reads {
  tag { sample_id }

  publishDir "${params.outdir}/sorted_bams",
    mode: 'copy'

  input:
  file ('reference.fasta')
  file reference_file
  set sample_id, file(for_reads), file(rev_reads), genome_size, read_count

  output:
  set sample_id, file("${sample_id}.sorted.bam")

  script:
  if (params.depth_cutoff  && read_count/genome_size > depth_cutoff.toInteger()){
    downsampling_factor = params.depth_cutoff.toInteger()/(read_count/genome_size)
    """
    mkdir downsampled_fastqs
    seqtk sample  ${for_reads} ${downsampling_factor} | gzip > downsampled_fastqs/${for_reads}
    seqtk sample  ${rev_reads} ${downsampling_factor} | gzip > downsampled_fastqs/${rev_reads}
    bwa mem reference.fasta downsampled_fastqs/${for_reads} downsampled_fastqs/${rev_reads} | samtools view -bS -F 4 - | samtools sort -O bam -o ${sample_id}.sorted.bam -
    """
  } else {
    """
    bwa mem reference.fasta ${for_reads} ${rev_reads} | samtools view -bS -F 4 - | samtools sort -O bam -o ${sample_id}.sorted.bam -
    """
  }
}

// call variants
process call_variants {
    tag { sample_id }

    input:
    file ('reference.fasta')
    set sample_id, file(sorted_bam)

    output:
    set sample_id, file("${sample_id}.bcf")

    script:
    """
    bcftools mpileup -a DP,AD,SP,ADF,ADR,INFO/AD,INFO/ADF,INFO/ADR -f reference.fasta ${sorted_bam} | bcftools call --ploidy 1 -m -Ob -o ${sample_id}.bcf
    """
}

// filter variants
process filter_variants {
    tag { sample_id }

    publishDir "${params.outdir}/filtered_bcfs",
      mode: 'copy'

    input:
    set sample_id, file(bcf_file)

    output:
    set sample_id, file("${sample_id}.filtered.bcf")

    script:
    """
    bcftools filter -s LowQual -e'${params.filtering_conditions}' ${bcf_file} -Ob -o ${sample_id}.filtered.bcf

    """
}

// produce pseudogenome
process create_pseudogenome {
  tag { sample_id }

  publishDir "${params.outdir}/pseudogenomes",
    mode: 'copy'

  input:
  file ('reference.fasta')
  set sample_id, file(filtered_bcf_file)

  output:
  file("${sample_id}.fas")

  script:
  """
  filtered_bcf_to_fasta.py  -r reference.fasta -b ${filtered_bcf_file} -o ${sample_id}.fas
  """
}

// create all sample pseudogenome alignment
process create_pseudogenome_alignment{
  tag { 'create pseudogenome alignment' }

  publishDir "${params.outdir}/pseudogenomes",
    mode: 'copy'

  input:
  file ('reference.fasta')
  file(pseudogenomes)

  output:
  file('aligned_pseudogenome.fas')
  file('low_quality_pseudogenomes.tsv')

  script:
  """
  touch low_quality_pseudogenomes.tsv
  for pseudogenome in ${pseudogenomes}
  do
    fraction_non_GATC_bases=`calculate_fraction_of_non_GATC_bases.py -f \$pseudogenome | tr -d '\\n'`
    if [[ 1 -eq "\$(echo "\$fraction_non_GATC_bases < ${params.non_GATC_bases_threshold}" | bc)" ]]; then
      cat \$pseudogenome >> aligned_pseudogenome.fas
    else
      echo "\$pseudogenome\t\$fraction_non_GATC_bases" >> low_quality_pseudogenomes.tsv
    fi
  done
  cat reference.fasta >> aligned_pseudogenome.fas
  """
}

//  remove non-informative positions with snp-sites

process create_variant_only_alignment{
  memory '15GB'
  cpus 4

  tag { 'create variant only pseudogenome alignment' }

  publishDir "${params.outdir}/pseudogenomes",
    mode: 'copy'

  input:
  file('aligned_pseudogenome')

  output:
  file('*.fas')

  script:
  if (params.remove_recombination){
    """
    run_gubbins.py --threads 4 -v -t hybrid aligned_pseudogenome
    snp-sites aligned_pseudogenome.filtered_polymorphic_sites.fasta -o aligned_pseudogenome.gubbins.variants_only.fas
    """

  } else {
    """
    snp-sites aligned_pseudogenome -o aligned_pseudogenome.variants_only.fas
    """
  }
}

// Build ML tree
process build_tree {
  memory { 15.GB * task.attempt }
  cpus 4

  tag {'build tree'}
  publishDir "${params.outdir}",
    mode: 'copy'

  input:
  file('aligned_pseudogenome.variants_only')

  output:
  file("*.treefile")
  file("*.contree")

  script:
  if (params.remove_recombination){
    """
    iqtree -s aligned_pseudogenome.variants_only -pre aligned_pseudogenome.gubbins.variants_only -m GTR+G -alrt 1000 -bb 1000 -nt AUTO -ntmax 4
    """
  } else {
    """
    iqtree -s aligned_pseudogenome.variants_only -m GTR+G -alrt 1000 -bb 1000 -nt AUTO -ntmax 4
    """
  }
}