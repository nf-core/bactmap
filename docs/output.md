# nf-core/bactmap: Output

This document describes the output produced by the pipeline. 

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:


* [Fetch from ENA](#fetch-from-ena) (Optional) Fetch reads from the ENA
* [Trim Reads](#trim-reads) Read trimmimg using trimmomatic 
* [Estimate genome size](#genome-size-estimation)
* [Downsample reads](#downsample-reads) 
* [Map reads](#map-reads)
* [Call variants](#call-variants)
* [Filter variant](#filter-variants)
* [Pseudogenome creation](#pseudogenome-creation)
* [Pseudogenome alignment creation](#pseudogenome-alignment-creation)
* [Recombination removal](#recombination-removal)(Optional)
* [Invariant site removal](#invariant-site-removal)
* [Phylogenetic tree creation](#phylogenetic-tree-removal) (Optional) 



## Fetch from ENA
This process will fetch reads from the ENA archive using the `enaDataGet` tool from [ENA Browser Tools](https://github.com/enasequence/enaBrowserTools)

## Trim Reads
Trim with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) reads based on the parameters ILLUMINACLIP:adapter_file.fas:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 and with MIN_LEN dynamically determined based on 30% of the read length
**Output directory: `<OUTPUT DIR>/trimmed_fastqs`**
Fastq files post trimming will be written here

## Genome Size Estimation
Estimate the size of the genome using [Mash](https://mash.readthedocs.io/en/latest/)

## Downsample reads
If the `--depth_cutoff` parameter is specified then reads will be downsampled using [seqtk](https://github.com/lh3/seqtk) to the specified depth


## Map reads
The reads will be mapped to the specified reference genome using [bwa mem](https://github.com/lh3/bwa)
**Output directory: `<OUTPUT DIR>/sorted_bams`**
Sorted bam files will be written here

## Call variants
Variants will be called using [samtools](http://www.htslib.org/doc/samtools.html)

## Filter variants
Variants will be filtered using [bcftools](http://www.htslib.org/doc/bcftools.html) in order to flag low quality SNPs using the default filter of `%QUAL<25 || FORMAT/DP<10 || MAX(FORMAT/ADF)<5 || MAX(FORMAT/ADR)<5 || MAX(FORMAT/AD)/SUM(FORMAT/DP)<0.9 || MQ<30 || MQ0F>0.1` 
**Output directory: `<OUTPUT DIR>/filtered_bcfs`**
Filtered vcf files will be written here
## Pseudogenome creation
A pseudogenome based on the variants called is created where missing positions are encoded as `-` characters and low quality positions as `N`. All other positions either match the reference or are encoded as a SNV of either G,A,T or C. The script [filtered_bcf_to_fasta.py](bin/filtered_bcf_to_fasta.py) is used.
**Output directory: `<OUTPUT DIR>/pseudogenomes`**
A pseudogenome for each sample will be written here

## Pseudogenome alignment creation
The pseudogenomes from the previous step are concatenanted to make a whole genome alignment
**Output directory: `<OUTPUT DIR>/pseudogenomes`**
The multi-sample pseudogenome alignment will be written here

## Recombination removal
Recombination is removed from the alignment using [gubbins](https://github.com/sanger-pathogens/gubbins/)

## Invariant sites
Invariant sites are removed using [snp-sites](https://github.com/sanger-pathogens/snp-sites)

## Phylogenetic tree creation
A Maximum likelihood tree is generated using [IQ-TREE](http://www.iqtree.org)
**Output directory: `<OUTPUT DIR>`**
The consensus tree `aligned_pseudogenome.variants_only.contree` including bootstrap values will be written here

# Software used within the pipeline
  - [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) A flexible read trimming tool for Illumina NGS data.
  - [mash](https://mash.readthedocs.io/en/latest/) Fast genome and metagenome distance estimation using MinHash.
  - [seqtk](https://github.com/lh3/seqtk) A fast and lightweight tool for processing sequences in the FASTA or FASTQ format.
  - [bwa mem](https://github.com/lh3/bwa) Burrow-Wheeler Aligner for short-read alignment
  - [samtools](http://www.htslib.org/doc/samtools.html) Utilities for the Sequence Alignment/Map (SAM) format
  - [bcftools](http://www.htslib.org/doc/bcftools.html) Utilities for variant calling and manipulating VCFs and BCFs
  - [filtered_bcf_to_fasta.py](bin/filtered_bcf_to_fasta.py) Python utility to create a pseudogenome from a bcf file where each position in the reference genome is included
  - [gubbins](https://github.com/sanger-pathogens/gubbins/) Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences
  - [snp-sites](https://github.com/sanger-pathogens/snp-sites) Finds SNP sites from a multi-FASTA alignment file 
  - [IQ-TREE](http://www.iqtree.org) Efficient software for phylogenomic inference

