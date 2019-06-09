# nf-core/bactmap
**A mapping-based pipeline for creating a phylogeny from bacterial whole genome sequences**

[![Build Status](https://travis-ci.org/nf-core/bactmap.svg?branch=master)](https://travis-ci.org/nf-core/bactmap)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.05.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/bactmap.svg)](https://hub.docker.com/r/nfcore/bactmap)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.


### Documentation
The nf-core/bactmap pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

### Pipeline description
This pipeline maps paired end short reads to a bacterial fasta reference sequence, calls qnd filters variants, produces a whole genome alignment from pseudogenomes derived from the variants and finally produces a robust maximum likelihood phylogentic tree.


#### Pipeline steps
The steps are:

1. Index a reference sequnce using bwa (the reference sequence must only contain the chromosome and no additional sequences such as plasmids).
2. (Optional) Fetch reads from the ENA
3. Trim reads using trimmomatic (dynamic MIN_LEN based on 33% of the read length)
4. Count number of reads and estimate genome size using Mash
5. Downsample reads if the `--depth_cutoff` argument was specified
6. Map reads to the specified reference genome with bwa mem
7. Call variants with samtools
8. Filter variants to flag low quality SNPs
9. Produce a pseudogenome based on the variants called. Missing positions are encoded as `-` characters and low quality positions as `N`
10. All pseudogenomes are concatenanted to make a whole genome alignment
11. (Optional) Recombination is removed from the alignment using gubbins
12. Invariant sites are removed using snp-sites
13. (Optional) Maximum likelihood tree generated using IQ-TREE


A sumary of this process is shown below in the diagram that was generated when running Nextflow using the -with-dag command

![workflow diagram](dag.png)

#### Pipeline outputs
These will be found in the directory specified by the `--output_dir` argument

  - (Optional) If accession numbers were used as the input source a directory called `fastqs` will contain the fastq file pairs for each accession number
  - A directory called `trimmed_fastqs` containing the reads after trimminb with TRIMMOMATIC
  - A directory called `sorted_bams` containing the alignmed sam files after mapping with bwa mem, conversion to bam and sorting
  - A directory called `filtered_bcfs` containing binary vcf files after filtering to flag low quality positions with LowQual in the FILTER column
  - A directory called `pseudogenomes` containing 
    - the pseudogenome from each sample
    - a whole genome alignment named `aligned_pseudogenome.fas` containing the concatenated sample pseudogenomes and the refrerence genome
    - a variant only alignment named `aligned_pseudogenome.variants_only.fas` with the invariant sites removed from `aligned_pseudogenome.fas` using snp-sites. If recombination removal was specified, the file will be named `aligned_pseudogenome.gubbins.variants_only.fas` with gubbins having been applied prior to invariant site removal.
  - Two newick tree files
   -  `aligned_pseudogenome.gubbins.variants_only.contree` If tree generation was specified, this file containing the consensus tree from IQ_TREE will be produced. The tree will possess assigned branch supports where branch lengths are optimized on the original alignment. If recombination removal was not specified the file will be named `aligned_pseudogenome.variants_only.contree`
   -  `aligned_pseudogenome.gubbins.variants_only.treefile` The original IQ-TREE maximum likelihood tree without branch supports. If recombination removal was not specified the file will be named `aligned_pseudogenome.variants_only.treefile`





### Credits
nf-core/bactmap was originally written by Anthony Underwood.

#### Software used within the workflow
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

