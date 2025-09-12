# nf-core/bactmap: Tutorials

This page provides a range of tutorials to help give you a bit more guidance on how to set up nf-core/bactmap runs yourself.

## Simple Tutorial

In this tutorial we will run you through a simple set up of a small nf-core/bactmap run.
It assumes that you have basic knowledge of read mapping input and output files.

### Preparation

#### Hardware

The datasets used should be small enough to run on your own laptop or a single server node.

If you wish to use a HPC cluster or cloud, and don't wish to use an 'interactive' session submitted to your scheduler, please see the [nf-core documentation](https://nf-co.re/docs/usage/configuration#introduction) on how to make a relevant config file.

You will need internet access and at least XXX GB of hardrive space.

#### Software

The tutorial assumes you are on a Unix based operating system, and have already installed Nextflow as well a software environment system such as [Conda](https://docs.conda.io/en/latest/miniconda.html), [Docker](https://www.docker.com/), or [Singularity/Apptainer](https://apptainer.org/).
The tutorial will use Docker, however you can simply replace references to `docker` with `conda`, `singularity`, or `apptainer` accordingly.

#### Data

First we will make a directory to run the whole tutorial in.

```bash
mkdir bactmap-tutorial
cd bactmap-tutorial/
```

We will use very small short-read FASTQ files used for testing. You can download these files, along with the reference genome, with the following commands:

```bash
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/prokaryotes/bacteroides_fragilis/illumina/fastq/test1_1.fastq.gz
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/prokaryotes/bacteroides_fragilis/illumina/fastq/test1_2.fastq.gz
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/prokaryotes/bacteroides_fragilis/illumina/fastq/test2_1.fastq.gz
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/prokaryotes/bacteroides_fragilis/illumina/fastq/test2_2.fastq.gz
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/prokaryotes/bacteroides_fragilis/nanopore/fastq/test.fastq.gz
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/prokaryotes/bacteroides_fragilis/reference/genome.fna.gz
```

### Preparing Input

#### Sample sheet

You provide the sequencing data FASTQ files to nf-core/bactmap via a input 'sample sheet' `.csv` file.
This is a 5 column table, that includes sample and library names, instrument platform, and paths to the sequencing data.

Open a text editor, and create a file called `samplesheet.csv`.
Copy and paste the following lines into the file and save it.

```csv title="samplesheet.csv"
sample,run_accession,instrument_platform,fastq_1,fastq_2
01,test_1,ILLUMINA,test1_1.fastq.gz,test1_2.fastq.gz
02,test_2,ILLUMINA,test2_1.fastq.gz,test2_2.fastq.gz
03,test_3,OXFORD_NANOPORE,test.fastq.gz,
```

Here we have specified three samples, two were sequenced on the Illumina platform, one on the Oxford Nanopore, and the paths to the FASTQ files.
If you had placed your FASTQ files elsewhere, you would give the full path (i.e., with relevant directories) to the `fastq_1` and `fastq_2` columns.

### Running the pipeline

Now that we have the sequencing reads (in FASTQ format) and a reference genome (FASTA, optionally gzipped), we can now run them with the pipeline. The following command will perform short read quality control, short- and long-read mapping, variant calling, create pseudogenomes and pseudogenome alignment, and finally read and mapping statistics.

```bash
nextflow run nf-core/bactmap -r 2.0.0 -profile docker \
--input samplesheet.csv --fasta genome.fna.gz --outdir ./results \
--perform_subsampling false
```

:::info
With all Docker containers pre-downloaded, this run took X minutes and X seconds on a laptop running XXX with XX GB RAM and XX CPUs.
If you are running nf-core/bactmap for the first time, expect this command to take longer as Nextflow will have to download each software container for each step of the pipeline.
:::

To break down each line of the command:

- Tell Nextflow to run nf-core/bactmap with the particular version and using the Docker container system
- Specify the input and outputs, i.e., paths to the `samplesheet.csv`, `genome.fna.gz`, and directory where to save the results
- Turn off subsampling, which is a step that randomly selects a subset of reads from the FASTQ files. This is useful for testing the pipeline on small datasets, but not needed for real data.
- (Optional) provide a _cap_ to the maximum amount of resources each step/job of the pipeline can use

:::warning
The pipeline runs occasionally fail due to a particular step of the pipeline requesting more resources than you have on your system. To avoid these failures, you can tell Nextflow to set a cap pipeline-step resource requests against a list called `resourceLimits` specified in Nextflow config file. These should represent the maximum possible resources of a machine or node. To learn how to increase computational resource to the pipeline, see the central [nf-core documentation](https://nf-co.re/docs/usage/configuration).
:::

### Output

In the resulting directory `results/` you will find a range of directories.

<!--TODO: add a diagram of the output structure -->

```tree
results/
├── alignpseudogenomes
├── bcftools
├── bowtie2
├── fastp
├── fastqc
├── minimap2
├── multiqc
├── nanoq
├── pipeline_info
├── porechop
├── pseudogenomes
├── rasusa
├── read_stats
├── samtools
├── snp-sites
```

<!--TODO: add a diagram of the output structure -->

To follow the same order as the command construction above

- Pipeline run report is found in `multiqc/` and resource statistics in `pipeline_info`
- Short-read QC results are found in `fastqc/` and `fastp/`
- Short-read mapping results are found in `bowtie2/`
- Long-read QC results are found in `fastqc/`, `porechop/` and `nanoq/`
- Read merging results are found in `run_merging/`
- Long-read mapping results are found in `minimap2/`
- Variant calling results are found in `samtools/`
- Variant calling statistics are found in `samtools/` and `bcftools/`
- Consensus pseudogenome results are found in `pseudogenomes/`
- Pseudogenome alignment results are found in `alignpseudogenomes/`
- Variant sites are found in `snp-sites/`

:::info
For read-preprocessing steps, only log files are stored in the `results/` directories by default. Refer to the parameters tab of the [nf-core/bactmap documentation](https://nf-co.re/bactmap/) for more options.
:::

The general 'workflow' of going through the results will typically be reviewing the `multiqc/multiqc_report.html` file to get general statistics of the entire run, particularly of the preprocessing.

Detailed descriptions of all results files can be found in the output tab of the [nf-core/bactmap documentation](https://nf-co.re/bactmap/).

### Clean up

Once you have completed the tutorial, you can run the following command to delete all downloaded and output files.

```bash
rm -r bactmap-tutorial/
```

:::warning
Don't forget to change out of the directory above before trying to delete it!
:::
