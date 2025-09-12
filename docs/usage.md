# nf-core/bactmap: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/bactmap/usage](https://nf-co.re/bactmap/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

**nf-core/bactmap** is a bioinformatics best-practice analysis pipeline for mapping short (Illumina) and long reads (Oxford Nanopore) from bacterial WGS to a reference sequence, creating filtered VCF files and making pseudogenomes based on high quality positions in the VCF files.

In addition to this page, you can find additional usage information on the following pages:

- [Tutorials](usage/tutorials.md)
- [FAQ and Troubleshooting](usage/faq-troubleshooting.md)

## General Usage

To run nf-core/bactmap, at a minimum two you require two inputs:

- a sequencing read samplesheet
- a reference genome

The samplesheet contains metadata and paths to the data of your input samples.

nf-core/bactmap includes optional pre-processing (adapter clipping, merge running etc.) or post-processing (visualisation) steps. These are opt in with a `--perform_<step>` flag. In some cases, the pre- and post-processing steps may also require additional files. Please check the parameters tab of this documentation for more information.

Please see the rest of this page for information about how to prepare input samplesheets and databases and how to run Nextflow pipelines. See the [parameters](https://nf-co.re/bactmap/parameters) documentation for more information about specific options the pipeline also offers.

## Samplesheet inputs

nf-core/bactmap can accept as input raw or preprocessed single- or paired-end short-read (e.g. Illumina) FASTQ files and long-read FASTQ files (e.g. Oxford Nanopore).

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 5 columns, and a header row as shown in the examples below.

This samplesheet is then specified on the command line as follows:

```console
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate different runs FASTQ files of the same sample before performing mapping and variant calling, when `--perform_runmerging` is supplied. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,run_accession,instrument_platform,fastq_1,fastq_2
2612,lane1,ILLUMINA,2612_lane1_R1.fq.gz,ILLUMINA,2612_lane1_R2.fq.gz
2612,lane2,ILLUMINA,2612_lane2_R1.fq.gz,ILLUMINA,2612_lane2_R2.fq.gz
2612,lane3,ILLUMINA,2612_lane3_R1.fq.gz,
```

::: info
Please note that the column name `run_accession` follows the definition of an ENA 'run'.
A 'run' corresponds to a single or paired-end set of demultiplexed FASTQs.
Given that demultiplexing of a given library happens per lane, each sequencing pair from each lane is a 'run'.
Therefore, for each sample, you may get multiple 'runs' consisting of _both_ lanes (of the same library) _and_ sequencing libraries.
Therefore ensure that each `run_accession` ID is unique, even if from the same sample!
:::

:::warning
Runs of the same sample sequenced on Illumina platforms with a combination of single and paired-end data will **not** be run-wise concatenated, unless pair-merging is specified. In the example above, `run3` will be profiled independently of `run1` and `run2` if pairs are not merged.
:::

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 5 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 5 samples, where `2612` has been sequenced twice.

```csv title="samplesheet.csv"
sample,run_accession,instrument_platform,fastq_1,fastq_2
2612,ERR5766176,ILLUMINA,/<path>/<to>/fastq/ERX5474932_ERR5766176_1.fastq.gz,/<path>/<to>/fastq/ERX5474932_ERR5766176_2.fastq.gz
2612,ERR5766180,ILLUMINA,/<path>/<to>/fastq/ERX5474936_ERR5766180_1.fastq.gz,
2613,ERR5766181,ILLUMINA,/<path>/<to>/fastq/ERX5474937_ERR5766181_1.fastq.gz,/<path>/<to>/fastq/ERX5474937_ERR5766181_2.fastq.gz
ERR3201952,ERR3201952,OXFORD_NANOPORE,/<path>/<to>/fastq/ERR3201952.fastq.gz,
```

:::warning
Input FASTQ _must_ be gzipped.
:::

:::warning
While one can include both short-read and long-read data in one run, we recommend that you split these across _two_ pipeline runs. This will make MultiQC run-reports more readable (due to run statistics having vary large number differences).
:::

| Column                | Description                                                                                                                                                                           |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`              | Unique sample name [required].                                                                                                                                                        |
| `run_accession`       | Run ID or name unique for each (pairs of) file(s) .Can also supply sample name again here, if only a single run was generated [required].                                             |
| `instrument_platform` | Sequencing platform reads generated on, selected from the EBI ENA [controlled vocabulary](https://www.ebi.ac.uk/ena/portal/api/controlledVocab?field=instrument_platform) [required]. |
| `fastq_1`             | Path or URL to sequencing reads or for Illumina R1 sequencing reads in FASTQ format. GZipped compressed files accepted. Can be left empty if data in FASTA is specified.              |
| `fastq_2`             | Path or URL to Illumina R2 sequencing reads in FASTQ format. GZipped compressed files accepted. Can be left empty if single end data.                                                 |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Reference genome

The reference genome must be provided in FASTA format. The pipeline will automatically index the reference genome using `bwa index` or `bowtie2 build` and `samtools faidx` if not already indexed.
The reference genome can be provided as a local file or a URL. The pipeline will automatically download the reference genome if a URL is provided.
The reference genome can be specified on the command line as follows:

```console
--fasta '[path to reference genome]'
```

The reference genome can be a single FASTA file or a multi-FASTA file. The reference genome can also be a gzipped FASTA file. The pipeline will automatically unzip the file if it is gzipped.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/bactmap --input ./samplesheet.csv --outdir ./results --fasta reference.fasta -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/bactmap -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
reference: 'reference.fasta'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Sequencing quality control

[`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. nf-core bactmap offers [`falco`](https://github.com/smithlabcode/falco) as an drop-in replacement, with supposedly better improvement particularly for long reads.

### Preprocessing Steps

nf-core/bactmap offers two main preprocessing steps for preprocessing raw sequencing reads:

- [**Read processing**](#read-processing): adapter clipping and pair-merging.
- [**Run merging**](#run-merging): concatenation of multiple FASTQ chunks/sequencing runs/libraries of a sample.

:::info
You can save the 'final' reads used for classification/profiling from any combination of these steps with `--save_analysis_ready_fastqs`.
:::

#### Read Processing

Raw sequencing read processing in the form of adapter clipping and paired-end read merging can be activated via the `--perform_shortread_qc` or `--perform_longread_qc` flags.

It is highly recommended to run this on raw reads to remove artifacts from sequencing that can cause false positive identification of taxa (e.g. contaminated reference genomes) and/or skews in taxonomic abundance profiles. If you have public data, normally these should have been corrected for, however you should still check that these steps have indeed been already performed.

There are currently two options for short-read preprocessing: [`fastp`](https://github.com/OpenGene/fastp) or [`adapterremoval`](https://github.com/MikkelSchubert/adapterremoval).

For adapter clipping, you can either rely on the tool's default adapter sequences, or supply your own adapters (`--shortread_qc_adapter1` and `--shortread_qc_adapter2`)
By default, paired-end merging is not activated. If paired-end merging is activated you can also specify whether to include unmerged reads in the reads sent for mapping/variant calling (`--shortread_qc_mergepairs` and `--shortread_qc_includeunmerged`).
You can also turn off clipping and only perform paired-end merging, if requested. This can be useful when processing data downloaded from the ENA, SRA, or DDBJ (`--shortread_qc_skipadaptertrim`).
Both tools support length filtering of reads and can be tuned with `--shortread_qc_minlength`.

There are currently two options for long-read Oxford Nanopore processing: [`porechop`](https://github.com/rrwick/Porechop), [`porechop_abi`](https://github.com/bonsai-team/Porechop_ABI).

For both short-read and long-read preprocessing, you can optionally save the resulting processed reads with `--save_preprocessed_reads`.

#### Run Merging

For samples that may have been sequenced over multiple runs, or for FASTQ files split into multiple chunks, you can activate the ability to merge across all runs or chunks with `--perform_runmerging`.

For more information how to set up your input samplesheet, see [Multiple runs of the same sample](#multiple-runs-of-the-same-sample).

Activating this functionality will concatenate the FASTQ files with the same sample name _after_ the optional preprocessing steps and _before_ mapping/variant calling. Note that libraries with runs of different pairing types will **not** be merged and this will be indicated on output files with a `_se` or `_pe` suffix to the sample name accordingly.

You can optionally save the FASTQ output of the run merging with the `--save_runmerged_reads`.

### Subsampling reads

Some sequencing runs may be too large to process in a reasonable time. In these cases, you can use the `--perform_subsampling` flag to randomly subsample your reads to a specified depth of coverage per sample. This is done before the mapping steps. By default, this step is activated for all samples.

### Read mapping

The nf-core/bactmap pipeline provides two strategies to map reads to a reference genome:

- **Short-read mapping**: `bwa-mem2` or `bowtie2`
- **Long-read mapping**: `minimap2`

By default, the pipeline will use `bowtie2` for short-read mapping. You can change this with the `shortread_mapping_tool` parameter.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/bactmap
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/bactmap releases page](https://github.com/nf-core/bactmap/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
