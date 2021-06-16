# nf-core/bactmap: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/bactmap/usage](https://nf-co.re/bactmap/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

This pipeline maps short reads (usually Illumina) to a bacterial reference (usually of the same species and the closest high quality genome available). For this purpose you will require a set of paired read files (single read file version coming in future version) and a reference genome in fasta format. The read file pairs are specified in a sample sheet. Please note although the pipeline can handle multiple contigs within the reference sequence, it is recommended that plasmid records are removed leaving only the chromosomal records (usually one chromosome in most bacteria) since plasmids are often acquired horizontally and are evolving at a different rate to the chromosome.

## Samplesheet input

You will need to create a samplesheet file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the example below.

```bash
--input '[path to samplesheet file]'
```

```bash
sample,fastq_1,fastq_2
G18582004,fastqs/G18582004_1.fastq.gz,fastqs/G18582004_2.fastq.gz
G18756254,fastqs/G18756254_1.fastq.gz,fastqs/G18756254_2.fastq.gz
G18582006,fastqs/G18582006_1.fastq.gz,fastqs/G18582006_2.fastq.gz
```

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/bactmap --input samplesheet.csv --reference chromosome.fasta -profile docker
```

Please note that although the reference is named `chromosome.fasta` in this example command it could be named anything BUT it is recommended that the reference only include chromosomal sequence (see note in introduction).

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Optional parameters

By default the pipeline does **not** perform the following steps but they can be enabled by adding the appropriate parameter to the command line:

* trim reads `--trim`
* remove recombiination using gubbins `--remove_recombination`
* build a RapidNJ tree `--rapidnj`
* build a FastTree tree `--fasttree`
* build an IQ-TREE tree `--iqtree`
* build a RAxML-NG tree `--raxmlng`

By default the pipeline subsamples reads to an estimated depth of 100x.
This can be turned off using the `--subsampling_off` parameter.
The desired coverage depth can be changed using the `subsampling_depth_cutoff` parameter

The steps in the pipeline have sensible defaults but complete control of the arguments passed to the software tools can be achieved by overriding the software arguments found in the [`modules.config`](../conf/modules.config) file with a custom config. For example the default args for IQ-TREE which specify using ModelFinder to find the best fit model could be overridden to specify a specific model. In the  [`modules.config`](../conf/modules.config) file these are specified as:

```bash
    'iqtree' {
        args = '-alrt 1000 -B 1000 -m MFP -czb'
        publish_dir = 'iqtree'
    }
```

These could be overridden by specifying a config file by adding an argument to the command line such as t `-c user.config` ([user.config](../assets/user.config)).
Example contents:

```bash
params {
    modules {
        'iqtree' {
            args = '-m GTR+G -czb'
        }
    }
}
```

This will specify that IQ-TREE should no longer provide SH-aLRT and the ultrafast bootstrap branch support values and will use the GTR+G model. The other options (`build` and `publish_dir`) will remain the same. Therefore to build an IQ-TREE, this step would either need to be turned on adding the `--iqtree` parameter to the Nextflow command line or by adding `iqtree = true` to the user config file within the params block.

The steps are described in more detail below along with their default parameters

### Comprehensive description the steps

1. **Reference sequence indexing**: reference sequence is indexed using [`bwa index`](https://github.com/lh3/bwa)
2. **Read trimming** (Optional if `params.trim` is set): reads are trimmed using [`fastp`](https://github.com/OpenGene/fastp). The default process configuration is found in the module.config and can be overridden as described above:

    ```bash
    'fastp' {
        args              = '--cut_front --cut_tail --trim_poly_x --cut_mean_quality 30 --qualified_quality_phred 30 --unqualified_percent_limit 10 --length_required 50'
        adapter_fasta     = "${baseDir}/assets/adapters.fas"
        publish_files     = ['json':'', 'html':'', 'log': 'log']
    }
    ```

    Please note that the default adapters are found in the [`adapters.fas`](../assets/adapters.fas) file. These can be supplemented by specifying the path to a different file and updating the `adapter_fasta` value in a config file specified using `-c` as described above.

3. **Subsample reads** (Optionally if `params.subsampling_off` is **not** set): reads are subsampled based on estimated genome size using [`mash sketch`](https://github.com/marbl/Mash) and the required depth of coverage as set by the `--subsampling_depth_cutoff` parameter using [`rasusa`](https://github.com/mbhall88/rasusa). The default process configurations are found in the module.config and can be overridden as described above.

    ```bash
        'mash_sketch' {
            args = '-k 32 -m 3'
        }
        'rasusa' {
            args = '--seed 23032021'
        }
    ```

4. **Map reads**:  reads are mapped to the indexed reference genome using [`bwa mem`](https://github.com/lh3/bwa) to produce a bam file.
The default process configuration is found in the module.config and can be overridden as described above:

    ```bash
        'bwa_mem' {
            args = ''
            args2 = '-F 4' // samtools view options discarding unmapped reads
            publish_files = false
        }
    ```

5. **Sort reads**: bam files are sorted using [`samtools`](http://www.htslib.org/doc/samtools.html)
6. **Call variants**: variants are called using [`bcftools mpileup`](http://samtools.github.io/bcftools/bcftools.html). A minimum base quality of 20 is used for pre-filtering and the following fields are included in the resulting VCF file `FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR`. A haploid and multiallelic model is assumed. These defaults are found in the module.config and can be overridden as described above:

    ```bash
        'bcftools_mpileup' {
            args          = '--min-BQ 20 --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR'
            args2         = '--ploidy 1 --multiallelic-caller'
            args3         = ''
            publish_files = false
            publish_dir   = 'variants'
        }
    ```

7. **Filter variants**: variants in the VCF file are filtered using [`bcftools filter`](http://samtools.github.io/bcftools/bcftools.html). The default process configuration is found in the module.config and can be overridden as described above:

    ```bash
        'bcftools_filter' {
            args = '--soft-filter LowQual --exclude "%QUAL<25 || FORMAT/DP<10 || MAX(FORMAT/ADF)<2 || MAX(FORMAT/ADR)<2 || MAX(FORMAT/AD)/SUM(FORMAT/DP)<0.9 || MQ<30 || MQ0F>0.1" --output-type z'
            suffix = '.filtered'
            publish_dir   = 'variants'
        }
    ```

    After this process a filtered VCF file will be produced that has a row for each position in the reference genome. Each of these will either have a value in the FILTER column of either `PASS` or `LowQual`.

8. **Create pseudogenome**: the filtered VCF is used to create a pseudogenome based on the reference genome using the [`vcf2pseudogenome.py script`](https://github.com/nf-core/bactmap/blob/dev/bin/vcf2pseudogenome.py). The base in a sample at a position where the VCF file row that has `PASS` in FILTER will be either ref or alt and the appropriate base will be encoded at that position. The base in a sample at a position where the VCF file row that has `LowQual` in FILTER is uncertain will be encoded as a `N` character. Missing data will be encoded as a `-` character.
9. **Align pseudogenomes** All samples pseudogenomes and the original reference sequence are concatenated together to produce a flush alignment where all the sequence for all samples at all positions in the original reference sequence will be one of `{G,A,T,C,N,-}`. Only those sequences that are high quality (based on the number of non GATC bases will be included.
The threshold for this is set in the default process configuration found in the module.config and can be overridden as described above.

    ```bash
    'alignpseudogenomes' {
        non_GATC_threshold = 0.5
        publish_dir = 'pseudogenomes'
    }
    ```

    This alignment can be used for other downstream processes such as phylogeny generation.
10. **Determine number of constant sites** Optionally this alignment can be processed to produce a phylogeny as part of this pipeline. The number of constant sites in the alignment will be determined using [`snp-sites`](https://github.com/sanger-pathogens/snp-sites).
11. **Remove recombination** (Optionally if `params.remove_recombination` is set): remove regions likely to have been acquired by horizontal transfer and recombination and therefore perturb the true phylogeny using [`gubbins`](https://sanger-pathogens.github.io/gubbins/). This should only be run on sets of samples that are closely related and not for example on a set of samples that have diversity spanning that of the entire species.
The default process configuration is found in the module.config and can be overridden as described above.

    ```bash
    'alignpseudogenomes' {
        non_GATC_threshold = 0.5
        publish_dir = 'pseudogenomes'
    }
    ```

12. **Build tree(s)**: Depending on the params set, run 0 - 4 tree building algorithms.
    * `--rapidnj` Build a neighbour-joining pylogeny using [`rapidnj`](https://birc.au.dk/software/rapidnj)
    The default process configuration is found in the module.config and can be overridden as described above.

    ```bash
        'rapidnj' {
            args = '-t d -b 1000 -n'
            publish_dir = 'rapidnj'
        }
    ```

    * `--fasttree` Build an approximately-maximum-likelihood phylogeny using [`FastTree`](http://www.microbesonline.org/fasttree)
    The default process configuration is found in the module.config and can be overridden as described above.

        ```bash
        'fasttree' {
            args = '-gtr -gamma -fastest'
            publish_dir = 'fasttree'
        }
        ```

    * `--iqtree` Build a maximum-likelihood phylogeny using [`IQ-TREE`](http://www.iqtree.org)
    The default process configuration is found in the module.config and can be overridden as described above.

        ```bash
        'iqtree' {
            args = '-alrt 1000 -B 1000 -m MFP -czb'
            publish_dir = 'iqtree'
        }
        ```

    * `--raxmlng` Build a maximum-likelihood phylogeny using [`RAxML Next Generation`](https://github.com/amkozlov/raxml-ng)
    The default process configuration is found in the module.config and can be overridden as described above.

        ```bash
        'raxmlng' {
            args = '--model GTR+G --bs-trees 1000'
            publish_dir = 'raxmlng'
        }
        ```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/bactmap
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/bactmap releases page](https://github.com/nf-core/bactmap/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
    * A generic configuration profile to be used with [Docker](https://docker.com/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
* `podman`
    * A generic configuration profile to be used with [Podman](https://podman.io/)
* `shifter`
    * A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
* `charliecloud`
    * A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
* `conda`
    * A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if a pipeline , `nf-core/rnaseq` is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue. For an in-depth explanation of how to fix these types of errors please see [here](fixing_resource_errors.md)

### Tool-specific options

For the ultimate flexibility, we have implemented and are using Nextflow DSL2 modules in a way where it is possible for both developers and users to change tool-specific command-line arguments (e.g. providing an additional command-line argument to a process) as well as publishing options (e.g. saving files produced by the process that aren't saved by default by the pipeline). In the majority of instances, as a user you won't have to change the default options set by the pipeline developer(s), however, there may be edge cases where creating a simple custom config file can improve the behaviour of the pipeline if for example it is failing due to a weird error that requires setting a tool-specific parameter to deal with smaller / larger genomes.

For an in-depth explanation of how to make tool-specific module changes see [here](tool-specific_option_configuration.md)

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
