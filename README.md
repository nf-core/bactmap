# ![nf-core/bactmap](docs/images/nf-core-bactmap_logo.png)

[![GitHub Actions CI Status](https://github.com/nf-core/bactmap/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/bactmap/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/bactmap/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/bactmap/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/bactmap/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A520.11.0--edge-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23bactmap-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/bactmap)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/bactmap** is a bioinformatics best-practise analysis pipeline for mapping short reads from bacterial WGS to a reference sequence, creating filtered VCF files, making pseudogenomes based on high quality positions in the VCF files and optionally reating a phylogeny from an alignment of the pseudogenomes.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker / Singularity containers making installation trivial and results highly reproducible.

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/bactmap/results).

## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Index the reference sequence using [`bwa index`](https://github.com/lh3/bwa)
2. (Optionally if `params.trim` is set) trim reads using [`fastp`](https://github.com/OpenGene/fastp). The default trimming parameters are

   ```console
   --cut_front --cut_tail --trim_poly_x --cut_mean_quality 30 --qualified_quality_phred 30 --unqualified_percent_limit 10 --length_required 50
   ```

3. (Optionally if `params.depth_cutoff` is set) downsample to the reads based on estimated genome size using [`mash sketch`](https://github.com/marbl/Mash) and the required depth of coverage as set by the `depth_cutoff` params using [`rasusa`](https://github.com/mbhall88/rasusa)
4. Map reads to the indexed reference genome using [`bwa mem`](https://github.com/lh3/bwa) to produce a bam file
5. Sort the bam file using [`samtools`](http://www.htslib.org/doc/samtools.html)
6. Call variants using [`bcftools mpileup`](http://samtools.github.io/bcftools/bcftools.html). A minimum base quality of 20 is used for pre-filtering and the following fields are included in the resulting VCF file `FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR`. A haploid and multiallelic model is assumed. These defaults can be overridden using the `modules.bcftools_mpileup` params. The defaults are:

    ```console
    args  = '--min-BQ 20 --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR'
    args2 = '--ploidy 1 --multiallelic-caller'
    ```

7. Variants are filtered using [`bcftools filter`](http://samtools.github.io/bcftools/bcftools.html). The defaults params in `modules.bcftools_filter` are:

    ```console
    args = '--soft-filter LowQual --exclude "%QUAL<25 || FORMAT/DP<10 || MAX(FORMAT/ADF)<2 || MAX(FORMAT/ADR)<2 || MAX(FORMAT/AD)/SUM(FORMAT/DP)<0.9 || MQ<30 || MQ0F>0.1" --output-type z'
    ```

    There will be a row in filtered VCF file for each position in the reference genome. The position will either have value in the FILTER column of either `PASS` or `LowQual`.
8. Use the filtered VCF to create a pseudogenome based on the reference genome using the [`vcf2pseudogenome.py script`](https://github.com/nf-core/bactmap/blob/dev/bin/vcf2pseudogenome.py). The base in a sample at a position where the VCF file row that has `PASS` in FILTER will be either ref or alt and the appropriate base will be encoded at that position. The base in a sample at a position where the VCF file row that has `LowQual` in FILTER is uncertain will be encoded as a `N` character. Missing data will be encoded as a `-` character.
9. All samples pseudogenomes and the original reference sequence will concatenated together to produce a flush alignment where all the sequence for all samples at all positions in the original reference sequence will be one of `{G,A,T,C,N,-}`. This alignment can be used for other downstream processes such as phylogeny generation.
10. Optionally this alignment can be processed to produce a phylogeny as part of this pipeline. The number of constant sites in the alignment will be determined using [`snp-sites`](https://github.com/sanger-pathogens/snp-sites).
11. (Optionally if `params.remove_recombination` is set) remove regions likely to have been acquired by horizontal transfer and recombination and therefore perturb the true phylogeny using [`gubbins`](https://sanger-pathogens.github.io/gubbins/). This should only be run on sets of samples that are closely related and not for example on a set of samples that have diversity spanning that of the entire species.
12. Depending on the params set, run 0 - 4 tree building algorithms.
    * `params.modules.rapidnj.build = true` Build a neighbour-joining pylogeny using [`rapidnj`](https://birc.au.dk/software/rapidnj)
    * `params.modules.fasttree.build = true` Build an approximately-maximum-likelihood phylogeny using [`FastTree`](http://www.microbesonline.org/fasttree)
    * `params.modules.iqtree.build = true` Build a maximum-likelihood phylogeny using [`IQ-TREE`](http://www.iqtree.org)
    * `params.modules.raxmlng.build = true` Build a maximum-likelihood phylogeny using [`RAxML Next Generation`](https://github.com/amkozlov/raxml-ng)

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Podman`](https://podman.io/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/bactmap -profile test,<docker/singularity/podman/conda/institute>
    ```

    * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    * If you are using `singularity` then the pipeline will auto-detect this and attempt to download the Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. It is also highly recommended to use the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) settings to store the images in a central location for future pipeline runs.
    * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

    <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

    ```bash
    nextflow run nf-core/bactmap -profile <docker/singularity/podman/conda/institute> --input samplesheet.csv --reference chromosome.fasta
    ```

See [usage docs](https://nf-co.re/bactmap/usage) for all of the available options when running the pipeline.

## Documentation

The nf-core/bactmap pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/bactmap/usage) and [output](https://nf-co.re/bactmap/output).

## Credits

nf-core/bactmap was originally written by [Anthony Underwood](https://github.com/aunderwo), [Adries Van Tonder](https://github.com/avantonder) and [Thanh Le Viet](https://github.com/thanhleviet).

We thank the following people for their extensive assistance in the development
of this pipeline:

* [Harshil Patel](https://github.com/drpatelh)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#bactmap` channel](https://nfcore.slack.com/channels/bactmap) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below, update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/bactmap for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->
An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
