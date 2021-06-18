# nf-core/bactmap: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.0.0 - Aluminium Spider - 2021-06-18

Initial release of nf-core/bactmap, created with the [nf-core](https://nf-co.re/) template.

The pipeline is composed of the following steps:

1. Index reference fasta file ([`BWA index`](https://github.com/lh3/bwa))
2. Trim reads for quality and adapter sequence (Optional) ([`fastp`](https://github.com/OpenGene/fastp))
3. Estimate genome size ([`mash sketch`](https://mash.readthedocs.io/en/latest/index.html))
4. Downsample fastq files (Optional) ([`Rasusa`](https://github.com/mbhall88/rasusa))
5. Variant calling
    1. Read mapping ([`BWA mem`](https://github.com/lh3/bwa))
    2. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    3. Call and filter variants ([`BCFtools`](http://samtools.github.io/bcftools/bcftools.html))
    4. Convert filtered bcf to pseudogenome fasta ([`vcf2pseudogenome.py`](https://github.com/nf-core/bactmap/blob/dev/bin/vcf2pseudogenome.py))
6. Create alignment from pseudogenome by concatenating fasta files having first checked that the sample sequences are high quality([`calculate_fraction_of_non_GATC_bases.py`](https://github.com/nf-core/bactmap/blob/dev/bin/))
7. Remove recombination (Optional) ([`Gubbins`](https://sanger-pathogens.github.io/gubbins/))
8. Extract variant sites from alignment ([`SNP-sites`](https://github.com/sanger-pathogens/snp-sites))
9. Construct phylogenetic tree (Optional)
    1. Fast/less accurate
        * neighbour joining [`RapidNJ`](https://birc.au.dk/software/rapidnj/)
        * approximate maximum likelihood [`FastTree2`](http://www.microbesonline.org/fasttree/))
    2. Slow/more accurate, maximum likelihood
        * [`IQ-TREE`](http://www.iqtree.org/),
        * [`RAxML-NG`](https://github.com/amkozlov/raxml-ng)
