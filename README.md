
# TBtypeR <img src="misc/TBtypeR_hex_logo.png" align="right" height="138"/>

[![Published in Communications
Biology](https://img.shields.io/badge/Published%20in-Communications%20Biology-green)](https://doi.org/10.1038/s42003-025-07705-9)
[![DOI](https://zenodo.org/badge/794900263.svg)](https://doi.org/10.5281/zenodo.14715907)

`TBtypeR` is an R package for accurate and sensitive quantification of
mixtures of *M. tuberculosis* (MTB) strains from whole genome sequencing
data. `TBtypeR` excels as detecting low-frequency mixed infections that
other tools struggle to detect, maintaining a sensitivity of ~85% for
minor strain frequencies of 2.5%, and ~55% for minor strain frequencies
of 1%. `TBtypeR` is implemented as a standalone R package and as part of
a end-to-end Nextflow pipeline, `TBtypeNF`.

### Performance

Extensive benchmarking of TBtypeR against available tools for MTB
mixture detection is detailed in our [publication in Communications
Biology](https://www.nature.com/articles/s42003-025-07705-9). TBtypeR
has the highest accuracy in prediction of minor strain fractions, with
other tools unable to accurately detect or quantify mixtures below 5%:
<img src="publication/10_data/fig_3_A.png" align="center"/>

### Implementaion

TBtypeR identifies mixtures by modelling reference and alternative
allele counts at polymorphic SNP sites with the binomial distribution.
These sites comprise a phylogenetic SNP “barcode” of over 10,000 sites
which has been derived and harmonised from the following publications:

| Publication | Lineages | Unique SNPs | \# Sublineages |
|----|----|----|----|
| [Napier et al.](https://doi.org/10.1186/s13073-020-00817-3) | L1-L9,La1-La3 | 7572 | 78 |
| [Zwyer et al.](https://doi.org/10.12688/openreseurope.14029.2) | La1-La3 | 1323 | 19 |
| [Thawornwattana et al.](https://doi.org/10.1099/mgen.0.000697) | L2 | 728 | 42 |
| [Coscolla et al.](https://doi.org/10.1099/mgen.0.000477) | L5,L6 | 643 | 16 |
| [Shuaib et al.](https://doi.org/10.3390/genes13060990) | L3 | 637 | 12 |

## The Nextflow Pipeline: TBtypeNF

`TBtypeNF` is a [Nextflow](https://www.nextflow.io/index.html) pipeline
for `TBtypeR` that takes FASTQ files as input and performs FASTQ
preprocessing with [fastp](https://github.com/OpenGene/fastp), read
alignment with [BWA-MEM](https://github.com/lh3/bwa), variant calling
with [BCFtools](https://samtools.github.io/bcftools/bcftools.html), and
quality control report generation using
[SAMtools](https://www.htslib.org/) and
[mosdepth](https://github.com/brentp/mosdepth). Variant calls from
BCFtools are then passed to TBtypeR to identify MTBC lineages and
mixtures. The output is an HTML report with detected MTBC strains and
mixtures frequencies. The Nextflow pipeline is the recommended way to
run TBtypeR because it ensures that variant allele counts are calculated
for all barcode sites, improving accuracy and sensitivity.

`TBtypeNF` requires a sample manifest in TSV format with column names
“sample”, “fastq1” and “fastq2” - see [example
manifest](TBtypeNF/resources/lung_example_manifest.tsv).

### Requirements

- [Nextflow](https://www.nextflow.io/) (≥ 22.03.0)
- Singularity/Apptainer or Docker

### Usage Example

``` bash
# download example manifest
wget https://raw.githubusercontent.com/bahlolab/TBtypeR/main/TBtypeNF/resources/lung_example_manifest.tsv -O my_manifest.tsv
# run the nextflow pipeline
nextflow run bahlolab/TBtypeR/TBtypeNF/main.nf -r main -profile singularity --manifest my_manifest.tsv
```

The above example assumes you are using Singularity/Apptainer as the
container engine, to use docker simply swictch `-profile singularity` to
`-profile docker`

### TBtypeNF Parameters

| Parameter | Description | Default Value |
|----|----|----|
| manifest | Input sample manifest | null |
| id | Run identifier, for naming output files | ‘TBtypeNF-run’ |
| outdir | Output files directory | ‘output’ |
| publish_bams | Save BAM files to output directory | false |
| fast | Run FastTBtypeNF workflow | false |
| max_mix | Maximum number of strains in a mixture to be detected | 3 |
| min_mix_prop | Minimum mixture proportion to be detected | 0.005 |

#### FastTBtypeNF

A faster workflow which skips alignment and quality control reporting is
implemented by using [Fastlin](https://github.com/rderelle/fastlin) to
generate allele counts for use by TBtypeR. This is generally 10x faster
to run and gives near identical results. To use this mode supply the
parameter `--fast` on the command line, e.g.:

``` bash
nextflow run bahlolab/TBtypeR/TBtypeNF/main.nf -r main -profile singularity --manifest my_manifest.tsv --fast
```

## The R package: TBtypeR

The easiest way to use `TBtypeR` is through the `TBtypeNF` pipeline.
However, additional parameters and customisation is available by using
the R package directly. `TBtypeR` can be installed with `devtools` as
follows:

``` r
devtools::install_github("bahlolab/TBtypeR")
```

### Data Requirements

The sensitivity of TBtypeR to detect low-frequency mixtures is dependant
on sequencing coverage. Based on our testing, we recommended to have an
average coverage of at least 20x for detection of 5% mixtures and above,
of at least 40x for mixtures down to 2.5%, and of at least 60x for
mixtures down to 1%.

It is recommended to either use TBtypeNF or [BCFtools
Call](https://samtools.github.io/bcftools/bcftools.html#call) to
generate VCF files for TBtypeR. VCF files generated by other software
may be compatible if the following conditions are met:  
1) **Reference Genome**: TBtypeR expects variants to be called angainst
the [H37Rv reference
genome]('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2/GCA_000195955.2_ASM19595v2_genomic.fna.gz').
The chromosome must be named either “AL123456.3” or “NC_000962.3”.  
2) **AD Format Field**: TBtypeR requires the allelic depth to be stored
in the VCF format field AD, consistent with BCFtools call output.  
3) **Barcode Sites**: Coverage of the majority of TBtypeR SNP barcode
sites for greater sensitivity. Many variant callers will drop non
polymorphic sites within a dataset and will be of limited use for
TBtypeR.

This can be achieved with BCFtools as follows:

``` bash
# download TBtypeR SNP barcode
wget https://github.com/bahlolab/TBtypeR/raw/refs/heads/main/TBtypeNF/resources/tbt_panel.tsv.bz2

# reformat for BCFtools call
bzcat tbt_panel.tsv.bz2 | tail -n+2 | awk 'BEGIN { FS = OFS = "\t" } { print $1, $2, $3 "," $4 }' > tbtyper_targets.tsv

# call variants with BCFtools
bcftools mpileup <SAMPLE_1>.bam -q 1 -Q 10 -d 200 -f <REFERENCE_FASTA> -a FMT/AD -Ou \
  | bcftools call -A -m --prior 1e-2 -C alleles -T tbtyper_targets.tsv -Ou \
  | bcftools annotate -x INFO,^FORMAT/GT,^FORMAT/AD -Oz -o <SAMPLE_1>.vcf.gz
  
# optionally merge samples with BCFtools merge
bcftools merge <SAMPLE_1>.vcf.gz <SAMPLE_2>.vcf.gz <...> --no-index -Oz -o merged.vcf.gz
```

Example usage of `TBtypeR`:

``` r
library(tidyverse)
library(TBtypeR)

# replace with path to your VCF file
vcf_filename <- system.file('vcf/example.vcf.gz', package = 'TBtypeR')

tbtype_result <- 
  # generate TBtypeR results
  tbtype(vcf = vcf_filename) %>% 
  # filter TBtypeR results
  filter_tbtype(max_phylotypes = 3) %>%
  # unnest data so there is 1 row per identified Mtb strain in each sample
  unnest_mixtures()

tbtype_result %>% 
  select(sample_id, n_phy, mix_phylotype, mix_prop) %>% 
  knitr::kable()
```

| sample_id   | n_phy | mix_phylotype | mix_prop |
|:------------|------:|:--------------|---------:|
| SRR13312530 |     2 | 4.2.1         |   0.8579 |
| SRR13312530 |     2 | 4.3.3         |   0.1421 |
| SRR13312531 |     1 | 4.3.3         |   1.0000 |
| SRR13312533 |     2 | 4.3.3         |   0.9192 |
| SRR13312533 |     2 | 4.2.1         |   0.0808 |

Visualise mixtures:

``` r
tbtype_result %>% 
  ggplot(aes(x = sample_id,
             y = mix_prop,
             fill = mix_phylotype)) +
  geom_col() +
  coord_flip() +
  labs(x = 'Sample ID',
       y = 'Minor Strain Fraction (%)', 
       fill = 'Sublineage') +
  theme(text = element_text(size = 6))
```

![](README_files/figure-gfm/visualise_result-1.png)<!-- -->

Detailed usage guides for the `tbtype` and `filter_tbtype` functions are
available in the package documentation by running `help(tbtype)` or
`help(filter_tbtype)`.

## Citation

TBtypeR is published in Communications Biology:

- [Munro, J. E., Coussens, A. K., & Bahlo, M. (2025). TBtypeR: Sensitive
  detection and sublineage classification of Mycobacterium tuberculosis
  complex mixed-strain infections. *Communications Biology*, 8(1),
  260.](https://doi.org/10.1038/s42003-025-07705-9)
