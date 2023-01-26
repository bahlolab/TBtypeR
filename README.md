TBtypeR
================

**TBtypeR** is an R package that performs accurate and sensitive
quantification of mixtures of *M. tuberculosis* (MTB) strains from whole
genome sequencing (WGS) data. In addition TBtypeR performs drug
resistance prediction based on the [2021 WHO drug resistance
catalogue](https://doi.org/10.1016/S2666-5247(21)00301-3). `NF-TBtypeR`
is an end-to-end pipeline for processing MTB WGS data

## Runinning the Pipeline

## Installation

`TBtypeR` can be installed with `devtools` as follows:

``` r
devtools::install_github("bahlolab/TBtypeR")
```

## Preparing input data

It is recommended to process data with the included Nextflow pipeline,
NF-TBtypeR. The pipeline will take FASTQ files as input, take them
though alignment and variant calling with bwa mem, samtools and
BCFtools, and output VCF & GDS format variants as well as optionally
TBtypeR results.

## Running TBtypeR

## Visualising Results

# TBtypeR

### **Notice**

`TBtypeR` is still in the development stages and and should not be used
in production. Documentation is currently poor and functionality is
subject to change.

### Basic usage

1)  Convert VCF file to GDS file using `SeqArray`
2)  Extract allele counts from GDS
3)  Fit phylotypes using TBtypeR
4)  Filter results

``` r
library(tidyverse)
library(SeqArray)
library(TBtypeR)
library(future)

# set desired parallelism with future::plan
plan(multiprocess, workers = 4)

# 1) convert VCF file to GDS file
# Note: Variants should be called against the H37RV genome
vcf_fn <- '/path/to/my.vcf.gz'
gds_fn <- '/path/to/my.gds'
seqVCF2GDS(vcf_fn, gds_fn, storage.option = 'ZIP_RA', ignore.chr.prefix = NA_character_)
gds <- seqOpen(gds_fn, allow.duplicate = TRUE)

# 2) Extract allele counts from GDS
allele_counts <- get_allele_counts_gds(gds)

# 3) Fit phylotypes with TBtypeR
results <- fit_phylotypes(allele_counts)

# 4) Filter for best match per sample
filtered_results <- 
  results %>% 
  filter(p_val < 0.001, abs_diff > 10) %>%
  group_by(sample_id) %>%
  arrange(desc(likelihood)) %>%
  slice(1) %>%
  ungroup()
```
