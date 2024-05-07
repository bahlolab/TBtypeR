#!/usr/bin/env Rscript
# Script to derive the TBtypeR WHO mutation catalogue

library(tidyverse)

wd <- file.path(normalizePath('.'), 'data-raw', 'wd')
if (!dir.exists(wd)) { dir.create(wd) }

who_mut_cat_url <- 'https://apps.who.int/iris/bitstream/handle/10665/341906/WHO-UCN-GTB-PCI-2021.7-eng.xlsx'
who_mut_cat_fn <- file.path(wd, basename(who_mut_cat_url))
# download mutation catalogue if necessary
if (!file.exists(who_mut_cat_fn)) { download.file(who_mut_cat_url, who_mut_cat_fn) }

# read Genome_indices sheet
genome_indices <-
  readxl::read_excel(who_mut_cat_fn, sheet = 'Genome_indices', col_types = 'text') %>%
  janitor::clean_names() %>%
  type_convert(col_types = cols()) %>%
  select(variant,
         pos = final_annotation_position,
         ref = final_annotation_reference_nucleotide,
         alt = final_annotation_alternative_nucleotide) %>%
  distinct() %>%
  mutate(ref = str_to_upper(ref),
         alt = str_to_upper(alt),
         pos = as.integer(pos))


# read Mutation_catalogue sheet
mut_cat_raw <-
  readxl::read_excel(who_mut_cat_fn, sheet = 'Mutation_catalogue', skip = 1, col_types = 'text') %>%
  (function(x) {
    at <- which(str_detect(colnames(x), '^\\.\\.\\.'))
    top_names <-
      readxl::read_excel(who_mut_cat_fn, sheet = 'Mutation_catalogue', range = 'A1:AZ1') %>%
      colnames()
    magrittr::set_colnames(x, replace(colnames(x), at, top_names[at]))
  }) %>%
  janitor::clean_names() %>%
  rename(variant = variant_common_name) %>%
  mutate(variant = str_remove_all(variant, '\\s+\\(.+\\)')) %>%
  mutate(across(everything(), ~str_remove(., '%$'))) %>%
  readr::type_convert(col_types = cols()) %>%
  filter(str_detect(final_confidence_grading, 'Assoc w R'))

vcf_tbl <-
  left_join(
    mut_cat_raw,
    genome_indices,
    by = 'variant',

  ) %>%
  select(pos, ref, alt, variant, drug, ppv, final_confidence_grading) %>%
  filter(!is.na(pos)) %>%
  arrange_all() %>%
  distinct() %>%
  mutate(id = row_number()) %>%
  mutate(chrom = 'AL123456.3',
         qual = '.',
         filter = '.',
         info = str_c('DR=', drug, ';PPV=', round(ppv, 4))) %>%
  select(chrom, pos, id, ref, alt, qual, filter, info)



ref_url <- 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2/GCA_000195955.2_ASM19595v2_genomic.fna.gz'
ref_fn <- file.path(wd, basename(ref_url))

# download ref fasta if necessary
if (!file.exists(ref_fn)) { download.file(ref_url, ref_fn) }
ref_bgz <- str_c(str_remove(ref_fn, '\\.gz'), '.bgz')
# bgzip and index ref fasta
if (!file.exists(ref_bgz)) {
  Rsamtools::bgzip(ref_fn)
  Rsamtools::indexFa(ref_bgz)
}

h37rv_ref <-
  Biostrings::readDNAStringSet(ref_fn) %>%
  setNames('AL123456.3')

# NB - software must be available on your system, possibly as a module
bcftools <- system('which bcftools || module load bcftools && which bcftools', intern = TRUE)

vcf_fn <- file.path(wd, 'who_mut_cat.vcf.gz')
norm_vcf_fn <- file.path(wd, 'who_mut_cat.norm.vcf.gz')

# create VCF
vcf_tbl %>%
  unite(line, everything(), sep='\t') %>%
  pull(line) %>%
  (function(x) c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=AL123456.3,length=4411532>",
    "##INFO=<ID=DR,Number=A,Type=String,Description=\"Drug Resistance\">",
    "##INFO=<ID=PPV,Number=A,Type=String,Description=\"Drug Resistance PPV\">",
    "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
    x
  )) %>%
  write_lines(vcf_fn)

# normalise VCF
str_c(
  bcftools,
  'norm -Oz',
  '-f', ref_bgz,
  vcf_fn,
  '-o', norm_vcf_fn,
  sep = ' ') %>%
  system()

# read catalogue from norm VCF
who_mut_cat_2021 <-
  read_tsv(norm_vcf_fn, comment = '##', col_types = cols()) %>%
  janitor::clean_names() %>%
  select(pos, ref, alt, id, info) %>%
  mutate(drug = str_extract(info, '(?<=DR=)[^;]+'),
         ppv = str_extract(info, '(?<=PPV=).+')) %>%
  select(-info, -id) %>%
  distinct() %>%
  arrange_all() %>%
  chop(c(drug, ppv)) %>%
  mutate(drugs = map_chr(drug, ~ str_c(., collapse = ';')),
         ppv = map_chr(ppv, ~ str_c(., collapse = ';'))) %>%
  mutate(chrom = 'AL123456.3',
         reference = 'WHO TEAM Global Tuberculosis Programme, 2021') %>%
  select(chrom, pos, ref, alt, drugs, ppv, reference) %>%
  arrange_all()

# publish files

usethis::use_data(who_mut_cat_2021, overwrite = TRUE, internal = FALSE)

res_dir <- normalizePath(file.path(wd, '..', '..', 'TBtypeNF', 'resources'))

write_tsv(who_mut_cat_2021, file.path(res_dir, 'who_mut_cat_2021.tsv.bz2'))

# cleanup
unlink(wd, recursive = T)
