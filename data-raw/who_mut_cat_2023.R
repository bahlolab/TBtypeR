#!/usr/bin/env Rscript
# Script to derive the TBtypeR WHO mutation catalogue

library(tidyverse)

wd <- file.path(normalizePath('.'), 'data-raw', 'wd')
if (!dir.exists(wd)) { dir.create(wd) }

who_mut_cat_url <- 'https://github.com/GTB-tbsequencing/mutation-catalogue-2023/raw/refs/heads/main/Final%20Result%20Files/WHO-UCN-TB-2023.7-eng.xlsx'
who_mut_cat_fn <- file.path(wd, basename(who_mut_cat_url))
# download mutation catalogue if necessary
if (!file.exists(who_mut_cat_fn)) { download.file(who_mut_cat_url, who_mut_cat_fn) }

who_mut_cat_23 <-
  readxl::read_xlsx(
    who_mut_cat_fn,
    sheet = 1,
    skip = 2,
    col_types = 'text'
  ) %>%
  select(-ncol(.)) %>%
  (function(x) {
    colnames(x) <-
      colnames(x) %>%
      str_remove('\\.\\.\\.[0-9]+') %>%
      make.unique() %>%
      str_replace_all('\\.1$', '_WHO')
    return(x)
  }) %>%
  janitor::clean_names() %>%
  readr::type_convert() %>%
  select(drug, variant, ppv, final_confidence_grading) %>%
  distinct() %>%
  filter(str_detect(final_confidence_grading, 'Assoc w R')) %>%
  na.omit() %>%
  inner_join(
    # coordinates
    readxl::read_xlsx(
      who_mut_cat_fn,
      sheet = 2,
      col_types = 'text'
    ) %>%
      janitor::clean_names() %>%
      readr::type_convert(),
    relationship = 'many-to-many'
  ) %>%
  mutate(chrom = 'AL123456.3') %>%
  select(drug, ppv, chrom, pos = position, ref = reference_nucleotide, alt = alternative_nucleotide) %>%
  inner_join(
    tribble(
      ~abbr, ~drug,
      "LEV", "Levofloxacin",
      "MXF", "Moxifloxacin",
      "RIF", "Rifampicin",
      "STM", "Streptomycin",
      "LZD", "Linezolid",
      "AMI", "Amikacin",
      "CAP", "Capreomycin",
      "KAN", "Kanamycin",
      "ETH", "Ethionamide",
      "INH", "Isoniazid",
      "PZA", "Pyrazinamide",
      "DLM", "Delamanid",
      "EMB", "Ethambutol",
      "BDQ", 'Bedaquiline',
      "CFZ", "Clofazimine"
    )
  ) %>%
  mutate(drug = abbr) %>%
  select(-abbr) %>%
  distinct() %>%
  na.omit() %>%
  group_by(drug, chrom, pos, ref, alt) %>%
  summarise(ppv = max(ppv), .groups = 'drop') %>%
  arrange(chrom, pos, ref, alt, drug) %>%
  mutate(ppv = round(ppv, 4)) %>%
  chop(c(drug, ppv)) %>%
  mutate(
    drug = map_chr(drug, str_c, collapse = '|'),
    ppv =  map_chr(ppv, str_c, collapse = '|')
  )

vcf_tbl <-
  who_mut_cat_23 %>%
  mutate(
    id     = row_number(),
    qual   = '.',
    filter = '.',
    info   = str_c('DR=', drug, ';PPV=', ppv)
  ) %>%
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
who_mut_cat_2023 <-
  read_tsv(norm_vcf_fn, comment = '##', col_types = cols()) %>%
  janitor::clean_names() %>%
  select(chrom=number_chrom, pos, ref, alt, id, info) %>%
  mutate(drugs = str_extract(info, '(?<=DR=)[^;]+'),
         ppv = str_extract(info, '(?<=PPV=).+')) %>%
  select(-info) %>%
  distinct() %>%
  mutate(
    across(c(drugs, ppv),  ~ str_replace_all(., '\\|', ';')),
    reference = 'WHO TEAM Global Tuberculosis Programme, 2023') %>%
  select(chrom, pos, ref, alt, drugs, ppv, reference) %>%
  arrange_all()

# publish files

usethis::use_data(who_mut_cat_2023, overwrite = TRUE, internal = FALSE)

res_dir <- normalizePath(file.path(wd, '..', '..', 'TBtypeNF', 'resources'))

write_tsv(who_mut_cat_2023, file.path(res_dir, 'who_mut_cat_2023.tsv.bz2'))

# cleanup
unlink(wd, recursive = T)
