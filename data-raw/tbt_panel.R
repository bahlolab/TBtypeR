library(tidyverse)
jemr::fix_tidyverse_conflicts()


# set and create working directory
wd <- file.path(normalizePath('.'), 'data-raw')
if (!dir.exists(wd)) { dir.create(wd) }

# download and read h37rv_ref
ref_url <- 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2/GCA_000195955.2_ASM19595v2_genomic.fna.gz'
ref_fn <- file.path(wd, basename(ref_url))
if (!file.exists(ref_fn)) { download.file(ref_url, ref_fn) }

h37rv_ref <-
  Biostrings::readDNAStringSet(ref_fn) %>%
  setNames('Chromosome')


# function to filter SNP barcodes
filter_bc <- function(bc_df, ref_bs) {
  flt_bc_df <-
    bc_df %>%
    filter(ref %in% c('A', 'C', 'G', 'T'),
           alt %in% c('A', 'C', 'G', 'T')) %>%
    filter(ref == as.character(
      BSgenome::getSeq(ref_bs,
                       GenomicRanges::GRanges('Chromosome', IRanges::IRanges(pos, width = 1))))) %>%
    na.omit() %>%
    add_count(pos) %>% filter(n == 1) %>% select(-n)

  if (nrow(flt_bc_df) < nrow(bc_df)) {
    message('removed ',  nrow(bc_df) - nrow(flt_bc_df), ' of ',  nrow(bc_df), ' row(s).')
  }
  return(flt_bc_df)
}

# ---------- NAPIER ----------
# SNP barcode covering wide range of TB lineages, best for lineages 4 and 1
# https://doi.org/10.1186/s13073-020-00817-3
napier_url <- 'https://raw.githubusercontent.com/GaryNapier/tb-lineages/main/fst_results_clean_fst_1_for_paper.csv'
napier_fn <- file.path(wd, 'napier_bc.rds')

if (!file.exists(napier_fn)) {
  read_csv(napier_url, col_types = cols(.default = "c")) %>%
    readr::type_convert(col_types = cols()) %>%
    janitor::clean_names() %>%
    select(lineage=new_lin, pos, ref, alt) %>%
    mutate(reference = 'Napier et al. 2020',
           date = lubridate::dmy('14 December 2020')) %>%
    filter_bc(h37rv_ref) %>%
    saveRDS(napier_fn)
}
napier_bc <- read_rds(napier_fn)

# ---------- COSCOLLA ----------
# detailed SNP barcode for lineage 5, 6 and 9
# https://doi.org/10.1099/mgen.0.000477
# requires manual download of Supplementary material 1 to {wd}/000477_1.xlsx

coscolla_bc_l5_l6_fn <- file.path(wd, 'coscolla_bc_l5_l6.rds')
if(!file.exists(coscolla_bc_l5_l6_fn)) {
  readxl::read_xlsx(file.path(wd, '000477_1.xlsx'), sheet ='Table S7', skip = 2) %>%
    janitor::clean_names() %>%
    select_if(~!all(is.na(.))) %>%
    mutate(sublineage =
             str_remove(sublineage, '^L') %>%
             map(str_split, '', simplify = T) %>%
             map_chr(str_c, collapse = '.')) %>%
    select(lineage = sublineage, pos = genomic_position, ref = ancestral_base, alt = mutation) %>%
    filter_bc(h37rv_ref) %>%
    saveRDS(coscolla_bc_l5_l6_fn)
}
coscolla_bc_l5_l6 <- readRDS(coscolla_bc_l5_l6_fn)

coscolla_bc_l9_fn <- file.path(wd, 'coscolla_bc_l9.rds')
if (!file.exists(coscolla_bc_l9_fn)){
  readxl::read_xlsx(file.path(wd, '000477_1.xlsx'), sheet ='Table S2', skip = 3) %>%
    janitor::clean_names() %>%
    select_if(~!all(is.na(.))) %>%
    mutate(lineage = '9') %>%
    select(lineage, pos = genome_position, ref = reference_base, alt = alternative_base) %>%
    filter_bc(h37rv_ref) %>%
    saveRDS(coscolla_bc_l9_fn)
}
coscolla_bc_l9 <- readRDS(coscolla_bc_l9_fn)

coscolla_bc <-
  bind_rows(coscolla_bc_l5_l6,
            coscolla_bc_l9) %>%
  mutate(reference = 'Coscolla et al. 2021',
         date = lubridate::dmy('08 February 2021'))


# ---------- THAWORNWATTANA ----------
# Lineage 2
# https://doi.org/10.1099/mgen.0.000697
thaworn_url <- 'https://figshare.com/ndownloader/files/28249617'
thaworn_fn <- file.path(wd, 'thaworn_bc.rds')

if (!file.exists(thaworn_fn)) {
  read_csv(thaworn_url, col_types = cols()) %>%
    mutate(lineage = str_remove(lineage, 'L')) %>%
    separate('allele_change', c('ref', 'alt'), sep = '/') %>%
    select(lineage, pos = position, ref, alt) %>%
    filter(lineage != '2.2.1.Modern') %>%
    filter_bc(h37rv_ref) %>%
    mutate(reference = 'Thawornwattana et al. 2021',
           date = lubridate::dmy('17 November 2021')) %>%
    saveRDS(thaworn_fn)
}

thaworn_bc <- read_rds(thaworn_fn)


# ---------- SHUAIB ----------
# lineage 3
# https://doi.org/10.3390/genes13060990
# requires download of Supplementary material S4 to {wd}/genes13060990_S4.xlsx
shuaib_fn <- file.path(wd, 'shuaib_bc.rds')
if (!file.exists(shuaib_fn)) {
  readxl::read_excel(file.path(wd, 'genes13060990_S4.xlsx')) %>%
    janitor::clean_names() %>%
    select(lineage = number_group_id,
           pos,
           ref = reference_allel,
           alt = group_allel) %>%
    mutate(lineage = str_remove(lineage, 'group_')) %>%
    filter_bc(h37rv_ref) %>%
    mutate(reference = 'Shuaib et al. 2022',
           date = lubridate::dmy('31 May 2022')) %>%
    saveRDS(shuaib_fn)
}
shuaib_bc <- read_rds(shuaib_fn)

# ---------- COMBINE BARCODES -----------

# check overlaps/ disagreements
bind_rows(
  napier_bc,
  coscolla_bc,
  thaworn_bc,
  shuaib_bc) %>%
  mutate(reference = word(reference) %>% str_to_lower()) %>%
  add_count(pos, ref, alt) %>% filter(n > 1) %>%
  # select(pos, lineage, reference) %>%
  group_by(pos) %>%
  summarise(refs = list(sort(unique(reference))),
            lins = list(lineage[match(refs[[1]], reference)])) %>%
  mutate(across(c(refs, lins), map_chr, \(x) str_c(x, collapse = ' -- '))) %>%
  count(refs, lins) %>%
  arrange(desc(n))


combined_bc <-
  bind_rows(
    napier_bc %>%
      filter(! str_starts(lineage, '2\\.'), # prefer thaworn_bc
             ! str_starts(lineage, '3\\.')), # prefer shuaib_bc
    coscolla_bc,
    thaworn_bc,
    shuaib_bc) %>%
  arrange(pos, ref, alt) %>%
  group_by(pos) %>%
  slice(which.min(date)) %>%
  ungroup() %>%
  mutate(genotype = if_else(lineage %in% c('4', '4.9'), 0L, 1L))

# clean up and add parents
tbt_panel <-
  combined_bc %>%
  select(lineage) %>%
  distinct() %>%
  mutate(parent = case_when(
    str_detect(lineage, '^[0-9]$') ~ 'root',
    str_detect(lineage, '^M\\.')   ~ 'root',
    lineage %in% c('2.2.A', '2.2.AA2', '2.2.AA3', '2.2.AA4', '2.2.B', '2.2.C', '2.2.D',
                   '2.2.E', '2.2.M1', '2.2.M2', '2.2.M3', '2.2.M4', '2.2.M5', '2.2.M6')
    ~ '2.2.1',
    {
      p <- str_extract(lineage, '^.+(?=\\.[^.]+$)')
      p %in% lineage               ~ p
    },
    {
      p2 <- str_extract(lineage, '^.+(?=(\\.[^.]+){2}$)')
      p2 %in% lineage               ~ p2
    },
  ))  %>%
  left_join(combined_bc, by = 'lineage') %>%
  mutate(chrom = 'AL123456.3') %>%
  select(chrom, pos, ref, alt, genotype, phylotype = lineage, parent_phylotype = parent, reference) %>%
  arrange_all()

usethis::use_data(tbt_panel, overwrite = TRUE, internal = FALSE)
