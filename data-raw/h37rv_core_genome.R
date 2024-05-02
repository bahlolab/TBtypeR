
# List of core genes from https://doi.org/10.3390/antibiotics10050500

wd <- file.path(normalizePath('.'), 'data-raw', 'wd')
if (!dir.exists(wd)) { dir.create(wd) }

core_genes_fn <- file.path(wd, 'core_genes.rds')

if (!file.exists(core_genes_fn)) {
  tab_3_fn <- file.path(wd, 'Table_3.XLSX')
  core_genome_zip_url <- 'https://www.mdpi.com/article/10.3390/antibiotics10050500/s1'
  core_genome_zip_fn <- file.path(wd, 'antibiotics-10-00500-s001.zip')
  download.file(core_genome_zip_url, core_genome_zip_fn)
  unzip(core_genome_zip_fn,
        files = 'antibiotics-1180064-supplementary/Table_3.XLSX',
        junkpaths = TRUE,
        exdir = wd)
  unlink(core_genome_zip_fn)

  core_genes <-
    readxl::read_excel(tab_3_fn) %>%
    pull(1) %>%
    str_extract('^[^,]+')

  file.remove(tab_3_fn)

  saveRDS(core_genes, core_genes_fn)
}
core_genes <- readRDS(core_genes_fn)


gff_url <- 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2/GCA_000195955.2_ASM19595v2_genomic.gff.gz'

h37rv_core_genome <-
  rtracklayer::import.gff(gff_url) %>%
  (function(x) x[x$locus_tag %in% core_genes & x$type == 'CDS'] ) %>%
  (function(x) GenomicRanges::GRanges(GenomicRanges::seqnames(x), x@ranges))

usethis::use_data(h37rv_core_genome, overwrite = TRUE, internal = FALSE)

# create bed file for TBtypeNF
res_dir <- normalizePath(file.path(wd, '..', 'TBtypeNF', 'resources'))
bed_fn <- file.path(res_dir, 'h37rv_core_genome.bed.bz2')
rtracklayer::export.bed(core_reg, bed_fn )
read_tsv(bed_fn, col_names = F) %>% select(1:3) %>% write_tsv(bed_fn, col_names = F)

