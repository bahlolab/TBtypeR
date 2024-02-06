
#' @importFrom readr write_lines
#' @importFrom readr write_tsv
#' @importFrom stringr str_remove str_remove_all str_sub str_detect
make_fastlin_bc <- function(
    output = 'TBtypeR.fastlin_bc.tsv',
    mode = c('TBtypeR', 'fastlin'),
    ref_fasta = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2/GCA_000195955.2_ASM19595v2_genomic.fna.gz',
    panel = TBtypeR::tbt_panel,
    genome_size = 4411532,
    window = 50L
)
{

  mode <- match.arg(mode)

  if (!is.data.frame(panel)) {
    panel <- read_panel(panel)
  }

  ref_genome <-
    Biostrings::readDNAStringSet(ref_fasta) %>%
    setNames('Chromosome')

  write_lines(
    str_c('genome_size\t', genome_size),
    output)

  panel %>%
    (\(x) `if`(mode == 'fastlin', rename_panel_fastlin(x), x)) %>%
    filter(pos > window, pos < genome_size - window) %>% # must allow +/-50 bp
    mutate(
      left =
        BSgenome::getSeq(
          ref_genome,
          GenomicRanges::GRanges(
            'Chromosome',
            IRanges::IRanges(pos-window, width = window))
        ) %>%
        as.character(),
      right =
        BSgenome::getSeq(
          ref_genome,
          GenomicRanges::GRanges(
            'Chromosome',
            IRanges::IRanges(pos+1, width = window))
        ) %>%
        as.character()) %>%
    (function(x) {
      if (mode == 'TBtypeR') {
        select(x, pos, ref, alt, left, right) %>%
          (\(x) bind_rows(
            mutate(x, site = str_c('p', pos, '.', ref)) %>%
              select(site, left, snp = ref, right),
            mutate(x, site = str_c('p', pos, '.', alt)) %>%
              select(site, left, snp = alt, right)
          ))
      } else {
        bind_rows(
          filter(x, genotype == 1) %>%
            select(phylotype, left, snp = alt, right) %>%
            # remove redundant 25-mers
            mutate(
              kmer = str_c(str_sub(left, start = window - 12 + 1), snp, str_sub(right, end = 12)),
              n = Biostrings::countPDict(Biostrings::PDict(kmer), ref_genome[[1]]),
              n = n + Biostrings::countPDict(Biostrings::PDict(kmer), Biostrings::reverseComplement(ref_genome[[1]]))
            ) %>%
            filter(n == 0) %>%
            select(-kmer, -n)
          ,
          filter(x, genotype == 0) %>%
            select(phylotype, left, snp = ref, right),
        )
      }
    }) %>%
    arrange_all() %>%
    write_tsv(output, append = T)
}

fastlin_allele_counts <- function(
    fl_output_tsv,
    panel = TBtypeR::tbt_panel,
    sample_id = 'sample'
)
{
  fl_cov <-
    read_tsv(fl_output_tsv, col_types = cols()) %>%
    select(bc = log_barcodes) %>%
    separate_rows(bc, sep = ',') %>%
    transmute(
      pos =
        str_extract(bc, '(?<=p)[0-9]+') %>%
        as.integer(),
      snp = str_extract(bc, '(?<=\\.)[ACGT]'),
      cov =
        str_extract(bc, '(?<=\\()[0-9]+') %>%
        as.integer())

  ac_tbl <-
    panel_with_vid(tbt_panel) %>%
    select(vid, pos, ref, alt) %>%
    left_join(select(fl_cov, pos, ref = snp, ac_ref = cov),
              by = c('pos', 'ref')) %>%
    left_join(select(fl_cov, pos, alt = snp, ac_alt = cov),
              by = c('pos', 'alt')) %>%
    filter(!(is.na(ac_ref) & is.na(ac_alt))) %>%
    mutate(across(starts_with('ac_'), \(x) replace_na(x, 0L)))

  ac_array <-
    abind::abind(
      matrix(ac_tbl$ac_ref, nrow = 1),
      matrix(ac_tbl$ac_alt, nrow = 1),
      along = 3)

  dimnames(ac_array) <- list(sample = sample_id, variant = ac_tbl$vid, allele = c('ref', 'alt'))
  return(ac_array)
}

rename_panel_fastlin <- function(panel) {
  # fastlin expect parent to be strict prefix
  orig <-
    panel %>%
    select(phylo = phylotype, parent = parent_phylotype) %>%
    distinct()

  new <- with(orig, str_remove(phylo, str_c('^', parent, '\\.')))
  new <- str_c(new, '[', seq_along(new), ']')
  names(new) <- orig$phylo
  mod <-
    orig %>%
    mutate(phylo = unname(new[phylo]),
           parent = if_else(!is.na(new[parent]),
                            unname(new[parent]),
                            parent))

  parents <- filter(mod, parent == 'root') %>% pull(phylo)
  while (TRUE) {
    childs <- filter(mod, parent %in% parents) %>% pull(phylo)
    if (length(childs) == 0) {
      break
    }
    new <-
      mod %>%
      with(if_else(
        phylo %in% childs,
        str_c(parent, '.', phylo),
        phylo
      ))
    names(new) <- mod$phylo

    mod <-
      mod %>%
      mutate(phylo = new[phylo],
             parent = if_else(!is.na(new[parent]), new[parent], parent))

    parents <- unname(new[childs])
  }

  mod <-
    mod %>%
    mutate(across(everything(), \(x) str_remove_all(x, '\\[[0-9]+\\]')))

  old <- with(orig, unique(c(phylo, parent)))
  new <- with(mod, unique(c(phylo, parent)))
  names(new) <- old

  panel %>%
    mutate(phylotype = new[phylotype],
           parent_phylotype = new[parent_phylotype])

}

