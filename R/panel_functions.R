
panel_to_vid <- function(panel) {
  tidyr::unite(panel, "vid", chrom, pos, ref, alt)
}

#' @importFrom dplyr bind_cols
panel_with_vid <- function(panel) {
  panel_to_vid(panel) %>%
    select(vid) %>%
    bind_cols(panel)
}

#' @importFrom tidyr pivot_wider unnest
panel_to_geno <- function(panel) {
  phylo <- panel_to_phylo(panel)

  panel_vid <-
    panel_to_vid(panel) %>%
    select(phylotype, vid, genotype)

  geno <-
    tidytree::as_tibble(phylo) %>%
    transmute(
      phylotype = label,
      parent = get_anc_list(phylo, node, as_label = TRUE) %>% map2(., label, c)
    ) %>%
    (function(x) {
      inner_join(
        unnest(x, parent),
        select(panel_vid, parent = phylotype, vid),
        by = "parent"
      ) %>%
        select(-parent)
    }) %>%
    mutate(gt = 1L) %>%
    pivot_wider(
      names_from = vid,
      values_from = gt,
      values_fill = 0L
    ) %>%
    arrange(phylotype) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("phylotype") %>%
    as.matrix() %>%
    (function(x) {
      x <- x[, panel_vid$vid]
      x[, panel_vid$genotype == 0L] <-
        abs(1L - x[, panel_vid$genotype == 0L])
      x
    })

  geno <- rbind(
    geno,
    tidytree::as_tibble(phylo) %>%
      as_tibble() %>%
      filter(parent == treeio::rootnode(phylo), node != parent) %>%
      with(colMeans(geno[label, , drop = F])) %>%
      (function(x) if_else(x > 0.5, 1L, 0L)) %>%
      matrix(nrow = 1, dimnames = list("root", NULL))
  )

  geno <- geno[phylo_labels(phylo), ]

  return(geno)
}

panel_to_phylo <- function(panel) {
  check_panel(panel)

  edges <-
    panel %>%
    select(from = parent_phylotype, to = phylotype) %>%
    distinct()

  # check acyclic
  gr <- igraph::graph_from_edgelist(as.matrix(edges), directed = T)
  assertthat::assert_that(igraph::is.dag(gr))

  root <- setdiff(edges$from, edges$to)
  tips <- setdiff(edges$to, edges$from)
  nodes <- c(root, setdiff(c(edges$from, edges$to), c(root, tips)))

  assertthat::assert_that(
    max(table(edges$to)) == 1,
    length(root) == 1
  )

  numbered <-
    c(
      setNames(seq_along(tips), tips),
      setNames(seq_along(nodes), nodes) + length(tips)
    )

  phylo <- list(
    edge = unname(cbind(numbered[edges$from], numbered[edges$to])),
    tip.label = tips,
    node.label = nodes,
    Nnode = length(nodes)
  )

  class(phylo) <- "phylo"
  phylo <- ape::ladderize(phylo, right = FALSE)
  check_phylo(phylo)

  return(phylo)
}


check_panel <- function(panel, mode = c("phylo", "none", "dr")) {
  mode <- match.arg(mode)

  assertthat::assert_that(
    is.data.frame(panel),
    all(c("chrom", "pos", "ref", "alt") %in% colnames(panel)),
    rlang::is_integerish(panel$pos),
    all(nchar(panel$ref) > 0),
    all(nchar(panel$alt) > 0),
    all(stringr::str_detect(panel$ref, "^[ACTG]+$")),
    all(stringr::str_detect(panel$alt, "^[ACTG]+$")),
    mode != "phylo" || all(c("genotype", "phylotype", "parent_phylotype") %in% colnames(panel)),
    mode != "phylo" || all(panel$genotype %in% c(0L, 1L)),
    mode != "dr" || "drugs" %in% colnames(panel)
  )

  invisible(TRUE)
}

check_phylo <- function(phy) {
  # Based on ape::checkValidPhylo

  assertthat::assert_that(
    is.character(phy$tip.label),
    length(phy$tip.label) > 0,
    is_scalar_integer(phy$Nnode),
    phy$Nnode > 0,
    is.matrix(phy$edge),
    is_integerish(phy$edge),
    all(phy$edge > 0),
    ncol(phy$edge) == 2
  )

  n <- length(phy$tip.label)
  m <- phy$Nnode
  tab <- tabulate(phy$edge)

  assertthat::assert_that(
    # some numbers in 'edge' are larger than 'n + m'
    length(tab) <= n + m,
    # each tip must appear once in 'edge'
    all(tab[seq_len(n)] == 1),
    # all nodes should appear at least twice in 'edge'
    all(tab[n + seq_len(m)] > 1),
    # nodes and tips should appear only once in the 2nd column of 'edge'
    all(table(phy$edge[, 2]) == 1),
    # the root node should not appear in the 2nd column of 'edge'
    !any(phy$edge[, 2] == n + 1L),
    # tips should not appear in the 1st column of 'edge'
    !any(phy$edge[, 1] <= n & phy$edge[, 1] > 0)
  )
}


#' @export
# check input files are compatible and correct
check_inputs <- function(ref_fasta,
                         panel,
                         dr_panel = NULL,
                         calling_regions = NULL
                         ) {

  # ref_fasta = '/vast/scratch/users/munro.j/runs/tbt_test/work/stage/9c/b9f8c9439d2a6d03335a09741dbadd/GCA_000195955.2_ASM19595v2_genomic.fna.gz'
  # panel = '~/pipelines/nf-tbtyper/resources/tbt_panel.tsv.gz'
  # calling_regions = '~/pipelines/nf-tbtyper/resources/ASM19595v2_core_genome.bed.gz'
  # dr_panel = '~/pipelines/nf-tbtyper/resources/who_dr_panel.tsv.gz'

  assert_that(
    is_string(ref_fasta) && file.exists(ref_fasta),
    is.data.frame(panel) ||
      (is_string(panel) && file.exists(panel)),
    is.null(calling_regions) ||
      inherits(calling_regions, 'GenomicRanges') ||
      (is_string(calling_regions) && file.exists(calling_regions)),
    is.null(dr_panel) ||
      is.data.frame(dr_panel) ||
      (is_string(dr_panel) && file.exists(dr_panel))
  )

  ref_bs <- Biostrings::readDNAStringSet(ref_fasta)
  ref_gr <- GenomicRanges::GRanges(str_extract(names(ref_bs), '^[^\\s]+'),
                                   IRanges::IRanges(start = 1, end =  Biostrings::width(ref_bs)))

  # check panel
  if (!is.data.frame(panel)) {
    panel <- read_panel(panel)
  }
  panel_gr <- with(panel, GenomicRanges::GRanges(chrom, IRanges::IRanges(start = pos, width = nchar(ref))))
  assert_that(all(GenomicRanges::countOverlaps(panel_gr,  ref_gr) > 0))

  # check dr_panel
  if (!is.null(dr_panel)) {
    if (!is.data.frame(dr_panel)) {
      dr_panel <- read_panel(dr_panel, phylo = FALSE)
    }
    dr_panel_gr <- with(dr_panel, GenomicRanges::GRanges(chrom, IRanges::IRanges(start = pos, width = nchar(ref))))
    assert_that(all(GenomicRanges::countOverlaps(dr_panel_gr,  ref_gr) > 0))
  }

  # check calling_regions
  if (!is.null(calling_regions)) {
    if (!inherits(calling_regions, 'GenomicRanges')) {
      calling_regions <- rtracklayer::import.bed(calling_regions)
    }
    assert_that(all(GenomicRanges::countOverlaps(calling_regions,  ref_gr) > 0))
  }

  return(TRUE)
}


#' @export
#' @importFrom dplyr bind_rows select
#' @importFrom tidyr chop
#' @importFrom purrr map2_chr
#' @importFrom stringr str_c str_split str_extract
create_calling_targets <- function(ref_fasta,
                                   panel,
                                   dr_panel = NULL,
                                   calling_regions = NULL,
                                   output_dir = getwd(),
                                   indel_pad = 100) {


  assert_that(
    is_scalar_integerish(indel_pad),
    indel_pad >= 0,
    is_string(output_dir),
    check_inputs(ref_fasta = ref_fasta,
                 panel = panel,
                 dr_panel = dr_panel,
                 calling_regions = calling_regions)
  )


  if (!dir.exists(output_dir)) {
    assert_that(dir.create(output_dir, recursive = T))
  }

  snp_targets_fn <- file.path(output_dir, "snp_targets.tsv.gz")
  indel_regions_fn <- file.path(output_dir, "indel_regions.bed.gz")

  snps <- c("A", "C", "G", "T")
  vars_set <-
    map_chr(snps, ~ str_c(union(., snps), collapse = ",")) %>%
    setNames(snps)


  if (is.data.frame(panel)) {
    comb_panel <- panel
  } else {
    comb_panel <- read_panel(panel)
  }
  if (!is.null(dr_panel)) {
    if (!is.data.frame(dr_panel)) {
      dr_panel <- read_panel(dr_panel, phylo = FALSE)
    }
    comb_panel <- bind_rows(comb_panel, dr_panel)
  }

  # create snp_targets
  snp_targets <-
    comb_panel %>%
    filter(nchar(ref) == 1, nchar(alt) == 1) %>%
    select(chrom, pos, ref) %>%
    distinct() %>%
    mutate(vars = vars_set[ref]) %>%
    select(chrom, pos, vars) %>%
    arrange_all()

  if (!is.null(calling_regions)) {

    if (!inherits(calling_regions, 'GenomicRanges')) {
      calling_regions <- rtracklayer::import.bed(calling_regions)
    }
    ref <- Biostrings::readDNAStringSet(ref_fasta)
    names(ref) <- str_extract(names(ref), '^[^\\s]+')
    seqs <- BSgenome::getSeq(ref, calling_regions)

    snp_targets <-
      tibble(chrom = as.character(GenomicRanges::seqnames(calling_regions)),
             pos = GenomicRanges::start(calling_regions),
             ref = as.character(seqs)) %>%
      mutate(ref = map(ref, ~ c(str_split(., '', simplify = T))),
             offset = map(ref, seq_along)) %>%
      unnest(c(ref, offset)) %>%
      mutate(pos = pos + offset - 1L) %>%
      mutate(vars = vars_set[ref]) %>%
      select(chrom, pos, vars) %>%
      bind_rows(snp_targets) %>%
      distinct() %>%
      arrange_all()
  }


  if (nrow(snp_targets)) {
    readr::write_tsv(snp_targets, snp_targets_fn, col_names = FALSE, progress = F)
    message("Created SNP targets file: ", snp_targets_fn)
  } else {
    message("No SNP targets")
  }

  # create indel regions
  indel_regions <-
    comb_panel %>%
    filter(!(nchar(ref) == 1 & nchar(alt) == 1)) %>%
    with(GenomicRanges::GRanges(chrom,
                                IRanges::IRanges(start = pos, width = nchar(ref)))) %>%
    (function(x) GenomicRanges::reduce(x + indel_pad)) %>%
    sort()

  if (length(indel_regions)) {
    rtracklayer::export.bed(indel_regions, indel_regions_fn)
    message("Created indel regions file: ", indel_regions_fn)
  } else {
    message("No indel regions")
  }
}

#' @export
read_panel <- function(filename,
                       phylo = TRUE) {

  if (phylo) {
    panel <-
      readr::read_tsv(filename, col_types = readr::cols(
        chrom = readr::col_character(),
        pos = readr::col_integer(),
        ref = readr::col_character(),
        alt = readr::col_character(),
        genotype = readr::col_integer(),
        phylotype = readr::col_character(),
        parent_phylotype = readr::col_character()
      ))

    check_panel(panel)

  } else {
    panel <-
      readr::read_tsv(filename, col_types = readr::cols(
        chrom = readr::col_character(),
        pos = readr::col_integer(),
        ref = readr::col_character(),
        alt = readr::col_character(),
        drugs = readr::col_character()
      ))

    check_panel(panel, 'dr')

  }

  return(panel)
}

