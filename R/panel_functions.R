
panel_to_vid <- function(panel) {
    tidyr::unite(panel, 'vid', chrom, pos, ref, alt)
}

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
    transmute(phylotype = label,
              parent = get_anc_list(phylo, node, as_label = TRUE) %>% map2(., label, c)) %>%
    (function(x) {
      inner_join(
        unnest(x, parent),
        select(panel_vid, parent = phylotype, vid),
        by = 'parent') %>%
        select(-parent)
    }) %>%
    mutate(gt = 1L) %>%
    pivot_wider(names_from = vid,
                values_from = gt,
                values_fill = 0L) %>%
    arrange(phylotype) %>%
    as.data.frame() %>%
    tibble::column_to_rownames('phylotype') %>%
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
      filter(parent == rootnode(phylo), node != parent) %>%
      with(colMeans(geno[label, , drop =F])) %>%
      (function(x) if_else(x > 0.5, 1L, 0L)) %>%
      matrix(nrow = 1, dimnames = list('root', NULL)))

  geno <- geno[phylo_labels(phylo), ]

  return(geno)
}

panel_to_phylo <- function(panel) {

  check_panel(panel)

  edges <-
    panel %>%
    select(from = parent_phylotype, to = phylotype) %>%
    distinct()

  root <- setdiff(edges$from, edges$to)
  tips <- setdiff(edges$to, edges$from)
  nodes <- c(root, setdiff(c(edges$from, edges$to), c(root, tips)))

  assertthat::assert_that(
    max(table(edges$to)) == 1,
    length(root) == 1)

  numbered <-
    c(setNames(seq_along(tips), tips),
      setNames(seq_along(nodes), nodes) + length(tips))

  phylo <- list(
    edge = unname(cbind(numbered[edges$from], numbered[edges$to])),
    tip.label = tips,
    node.label = nodes,
    Nnode = length(nodes))

  class(phylo) <- "phylo"
  phylo <- ape::ladderize(phylo, right=FALSE)
  check_phylo(phylo)

  return(phylo)
}


check_panel <- function(panel, phylotype = TRUE) {


  assertthat::assert_that(
    is.data.frame(panel),
    all(c('chrom', 'pos', 'ref', 'alt') %in% colnames(panel)),
    rlang::is_integerish(panel$pos),
    all(nchar(panel$ref) > 0),
    all(nchar(panel$alt) > 0),
    all(stringr::str_detect(panel$ref, '^[ACTG]+$')),
    all(stringr::str_detect(panel$alt, '^[ACTG]+$')),
    !phylotype || all(c('genotype', 'phylotype', 'parent_phylotype') %in% colnames(panel)),
    !phylotype || all(panel$genotype %in% c(0L, 1L)))


  invisible(TRUE)
}

check_phylo <-  function(phy) {
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
#' @importFrom dplyr bind_rows select
#' @importFrom tidyr chop
#' @importFrom purrr map2_chr
#' @importFrom stringr str_c
create_bcftools_targets <- function(dir,
                                    panel = bind_rows(
                                      TBtyper::tbt_panel,
                                      TBtyper::who_dr_panel),
                                    indel_pad = 100) {

  assert_that(check_panel(panel, phylotype = FALSE),
              is_scalar_integerish(indel_pad),
              indel_pad >= 0)

  if (!dir.exists(dir)) { assert_that(dir.create(dir, recursive = T)) }

  snp_targets_fn <- file.path(dir, 'snp_targets.tsv.gz')
  indel_regions_fn <- file.path(dir, 'indel_regions.bed.gz')

  snps <- c('A', 'C', 'G', 'T')
  # create snp_targets
  snp_targets <-
    panel %>%
    filter(nchar(ref) == 1, nchar(alt) == 1) %>%
    select(chrom, pos, ref) %>%
    distinct() %>%
    mutate(vars = map_chr(ref, ~ str_c(union(., snps), collapse = ','))) %>%
    select(chrom, pos, vars) %>%
    arrange_all()

  if (nrow(snp_targets)) {
    readr::write_tsv(snp_targets, snp_targets_fn, col_names = FALSE, progress = F)
    message('Created SNP targets file: ', snp_targets_fn)
  } else {
    message('No SNP targets')
  }


  # create indel regions
  indel_regions <-
    panel %>%
    filter(!(nchar(ref) == 1 & nchar(alt) == 1)) %>%
    with(GenomicRanges::GRanges(chrom, IRanges::IRanges(start = pos, width = nchar(ref)))) %>%
    (function(x) GenomicRanges::reduce(x + indel_pad)) %>%
    sort()

  if (length(indel_regions)) {
    rtracklayer::export.bed(indel_regions, indel_regions_fn)
    message('Created indel regions file: ', indel_regions_fn)
  } else {
    message('No indel regions')
  }

}
