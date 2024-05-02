
#' @export
plot_phylotype_distribution <- function(tbt_res,
                                        phy_pal = NULL) {

  if (is.null(phy_pal)) {
    phy_pal <- tbt_palette(unique(unlist(tbt_res$phylotype)), sort=T)
  }

  tbt_res %>%
    unnest_mixtures(warn=F) %>%
    select(phylotype = mix_phylotype, mix_prop) %>%
    group_by(phylotype) %>%
    summarise(mix_prop = sum(mix_prop, na.rm = T)) %>%
    arrange(mix_prop) %>%
    mutate(phylotype = forcats::as_factor(phylotype)) %>%
    ggplot(aes(phylotype, mix_prop, fill = phylotype)) +
    geom_col() +
    scale_fill_manual(values = phy_pal) +
    coord_flip() +
    geom_text(aes(label = phylotype, y = 0, hjust = 0), nudge_y = 0.02) +
    guides(fill = 'none') +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text.y.left = element_text(angle = 0),
          panel.grid.minor.y = element_blank()) +
    labs(y = 'total')
}

#' @export
plot_phylotype_heatmap <- function(tbt_res, snp_dist,
                                   phy_pal = NULL,
                                   log_scale = FALSE) {
  assert_that(
    is.data.frame(tbt_res),
    inherits(snp_dist, 'dist'),
    is_bool(log_scale)
  )

  require(ggplot2)

  # TODO: check inputs

  if (log_scale) {
    snp_dist <- log1p(snp_dist)
  }

  levels <- intersect(labels(snp_dist),with(tbt_res, sample_id[!failed]))
  dist_matrix <- as.matrix(snp_dist)[levels, levels]
  clustering <- hclust(as.dist(dist_matrix), method = 'average')
  levels <- with(clustering, labels[order])



  treep <-
    ggtree::ggtree(clustering) +
    coord_flip() +
    scale_x_reverse() +
    theme(plot.margin = unit(c(0,0,0,0), units = "lines" ))

  phyp_data <-
    tbt_res %>%
    unnest_mixtures(warn = FALSE) %>%
    filter(sample_id %in% levels) %>%
    select(sample_id, phylotype = mix_phylotype, mix_prop) %>%
    mutate(sample_id = factor(sample_id, levels)) %>%
    arrange(sample_id) %>%
    mutate(rank = row_number()) %>%
    group_by(phylotype) %>%
    mutate(rank = mean(rank)) %>%
    ungroup() %>%
    arrange(rank, sample_id) %>%
    mutate(phylotype = as.factor(phylotype))

  if (is.null(phy_pal)) {
    phy_pal <- tbt_palette(levels(phyp_data$phylotype), sort=F)
  }

  phyp <-
    ggplot(phyp_data) +
    geom_col(aes(sample_id, mix_prop, fill = phylotype)) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = unit(c(0,0,0,0), units = "lines" )) +
    scale_fill_manual(values = phy_pal)

  hm <-
    dist_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column('sample_id') %>%
    tidyr::pivot_longer(-sample_id,
                        names_to = 'sample_id_2',
                        values_to = 'dist') %>%
    mutate(sample_id = factor(sample_id, levels),
           sample_id_2 = factor(sample_id_2, rev(levels))) %>%
      ggplot(aes(sample_id, sample_id_2, fill = dist)) +
    geom_tile() +
    scale_fill_viridis_c(option = 'E') +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = unit(c(0,0,0,0), units = "lines" ))


  plots <-
    cowplot::plot_grid(treep, phyp + guides(fill = 'none'), hm + guides(fill = 'none'),
                       ncol = 1, align = 'v', axis='lr', rel_heights = c(1,1,8))

  phyp_leg <- cowplot::get_legend(phyp)
  hm_leg <- cowplot::get_legend(hm)

  legends <-
    cowplot::plot_grid(phyp_leg, hm_leg, ncol = 1, rel_heights = c(3,1))

  cowplot::plot_grid(plots, legends, nrow = 1, rel_widths = c(10, 2))
}


#' @export
plot_sample_comparison <- function(tbt_dr_res,
                                   snp_dist,
                                   sample_order = NULL,
                                   phy_pal = NULL,
                                   MAX_DIST = NULL,
                                   title = NULL) {

  assert_that(
    is.data.frame(tbt_dr_res),
    inherits(snp_dist, 'dist'),
    is.null(sample_order) || is.character(sample_order),
    is.null(phy_pal) ||
      (is.character(phy_pal) &&
         !is.null(names(phy_pal)) &&
         all(unlist(tbt_dr_res$mix_phylotype) %in% names(phy_pal))),
    all(tbt_dr_res$sample_id %in% labels(snp_dist))
  )

  if (is.null(sample_order)) {
    sample_order <- unique(tbt_dr_res$sample_id) %>% sort()
  }

  data <-
    tbt_dr_res %>%
    unnest_mixtures(warn = F) %>%
    select(sample_id, phylotype = mix_phylotype, mix_prop, mix_drug_res) %>%
    mutate(sample_id = factor(sample_id, sample_order)) %>%
    unnest(mix_drug_res) %>%
    select(-variant_data) %>%
    nest(drug_data = c(drug, status)) %>%
    (function(x) {
      unique(x$sample_id) %>%
        as.character() %>%
        { as.matrix(snp_dist)[.,.] } %>%
        as.data.frame() %>%
        rownames_to_column('sample_id') %>%
        pivot_longer(-sample_id, names_to = 'sample_id_2', values_to = 'dist') %>%
        mutate(sample_id = factor(sample_id, sample_order),
               sample_id_2 = factor(sample_id_2, sample_order)) %>%
        { left_join(x, ., by = 'sample_id')} %>%
        nest(dist_data = c(sample_id_2, dist))
    }) %>%
    arrange(sample_id, desc(mix_prop)) %>%
    group_by(sample_id) %>%
    mutate(phy_ord = seq_along(phylotype)) %>%
    ungroup() %>%
    mutate(phy_ord = as.ordered(phy_ord))

  if (is.null(phy_pal)) {
    phy_pal <- tbt_palette(data$phylotype)
  }

  if (is.null(MAX_DIST)) {
    MAX_DIST <- data$dist_data %>% map('dist') %>% unlist() %>% max(na.rm=T)
  }


  require(ggplot2)

  p1 <-
    data %>%
    select(sample_id, phylotype, phy_ord, mix_prop) %>%
    distinct() %>%
    ggplot(aes(phy_ord, mix_prop)) +
    geom_col(aes(fill = phylotype), width = .85, col = 'black', size = 0.25) +
    geom_text(aes(label = phylotype, y = 0, hjust = 0), nudge_y = 0.02) +
    coord_flip() +
    scale_fill_manual(values = phy_pal) +
    facet_grid(sample_id ~ .,  scales = 'free_y', space = 'free_y', switch = 'y') +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle=45, vjust = 1, hjust=1),
          strip.text.y.left = element_text(angle = 0),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank()) +
    guides(fill = 'none') +
    labs(y = 'mix. prop.')

  if (!is.null(title)) {
    p1 <- p1 + ggtitle(title)
  }

  p2 <-
    data %>%
    select(sample_id, phy_ord, dist_data) %>%
    unnest(dist_data) %>%
    distinct() %>%
    mutate(dist = if_else(phy_ord == 1, dist, NA_real_)) %>%
    mutate(lab = signif(dist, digits = 2)) %>%
    ggplot(aes(phy_ord, sample_id_2, fill = dist)) +
    geom_tile(width = 0.85, col = 'black', size = 0.25) +
    ggfittext::geom_fit_text(aes(label = lab), size = 8, min.size = 2, contrast = T) +
    coord_flip() +
    facet_grid(sample_id ~ .,  scales = 'free_y', space = 'free_y', switch = 'y') +
    guides(fill = 'none') +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle=45, vjust = 1, hjust=1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.text.y = element_blank() ,
          strip.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), units = "lines" )) +
    scale_fill_viridis_c(trans = 'log1p',
                         option = 'A',
                         begin = 0.15,
                         limits = c(0, MAX_DIST),
                         na.value = NA) +
    labs(y = 'SNP dist.')

  p3 <-
    data %>%
    select(sample_id, phy_ord, drug_data) %>%
    unnest(drug_data) %>%
    ggplot(aes(phy_ord, drug)) +
    geom_tile(aes(fill = status), col = 'black', size = 0.25, width = 0.85) +
    coord_flip() +
    scale_fill_ordinal(na.value = 'gray50', drop = FALSE) +
    facet_grid(sample_id ~ .,  scales = 'free_y', space = 'free_y', switch = 'y') +
    theme(axis.title.y = element_blank(),
          axis.text.x = element_text(angle=45, vjust = 1, hjust=1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.y = element_blank(),
          strip.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), units = "lines" ),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.title = element_blank(),
          legend.position = 'top',
          legend.key.size = unit(4, 'mm')) +
    labs(y = 'pred. res.')

  NMP <- 10
  NSM <- dplyr::n_distinct(data$sample_id) + 2
  NDR <- (map(data$drug_data, 'drug') %>% unlist() %>% dplyr::n_distinct()) + 2
  cowplot::plot_grid(p1, p2, p3,
                     nrow = 1,
                     rel_widths = c(NMP, NSM, NDR),
                     align='h',
                     axis = 'tb')
}

#' @export
tbt_palette <- function(x, sort=TRUE) {
  unique(x) %>%
    { `if`(sort, sort(.), .) } %>%
    { setNames(hue_tonal(length(.)), .) }
}


hue_tonal <- function(n, l1=70, l2=40, h=c(0, 360)+15, phase=1) {
  h1 <- scales::hue_pal(h=h, l=l1)(n)
  h2 <- scales::hue_pal(h=h, l=l2)(n)
  r <- character(n)
  if(phase == 1){
    r[seq.int(1, n, 2)] <- h1[seq.int(1, n, 2)]
    if (n > 1){
      r[seq.int(2, n, 2)] <- h2[seq.int(2, n, 2)]
    }
  } else {
    r[seq.int(1, n, 2)] <- h1[seq.int(2, n, 2)]
    if (n > 1) {
      r[seq.int(2, n, 2)] <- h2[seq.int(1, n, 2)]
    }
  }
  r
}

