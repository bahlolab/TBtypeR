
#' @export
plot_phylotype_heatmap <- function(tbt_res, snp_dist,
                                   log_scale = FALSE) {
  assert_that(
    is.data.frame(tbtype_results),
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

  phyp <-
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
    mutate(phylotype = as.factor(phylotype)) %>%
    (function(x) {
      pal <-
        levels(x$phylotype) %>%
        (function(y) setNames(hue_tonal(length(y), h = c(10, 350)), y))
      ggplot(x) +
        geom_col(aes(sample_id, mix_prop, fill = phylotype)) +
        theme(axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              plot.margin = unit(c(0,0,0,0), units = "lines" )) +
        scale_fill_manual(values = pal)
    })

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

plot_drug_res_heatmap <- function(tbt_dr_res,
                                  snp_dist,
                                  metadata = 'sample_id',
                                  order = NULL,
                                  pal = NULL) {

  # tbt_dr_res <- readRDS('~/scratch/tbt_dr_res.rds')
  # snp_dist <- readRDS('~/scratch/snp_dist.rds')
  # metadata <- c('culture_date', 'set', 'sample_id')

  assert_that(
    is.data.frame(tbt_dr_res),
    inherits(snp_dist, 'dist'),
    is.null(order) || is.character(order),
    is.character(metadata) && length(metadata) > 0 && all(metadata %in% colnames(tbt_dr_res)),
    is.null(pal) || (is.character(pal) && length(pal) > 0 && !is.null(names(pal)))
  )

  require(ggplot2)

  phyp <-
    tbt_dr_res %>%
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
    mutate(phylotype = as.factor(phylotype)) %>%
    (function(x) {
      if (is.null(pal)) {
        pal <-
          levels(x$phylotype) %>%
          (function(y) setNames(hue_tonal(length(y), h = c(10, 350)), y))
      }
      ggplot(x) +
        geom_col(aes(sample_id, mix_prop, fill = phylotype)) +
        theme(axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              plot.margin = unit(c(0,0,0,0), units = "lines" )) +
        scale_fill_manual(values = pal) +
        coord_flip()
    })



}

hue_tonal <- function(n, l1=70, l2=40, h=c(0, 360), phase=1) {
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


