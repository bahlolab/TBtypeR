
# allele_counts should be a three dimensional array
# rows == sample, cols == variant, dim 3-1 == ref_allele count, dim 3-2 == alt_allele count
# e.g allele_counts <- array (c(ref_count, alt_count), dim c (dim(ref_count), 2))
#' @export
#' @importFrom rlang is_integer is_scalar_integerish is_bool
#' @importFrom treeio Nnode2
#' @importFrom magrittr set_colnames set_rownames
fit_phylotypes <- function(allele_counts,
                           phylo = NULL,
                           geno = NULL,
                           max_phylotypes = 5L,
                           min_mix_prop = 0.005,
                           min_depth = 50L,
                           max_depth = 250L,
                           max_p_val_perm = 0.001,
                           max_p_val_wsrst = 0.001,
                           spike_in_p = 0.01,
                           reoptimise = TRUE,
                           reoptimise_max_iter = 5L,
                           min_sites = 1000L,
                           n_perm = 1000L,
                           exclude_parent = TRUE,
                           exclude_child = TRUE,
                           exclude_ancestor = TRUE,
                           exclude_descendant = FALSE,
                           exclude_distance = 10L,
                           exclude_inner = FALSE,
                           fuzzy_phylotypes = TRUE,
                           fuzzy_max_dist = 100L,
                           error_rate = 0.005
                           ) {

  if (is.null(phylo)) {
    phylo <- phylo_v1
  }

  if (is.null(geno)) {
    geno <- geno_v1
  }

  # check args
  stopifnot(
    is_integer(allele_counts) && is.array(allele_counts) && min(allele_counts, na.rm = T) >= 0L,
    length(dim(allele_counts)) == 3L && dim(allele_counts)[3] == 2L,
    is_integer(geno) && is.matrix(geno) && max(geno, na.rm = T) == 1L && min(geno, na.rm = T) == 0L,
    is_phylo(phylo),
    setequal(rownames(geno), c(phylo$tip.label, phylo$node.label)),
    all(colnames(allele_counts) %in% colnames(geno)),
    is_scalar_integerish(max_phylotypes) & max_phylotypes > 0,
    is_scalar_integerish(min_depth) & min_depth > 0,
    is_scalar_integerish(max_depth) & max_depth > 0,
    is_scalar_proportion(max_p_val_perm),
    is_scalar_proportion(max_p_val_wsrst),
    is_scalar_integerish(n_perm),
    is.null(error_rate) || is_scalar_double(error_rate))

  # subset genotypes to those in input allele counts
  geno_sub <- geno[, colnames(allele_counts)]
  message('Note: using ', ncol(geno_sub), ' of ', ncol(geno), ' genotypes')

  # subset phylotypes to those that can be distinguished over input genotypes
  phylo_sub <- collapse_phylotypes(phylo, geno_sub)
  message('Note: using ', Nnode2(phylo_sub), ' of ', Nnode2(phylo), ' phylotypes')

  # transpose genotypes
  t_geno <- t(geno_sub)[, node_to_label(phylo_sub, seq_len(Nnode2(phylo_sub)))]
  rate <- rowSums(t_geno) / ncol(t_geno)

  perm_geno <-
    vapply(seq_len(nrow(t_geno)),
           function (i) { rbinom(n_perm, 1, rate[i]) },
           integer(n_perm)) %>% t()

  readr::write_tsv(tibble(event = character(), index = integer()), '.tbt_log.tsv')
  results <-
    seq_len(nrow(allele_counts)) %>%
    furrr::future_map_dfr(function(i) {
      readr::write_tsv(tibble(event = 'start', index = i), '.tbt_log.tsv', append = T)
      fit_sample(
        phylo = phylo_sub,
        gts = t_geno,
        pgts = perm_geno,
        sm_allele_counts = allele_counts[i, ,],
        max_phylotypes = max_phylotypes,
        min_mix_prop = min_mix_prop,
        min_depth = min_depth,
        max_depth = max_depth,
        max_p_val_perm = max_p_val_perm,
        max_p_val_wsrst = max_p_val_wsrst,
        spike_in_p = spike_in_p,
        reoptimise = reoptimise,
        reoptimise_max_iter = reoptimise_max_iter,
        min_sites = min_sites,
        exclude_parent = exclude_parent,
        exclude_child = exclude_child,
        exclude_ancestor = exclude_ancestor,
        exclude_descendant = exclude_descendant,
        exclude_distance = exclude_distance,
        exclude_inner = exclude_inner,
        fuzzy_phylotypes = fuzzy_phylotypes,
        fuzzy_max_dist = fuzzy_max_dist,
        error_rate = error_rate) %>%
        mutate(sample_id = rownames(allele_counts)[i]) %T>% {
          readr::write_tsv(tibble(event = 'end', index = i), '.tbt_log.tsv', append = T)
        }
    }) %>%
    select(sample_id, tidyr::everything())

  return(results)
}

#' @importFrom tidyr replace_na chop
#' @importFrom magrittr "%<>%"
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter if_else bind_rows tibble first last
#' @importFrom purrr map_int
#' @importFrom stringr str_c
fit_sample <- function(phylo,
                       gts,
                       pgts,
                       sm_allele_counts,
                       max_phylotypes,
                       min_mix_prop,
                       min_depth,
                       max_depth,
                       max_p_val_perm,
                       max_p_val_wsrst,
                       spike_in_p,
                       reoptimise,
                       reoptimise_max_iter,
                       min_sites,
                       exclude_parent,
                       exclude_child,
                       exclude_ancestor,
                       exclude_descendant,
                       exclude_distance,
                       exclude_inner,
                       fuzzy_phylotypes,
                       fuzzy_max_dist,
                       error_rate) {

  data <-
    tibble(variant = rownames(sm_allele_counts),
           bac = sm_allele_counts[,'Alt'],
           dp = rowSums(sm_allele_counts),
           gts = gts,
           pgts = pgts) %>%
    filter(dp >= min_depth) %>%
    # round down samples with depth greater than max_depth
    mutate(dp_gt_mx = dp > max_depth,
           bac = replace(bac, dp_gt_mx, round(max_depth * (bac / dp)[dp_gt_mx]) %>% as.integer()),
           dp = replace(dp, dp_gt_mx, max_depth)) %>%
    select(-dp_gt_mx) %>%
    na.omit()

  phy_gts <- data$gts %>% set_colnames(node_to_label(phylo, seq_len(Nnode2(phylo)))) %>%
    set_rownames(data$variant)
  perm_gts <- data$pgts
  data %<>% select(-gts, -pgts)

  # arbitrary threshold on min sites
  if (nrow(data) < min_sites) {  return(tibble(note = "insufficient data"))  }

  # TODO - move this up 1 level to aviod recalculating
  node_dist <- phylo_geno_dist(phylo, t(phy_gts))

  # find the first phylotype
  res_1 <-
    score_phy_mix(data, phy_gts, perm_gts, error_rate = error_rate) %>%
    mutate(node = label_to_node(phylo, phylotype)) %>%
    arrange(desc(likelihood), p_val_perm) %>%
    filter(! node %in% exclusions(integer(), phylo, node_dist,
                                  exclude_inner = exclude_inner))
  res_1_top <-
    res_1 %>%
    filter(p_val_perm < max_p_val_perm) %>% # rho > min_rho) %>%
    slice(1)

  # store the best match in sample_fit
  if (nrow(res_1_top) > 0) {

    sample_fit <-
      res_1_top %>%
      mutate(mix_n = 1, mix_index = 1, mix_prop = 1, search = list(res_1),
             p_val_wsrst = 0, fuzzy_sites = list(integer()), abs_diff = Inf) %>%
      select(mix_n, mix_index, mix_prop, node, phylotype, mix_prop, likelihood,
             p_val_perm, p_val_wsrst, abs_diff, search, fuzzy_sites)

    if (fuzzy_phylotypes) {
      # find sites that differ in neighbouring nodes and may be incorrect in the model
      fuzz <- optim_fuzzy(data = data, node_optim = sample_fit$node[1], index_optim = 1, phy_gts = phy_gts,
                          mix_gts = phy_gts[, sample_fit$node[1], drop = FALSE], mix_prop = 1, phylo = phylo,
                          node_dist = node_dist, max_dist = fuzzy_max_dist, error_rate = error_rate)
      sample_fit[1,] %<>% mutate(likelihood = likelihood + fuzz$lh_delta,
                                 fuzzy_sites = list(fuzz$sites))
    }

  } else {
    return(tibble(note = "no matches passing filters"))
  }

  for (i in seq_len(max_phylotypes - 1)) {
    lh_last <- filter(sample_fit, mix_n == i) %>% pull(likelihood) %>% first()
    nodes_last <- filter(sample_fit, mix_n == i) %>% pull(node)
    mix_prop_last <- filter(sample_fit, mix_n == i) %>% pull(mix_prop)
    mix_gts_last <-
      filter(sample_fit, mix_n == i) %>%
      pmap(function(node, fuzzy_sites, ...) flip_at(phy_gts[, node, drop = FALSE], fuzzy_sites)) %>%
      do.call('cbind', .)

    # find next mixture component
    nodes_exclude <- exclusions(nodes_last, phylo, node_dist,
                                exclude_parent = exclude_parent,
                                exclude_child = exclude_child,
                                exclude_ancestor = exclude_ancestor,
                                exclude_descendant = exclude_descendant,
                                exclude_distance = exclude_distance,
                                exclude_inner = exclude_inner)

    if (ncol(phy_gts[, -nodes_exclude, drop = FALSE]) < 1) { break }

    res_next <-
      score_phy_mix(data = data,
                    phy_gts_search = phy_gts[, -nodes_exclude, drop = FALSE],
                    perm_gts = perm_gts,
                    error_rate = error_rate,
                    phy_gts_fixed = mix_gts_last,
                    fit = c_spike_in(mix_prop_last, spike_in_p)) %>%
      mutate(node = label_to_node(phylo, phylotype)) %>%
      arrange(desc(likelihood), p_val_perm)

    node_next <- res_next %>% filter(p_val_perm < max_p_val_perm) %>% pull(node) %>% first()
    # stop if we havent found a new candidate node
    if (is.na(node_next)) { break }
    # optimise mixture proportions
    optim_prop <- optim_phy_mix(data, mix_gts = cbind(mix_gts_last, phy_gts[, node_next]), error_rate = error_rate)
    # calculate wilcox signed rank sum p value
    # probably should allow fuzziness for this test too
    sites_diff <-
      map(seq_len(ncol(mix_gts_last)), ~ which(mix_gts_last[, .] != phy_gts[, node_next])) %>%
      unlist() %>% sort() %>% unique()
    lh0 <- mix_likelihood(data[sites_diff,],
                          mix_model(mix_gts_last[sites_diff, , drop=F],
                                    mix_prop_last),
                          error_rate = error_rate,
                          by_site = TRUE)
    lh1 <- mix_likelihood(data[sites_diff,],
                          mix_model(cbind(mix_gts_last[sites_diff, , drop=F], phy_gts[sites_diff, node_next]),
                                    optim_prop$prop),
                          error_rate = error_rate,
                          by_site = TRUE)
    p_val_wsrst <- suppressWarnings(
      tryCatch(wilcox.test(lh0, lh1, paired = TRUE, alternative = 'less')$p.value,
               error = function(e) { NA_real_ })
    )

    mix_fit <-
      sample_fit %>%
      filter(mix_n == i) %>%
      select(mix_index, node, fuzzy_sites) %>%
      bind_rows(tibble(mix_index = i + 1, node = node_next, fuzzy_sites = list(integer()))) %>%
      mutate(mix_n = i + 1,
             mix_prop = optim_prop$prop,
             search = map(node, ~ tibble(phylotype = character(),
                                         p_val_perm = numeric())),
             likelihood = optim_prop$likelihood,
             p_val_wsrst = p_val_wsrst) %>%
      select(mix_n, mix_index, mix_prop, node, likelihood, p_val_wsrst, search, fuzzy_sites)

    mix_fit$search[[i + 1]] <- res_next

    if (fuzzy_phylotypes) {
      mix_gts <-
        mix_fit %>%
        pmap(function(node, fuzzy_sites, ...) flip_at(phy_gts[, node, drop = FALSE], fuzzy_sites)) %>%
        do.call('cbind', .)
      fuzz <- optim_fuzzy(data = data, node_optim = node_next, index_optim = i + 1, phy_gts = phy_gts, mix_gts = mix_gts,
                          mix_prop = mix_fit$mix_prop, phylo = phylo, node_dist = node_dist, max_dist = fuzzy_max_dist,
                          error_rate = error_rate)
      mix_fit$likelihood %<>% { . + fuzz$lh_delta }
      mix_fit$fuzzy_sites[[i + 1]] <- fuzz$sites
    }

    # reoptimise existing mixture components in light of new component
    if (reoptimise && replace_na(p_val_wsrst < max_p_val_wsrst, FALSE)) {
      n_iter <- 0L
      n_alt <- 0L
      n_alt_last <- rep(-1L, i+1)
      n_alt_last[i+1] <- 0L
      while(n_iter < reoptimise_max_iter) {
        n_iter <- n_iter + 1L
        change_made <- FALSE
        for (j in seq_len(nrow(mix_fit))) {
          if (n_alt > n_alt_last[j]) {
            n_alt_last[j] <- n_alt
            # find highest scoring alternate components
            nodes_exclude <- exclusions(mix_fit$node[-j], phylo, node_dist,
                                        exclude_parent = exclude_parent,
                                        exclude_child = exclude_child,
                                        exclude_ancestor = exclude_ancestor,
                                        exclude_descendant = exclude_descendant,
                                        exclude_distance = exclude_distance,
                                        exclude_inner = exclude_inner)

            mix_gts <-
              mix_fit[-j, ] %>%
              pmap(function(node, fuzzy_sites, ...) flip_at(phy_gts[, node, drop = FALSE], fuzzy_sites)) %>%
              do.call('cbind', .)

            mix_fit$search[[j]] <-
              score_phy_mix(data = data,
                            phy_gts_search = phy_gts[, -nodes_exclude, drop = FALSE],
                            perm_gts = perm_gts,
                            error_rate = error_rate,
                            phy_gts_fixed = mix_gts,
                            fit = c(mix_fit$mix_prop[-j], mix_fit$mix_prop[j])) %>%
              mutate(node = label_to_node(phylo, phylotype)) %>%
              arrange(desc(likelihood), p_val_perm)

            node_next <- mix_fit$search[[j]] %>% filter(p_val_perm < max_p_val_perm) %>% pull(node) %>% first()

            if (length(node_next) == 1 && ! is.na(node_next)) {
              mix_gts <- insert_cols(mix_gts, phy_gts[, node_next], j)
              # find maximum likelihood proportion
              optim_prop <- optim_phy_mix(data, mix_gts = mix_gts, fit = mix_fit$mix_prop, error_rate = error_rate)
              lh_next <- optim_prop$likelihood

              if (fuzzy_phylotypes) {
                fuzz <- optim_fuzzy(data = data, node_optim = node_next, index_optim = j, phy_gts = phy_gts, mix_gts = mix_gts,
                                    mix_prop = optim_prop$prop, phylo = phylo, node_dist = node_dist, max_dist = fuzzy_max_dist,
                                    error_rate = error_rate)
                lh_next <- lh_next + fuzz$lh_delta
              }
              # accept is proportion is highed than current
              if (lh_next > mix_fit$likelihood[1]) {
                mix_fit$likelihood <- lh_next
                mix_fit$node[j] <- node_next
                mix_fit$mix_prop <- optim_prop$prop
                if (fuzzy_phylotypes) {
                  mix_fit$fuzzy_sites[[j]] <- fuzz$sites
                }
                change_made <- TRUE
                n_alt <- n_alt + 1L
              }
            }
          }
        }
        # exit if we couldn't find a better solution
        if (!change_made) { break }
      }
    }

    mix_fit %<>%
      mutate(phylotype = node_to_label(phylo, node),
             p_val_perm = map_dbl(search, ~ .$p_val_perm[1]))

    if (min(mix_fit$mix_prop) < min_mix_prop || replace_na(mix_fit$p_val_wsrst[1] > max_p_val_wsrst, TRUE)) { break }

    model_current <-
      mix_fit %>%
      pmap(function(node, fuzzy_sites, ...) flip_at(phy_gts[, node, drop = FALSE], fuzzy_sites)) %>%
      do.call('cbind', .) %>%
      { colSums(t(.) * mix_fit$mix_prop) }

    model_last <- colSums(t(mix_gts_last) * mix_prop_last)
    mix_fit$abs_diff <- sum(abs(model_last - model_current))

    sample_fit <- bind_rows(sample_fit, mix_fit)
  }

  return(chop(sample_fit, c(mix_prop, mix_index, node, phylotype, p_val_perm, search, fuzzy_sites)))
}

#' @importFrom purrr map map_dbl pmap_dbl map2_dbl
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter summarise
score_phy_mix <- function(data, phy_gts_search, perm_gts, error_rate,
                          phy_gts_fixed = matrix(0L, nrow(phy_gts_search), 0L),
                          fit = rep(1/(ncol(phy_gts_fixed) + 1), ncol(phy_gts_fixed) + 1)) {

  # check args
  stopifnot(is.data.frame(data),
            is.matrix(phy_gts_search),
            is.matrix(phy_gts_fixed),
            is.matrix(perm_gts),
            !is.null(colnames(phy_gts_search)),
            nrow(data) == nrow(phy_gts_search),
            nrow(data) == nrow(phy_gts_fixed),
            nrow(data) == nrow(perm_gts),
            is_double(fit) && all(fit >= 0 | fit <= 1) && (abs(1 - sum(fit)) < 2e-10),
            length(fit) == 1 + ncol(phy_gts_fixed))

  # calculate sitewise marginal likelihoods for gt=0 and gt=1
  if (length(fit) == 1L) {
    site_model_0 <- numeric(nrow(data))
    site_model_1 <- numeric(nrow(data))
    site_model_1[] <- 1
  } else {
    site_model_0 <-
      cbind(phy_gts_fixed, 0) %>% { colSums(t(.) * c(fit)) } %>%
      force_to_interval()
    site_model_1 <-
      cbind(phy_gts_fixed, 1) %>% { colSums(t(.) * c(fit)) } %>%
      force_to_interval()
  }
  site_lh_0 <- with(data, binom_likelihood(bac, dp, site_model_0, error_rate, by_site = TRUE))
  site_lh_1 <- with(data, binom_likelihood(bac, dp, site_model_1, error_rate, by_site = TRUE))

  # permutation test likelihoods
  perm_lhs <- map_dbl(seq_len(ncol(perm_gts)),
                      ~ sum(if_else(perm_gts[, .] == 1, site_lh_1, site_lh_0)))

  # score likelihoods for all search phylotypes
  result <-
    tibble(index = seq_len(ncol(phy_gts_search)),
           phylotype = colnames(phy_gts_search)) %>%
    mutate(lh = map_dbl(index, ~ sum(if_else(phy_gts_search[, .] == 1, site_lh_1, site_lh_0))),
           p_val = map_dbl(lh, ~ (sum(. < perm_lhs, na.rm = T)) / (length(perm_lhs)))) %>%
    arrange(desc(lh), p_val) %>%
    select(phylotype, likelihood = lh, p_val_perm = p_val)

  return(result)
}

# optim_fuzzy <- function(data, phy_gts, phylo, node_dist, max_dist, mix_nodes, mix_fit, index_optim = 1L)
optim_fuzzy <- function(data, node_optim, index_optim, phy_gts, mix_gts, mix_prop, phylo, node_dist, error_rate, max_dist) {

  # check args
  stopifnot(is.data.frame(data),
            is_phylo(phylo),
            is.matrix(phy_gts) && ncol(phy_gts) > 0,
            is.matrix(mix_gts) && ncol(mix_gts) > 0,
            nrow(data) == nrow(phy_gts),
            ncol(phy_gts) == Nnode2(phylo),
            is.matrix(node_dist) && ncol(node_dist) == nrow(node_dist),
            ncol(node_dist) == ncol(phy_gts),
            is_double(mix_prop) && all(mix_prop >= 0 | mix_prop <= 1) && (abs(1 - sum(mix_prop)) < 2e-10),
            ncol(mix_gts) == length(mix_prop),
            is_scalar_integerish(max_dist) && max_dist >= 0,
            is_scalar_integerish(index_optim) && index_optim > 0 && index_optim <= length(mix_prop),
            is_scalar_integerish(node_optim) && node_optim > 0 && node_optim <= ncol(phy_gts))

  # null result
  result <- list(sites = integer(),
                 optimised = FALSE,
                 neighbour = NA_integer_,
                 n_sites = 0,
                 lh_delta = 0)

  # find neighbouring nodes
  neighb_nodes <-
    node_dist_range(node = node_optim,
                    phylo = phylo,
                    node_dist = node_dist,
                    max_dist = max_dist,
                    inclusive = TRUE)

  if (length(neighb_nodes) > 0) {
    mix_gts <- cbind(phy_gts[, node_optim], mix_gts[, -index_optim, drop = FALSE])
    # find sites that differ from neighbours and flipping increases likelihood
    optim_sites <-
      optim_fuzzy_search(data,
                         mix_gts = mix_gts,
                         gts_neighb = phy_gts[, neighb_nodes, drop = FALSE],
                         dist_neighb = node_dist[neighb_nodes],
                         prop = c(mix_prop[index_optim], mix_prop[-index_optim]),
                         max_sites = max_dist,
                         error_rate = error_rate)

    if (length(optim_sites$n_index) == 1) {
      result <- list(sites = optim_sites$site_index,
                     optimised = TRUE,
                     neighbour = neighb_nodes[optim_sites$n_index],
                     n_sites = length(optim_sites$site_index),
                     lh_delta = optim_sites$lh_delta)
    }
  }

  return(result)
}

#' @importFrom dplyr count mutate filter bind_cols select starts_with
#' @importFrom tidyr pivot_longer
optim_fuzzy_search <- function(data, mix_gts, gts_neighb, dist_neighb, prop, max_sites, error_rate, fuzzy_conf = 1e-10) {
  # first column of mix_gts corresponds to node to be optimised
  # check args
  stopifnot(is.data.frame(data),
            is.matrix(mix_gts) && ncol(mix_gts) > 0,
            is.matrix(gts_neighb) && ncol(gts_neighb) > 0,
            nrow(data) == nrow(mix_gts),
            nrow(data) == nrow(gts_neighb),
            is_double(prop) && all(prop >= 0 | prop <= 1) && (abs(1 - sum(prop)) < 2e-10),
            is_scalar_integerish(max_sites))

  sites_optim <- rowSums(cbind(mix_gts[,1], gts_neighb)) %>% { which(. > 0 & . < ncol(gts_neighb) + 1) }
  data_sub <- data[sites_optim, ]
  mix_gts <- mix_gts[sites_optim, , drop =F]
  gts_neighb <- gts_neighb[sites_optim, , drop = FALSE]

  # find neighbour with highest likelihood delta after flipping up to max_dist sites
  # break ties by going with lower distance
  result <-
    tibble(site_index = sites_optim,
         gt_0 = mix_gts[, 1],
         site_lh_0 = {
           cbind(0, mix_gts[, -1, drop = FALSE]) %>%
             { colSums(t(.) * c(prop)) } %>%
             force_to_interval() %>%
             { with(data_sub, binom_likelihood(bac, dp, ., error_rate, by_site = TRUE)) }
         },
         site_lh_1 = {
           cbind(1, mix_gts[, -1, drop = FALSE])  %>%
             { colSums(t(.) * c(prop)) } %>%
             force_to_interval() %>%
             { with(data_sub, binom_likelihood(bac, dp, ., error_rate, by_site = TRUE)) }
         }) %>%
    mutate(site_lh_0 = if_else(gt_0 == 0, site_lh_0, site_lh_0 + log(fuzzy_conf)),
           site_lh_1 = if_else(gt_0 == 1, site_lh_1, site_lh_1 + log(fuzzy_conf)),
           gt_opt = if_else(site_lh_0 > site_lh_1, 0, 1)) %>%
    bind_cols(
      set_colnames(gts_neighb, str_c('N', seq_len(ncol(gts_neighb)))) %>% as_tibble()) %>%
    filter(gt_0 != gt_opt) %>%
    mutate(lh_delta = abs(site_lh_0 - site_lh_1)) %>%
    select(-gt_0, -site_lh_0, -site_lh_1) %>%
    pivot_longer(cols = starts_with('N'),
                 names_to = 'n_index',
                 names_prefix = 'N',
                 values_to = 'gt_n') %>%
    mutate(n_index = as.integer(n_index)) %>%
    filter(gt_opt == gt_n) %>%
    group_by(n_index) %>%
    arrange(desc(lh_delta)) %>%
    slice(seq_len(max_sites)) %>%
    summarise(site_index = list(site_index),
              num_alt = n(),
              lh_delta = sum(lh_delta),
              dist_n = dist_neighb[n_index[1]]) %>%
    arrange(desc(lh_delta), dist_n, n_index) %>%
    slice(1) %>%
    { `if`(nrow(.) > 0,
           with(., list(n_index = n_index, site_index = site_index[[1]], lh_delta = lh_delta)),
           list(n_index = integer(), site_index = integer(), lh_delta = 0))
    }

  return(result)
}

# Implements a binary search to optimise phylotype mixture
#' @importFrom purrr pmap map_dbl map_chr map_df map2_lgl
#' @importFrom tidyr expand_grid
#' @importFrom dplyr arrange_all transmute
#' @importFrom magrittr "%T>%"
optim_phy_mix <- function(data,
                          mix_gts,
                          error_rate,
                          fit = NULL,
                          resolution = 10000L,
                          max_iter = 1000L) {

  stopifnot(is.data.frame(data),
            is.matrix(mix_gts),
            nrow(mix_gts) == nrow(data),
            ncol(mix_gts) > 1,
            is.null(fit) || is_double(fit) && length(fit) == ncol(mix_gts),
            is_scalar_integerish(resolution) && resolution > ncol(mix_gts),
            is_scalar_integerish(max_iter))

  if (is.null(fit)) {
    fit <- rep(1 / ncol(mix_gts), ncol(mix_gts))
  }

  # sites that differ between mixture componenets
  mix_n <- ncol(mix_gts)
  sites_diff <- rowSums(mix_gts) %>% { which(. > 0 & . < mix_n) }
  data_sub <- data[sites_diff, ]
  mgt_sub <- mix_gts[sites_diff, , drop=FALSE]
  t_m_gt <- t(mix_gts[sites_diff, , drop=FALSE])

  # likelihood at sites that don't differ in mixture
  lh_same <-
    `if`(length(sites_diff) < nrow(data),
         mix_likelihood(data[-sites_diff, ], mix_gts[-sites_diff, 1], error_rate = error_rate, by_site = FALSE),
         0)

  top_fit <- as.integer(fit * resolution)
  resolution <- sum(top_fit)
  top_lh <- mix_likelihood(data_sub, mix_model(mgt_sub, top_fit / resolution), error_rate = error_rate, by_site = FALSE)

  # set of possible directions to search
  search_0 <-
    as.character(seq_len(mix_n)) %>%
    setNames(., .) %>%
    map(~ c(TRUE, FALSE)) %>%
    do.call('expand_grid', .) %>%
    mutate(state = seq_len(n())) %>%
    gather(-state, key = 'node', value = 'is_on') %>%
    arrange(state, node) %>%
    group_by(state) %>%
    filter(any(is_on)) %>%
    summarise(nodes = list(which(is_on))) %>%
    with(expand_grid(nodes_1 = nodes, nodes_2 = nodes)) %>%
    filter(purrr::map2_lgl(nodes_1, nodes_2, ~ length(intersect(.x, .y)) == 0)) %>%
    mutate(n_1 = lengths(nodes_1),
           n_2 = lengths(nodes_2),
           lcm = pracma::Lcm(n_1, n_2))

  # memoise likelihood computations using a hasmap
  hashes <- character()
  n_iter <- 0L
  coeff <- 1
  converged <- FALSE
  # lh_path <- rep(NA_real_, max_iter + 1)
  # new_states <- rep(0L, max_iter + 1)

  while(n_iter < max_iter) {
    n_iter <- n_iter + 1L
    # lh_path[n_iter] <- top_lh

    search_1 <-
      search_0 %>%
      mutate(max_delta = coeff * map_int(nodes_1, ~ min(top_fit[.]) * length(.)),
             delta = as.integer(max_delta - (max_delta %% lcm))) %>%
      filter(delta > 0)

    if (nrow(search_1) == 0) {
      # no remaining states to test
      converged <- TRUE
      break
    }

    search_2 <-
      search_1 %>%
      transmute(
        fit = pmap(., function(nodes_1, nodes_2, n_1, n_2, delta, ...) {
          x <- top_fit
          x[nodes_1] %<>% { . - delta / n_1 }
          x[nodes_2] %<>% { . + delta / n_2 }
          return(as.integer(x)) }),
        fit_hash = map_chr(fit, ~ digest::digest(.))) %>%
      filter(!fit_hash %in% hashes) %>%
      group_by(fit_hash) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(lh = map_dbl(fit, function(fit_) {
        mix_likelihood(data_sub, mix_model(mgt_sub, fit_ / resolution), error_rate = error_rate, by_site = FALSE)
      }))
    # add hashes to search set so we don't search again
    # new_states[n_iter] <- nrow(search_2)
    hashes <- union(hashes, search_2$fit_hash)

    search_2 %<>% filter(lh > top_lh) %>% arrange(desc(lh))

    if (nrow(search_2) == 0) {
      # try again, reducing coeff by half
      coeff <- coeff / 2
    } else {
      # accept new best solution, reset coeff
      top_fit <- search_2$fit[[1]]
      top_lh <- search_2$lh[[1]]
      coeff <- 1
    }
  }
  return(list(prop = top_fit / resolution,
              likelihood = top_lh + lh_same,
              n_iter = n_iter,
              converged = converged))
}


