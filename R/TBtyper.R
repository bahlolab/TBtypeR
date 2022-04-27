
#' @export
#' @importFrom rlang is_integer is_scalar_integerish is_bool is_scalar_character
#' @importFrom treeio Nnode2
#' @importFrom magrittr set_colnames set_rownames
#' @importFrom assertthat assert_that
#' @importFrom SeqArray seqOpen
tbtype <- function(gds,
                   panel = TBtyper::tbt_panel,
                   max_phylotypes = 5L,
                   min_mix_prop = 0.001,
                   min_depth = 50L,
                   max_depth = 250L,
                   max_p_val_perm = 0.01,
                   max_p_val_wsrst = 0.01,
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
                   error_rate = 0.005,
                   verbose = FALSE) {

  # TODO: allow VCF input as well as gds input
  # TODO: Cache Allele counts?
  # TODO: Canned visualisations?

  # check args
  assert_that(
    is_gds(gds),
    is_scalar_integerish(max_phylotypes) && max_phylotypes > 0,
    is_scalar_integerish(min_depth) && min_depth > 0,
    is_scalar_integerish(max_depth) && max_depth > 0,
    is_scalar_proportion(max_p_val_perm),
    is_scalar_proportion(max_p_val_wsrst),
    is_scalar_integerish(n_perm),
    is.null(error_rate) || is_scalar_double(error_rate))

  phylo <- panel_to_phylo(panel)
  geno <- panel_to_geno(panel)
  allele_counts <- get_allele_counts(gds, panel = panel, verbose = verbose)

  # subset genotypes to those in input allele counts
  geno_sub <- geno[, colnames(allele_counts), drop =F]
  message('Note: using ', ncol(geno_sub), ' of ', ncol(geno), ' genotypes')

  # subset phylotypes to those that can be distinguished over input genotypes
  phylo_sub <- collapse_phylotypes(phylo, geno_sub)
  message('Note: using ', Nnode2(phylo_sub), ' of ', Nnode2(phylo), ' phylotypes')

  # transpose genotypes
  t_geno <- t(geno_sub)[, phylo_labels(phylo_sub), drop=F]
  rate <- rowSums(t_geno) / ncol(t_geno)

  # "permuted" genotypes background
  perm_geno <-
    vapply(seq_len(nrow(t_geno)),
           function (i) { rbinom(n_perm, 1, rate[i]) },
           integer(n_perm)) %>%
    t() %>%
    set_rownames(rownames(t_geno))

  # readr::write_tsv(tibble(event = character(), index = integer()), '.tbt_log.tsv')
  results <-
    seq_len(nrow(allele_counts)) %>%
    furrr::future_map_dfr(function(i) {
      # readr::write_tsv(tibble(event = 'start', index = i), '.tbt_log.tsv', append = T)
      tbtype_sample(
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
        n_perm = n_perm,
        exclude_parent = exclude_parent,
        exclude_child = exclude_child,
        exclude_ancestor = exclude_ancestor,
        exclude_descendant = exclude_descendant,
        exclude_distance = exclude_distance,
        exclude_inner = exclude_inner,
        fuzzy_phylotypes = fuzzy_phylotypes,
        fuzzy_max_dist = fuzzy_max_dist,
        error_rate = error_rate) %>%
        mutate(sample_id = rownames(allele_counts)[i])
      # %T>% { readr::write_tsv(tibble(event = 'end', index = i), '.tbt_log.tsv', append = T) }
    }, .options = furrr::furrr_options(seed = NULL)) %>%
    select(sample_id, tidyr::everything())

  return(results)
}

#' @importFrom tidyr replace_na chop
#' @importFrom magrittr "%<>%"
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter if_else bind_rows tibble first last
#' @importFrom purrr map_int
#' @importFrom stringr str_c
tbtype_sample <- function(phylo,
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
                          n_perm,
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
           bac = sm_allele_counts[,'alt'],
           dp = rowSums(sm_allele_counts)) %>%
    filter(dp >= min_depth) %>%
    # round down samples with depth greater than max_depth
    mutate(dp_gt_mx = dp > max_depth,
           bac = replace(bac, dp_gt_mx, round(max_depth * (bac / dp)[dp_gt_mx]) %>% as.integer()),
           dp = replace(dp, dp_gt_mx, max_depth)) %>%
    select(-dp_gt_mx) %>%
    na.omit()

  # arbitrary threshold on min sites
  if (nrow(data) < min_sites) {  return(tibble(note = "insufficient data"))  }

  phy_gts <- gts[data$variant,]
  perm_gts <- pgts[data$variant,]

  # TODO - move this up 1 level to aviod recalculating
  node_dist <- phylo_geno_dist(phylo, t(phy_gts))

  # initial search over all phylotypes
  res_1 <-
    phy_lh(data, phy_gts = phy_gts[, -rootnode(phylo)], perm_gts = perm_gts, error_rate = error_rate) %>%
    mutate(node = label_to_node(phylo, phylotype))

  res_1_top <-
    res_1 %>%
    filter(p_val_perm < max_p_val_perm) %>%
    slice(1)

  # store the best match in sample_fit
  if (nrow(res_1_top) > 0) {

    optim_error <- optimise_error_rate(data, phy_gts[, res_1_top$node])

    sample_fit <-
      res_1_top %>%
      mutate(likelihood = optim_error$likelihood, error_rate = optim_error$error_rate,
             mix_n = 1, mix_index = 1, mix_prop = 1,  #search = list(res_1),
             p_val_wsrst = -1, fuzzy_sites = list(integer()), abs_diff = Inf) %>%
      select(mix_n, mix_index, mix_prop, node, phylotype, mix_prop, likelihood,
             error_rate, p_val_perm, p_val_wsrst, abs_diff, fuzzy_sites)

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
    model_last <-
      colSums(t(mix_gts_last) * mix_prop_last) %>%
      force_to_interval()

    # find next mixture component
    nodes_exclude <- exclusions(nodes_last, phylo, node_dist,
                                exclude_parent = exclude_parent,
                                exclude_child = exclude_child,
                                exclude_ancestor = exclude_ancestor,
                                exclude_descendant = exclude_descendant,
                                exclude_distance = exclude_distance,
                                exclude_inner = exclude_inner)

    res_next <-
      phy_mix_lh(data = data,
                 phy_gts_last = mix_gts_last,
                 phy_gts_search = phy_gts[, -nodes_exclude, drop = FALSE],
                 mix_prop_last = mix_prop_last,
                 mix_prop_search = spike_in_p,
                 error_rate = optim_error$error_rate,
                 n_perm = n_perm,
                 max_p_val_perm = max_p_val_perm,
                 max_p_val_wsrst = max_p_val_wsrst) %>%
      mutate(node = label_to_node(phylo, phylotype))

    res_next_top <-
      res_next %>%
      filter(p_val_perm < max_p_val_perm,
             p_val_wsrst < max_p_val_wsrst) %>%
      slice(1)

    node_next <- res_next_top$node

    # stop if we haven't found a new candidate node
    if (length(node_next) == 0) { break }

    mix_gts <- cbind(mix_gts_last, phy_gts[, node_next])
    # optimise mixture proportions
    optim_prop <- optim_phy_mix(data, mix_gts = mix_gts,
                                error_rate = optim_error$error_rate,
                                fit = c_spike_in(mix_prop_last, spike_in_p))

    model_curr <-
      colSums(t(mix_gts) * optim_prop$prop) %>%
      force_to_interval()

    optim_error <- optimise_error_rate(data, model_curr)

    mix_fit <-
      sample_fit %>%
      filter(mix_n == i) %>%
      select(mix_index, node, fuzzy_sites) %>%
      bind_rows(tibble(mix_index = i + 1, node = node_next, fuzzy_sites = list(integer()))) %>%
      mutate(mix_n = i + 1,
             mix_prop = optim_prop$prop,
             phylotype = node_to_label(phylo, node),
             likelihood = optim_error$likelihood,
             error_rate = optim_error$error_rate,
             p_val_perm = res_next_top$p_val_perm,
             p_val_wsrst = res_next_top$p_val_wsrst,
             abs_diff = sum(abs(model_last - model_curr))) %>%
      select(mix_n, mix_index, mix_prop, node, phylotype, likelihood, error_rate,
             p_val_perm, p_val_wsrst, abs_diff, fuzzy_sites)

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

    sample_fit <- bind_rows(sample_fit, mix_fit)
  }

  return(chop(sample_fit, c(mix_prop, mix_index, node, phylotype, fuzzy_sites)))
}


optimise_error_rate <-  function(data, model,
                                 min_error = 0.0001,
                                 max_error = 0.1,
                                 min_delta_e = 0.0005,
                                 n = 10L,
                                 max_iter = 10L) {

  min_e <- min_error
  max_e <- max_error
  n_iter <- 0L
  e_last <- NA_real_
  converged <- FALSE

  while (n_iter < max_iter) {
    n_iter <- n_iter + 1L
    grid <- seq(min_e, max_e, length.out = n)
    lh <- map_dbl(grid, function(e) {
      with(data, binom_likelihood(bac, dp, model, e, by_site = F))
    })
    # message('LH=', max(lh), ' @ ', grid[which.max(lh)])
    mx <- which.max(lh)
    max_lh <- max(lh)
    e_best <- grid[mx]
    if (!is.na(e_last) && abs(e_last - e_best) < min_delta_e) {
      converged <- TRUE
      break
    }
    e_last <- e_best
    min_e <- grid[max(1L, mx - 1L)]
    max_e <- grid[min(n, mx+1L)]
  }

  return(list(error_rate = e_best, likelihood = max_lh, n_iter = n_iter, converged = converged))
}

#' @importFrom purrr map map_dbl
#' @importFrom dplyr tibble arrange desc
phy_lh <- function(data, phy_gts_search, error_rate, perm_gts = NULL) {

  # check args
  stopifnot(is.data.frame(data),
            is.matrix(phy_gts_search),
            !is.null(colnames(phy_gts_search)),
            is_scalar_proportion(error_rate),
            is.null(perm_gts) || is.matrix(perm_gts),
            is.null(perm_gts) || nrow(perm_gts) == nrow(phy_gts_search))

  # calculate sitewise marginal likelihoods for gt=0 and gt=1
  site_lh_0 <- with(data, binom_likelihood(bac, dp, 0, error_rate, by_site = TRUE))
  site_lh_1 <- with(data, binom_likelihood(bac, dp, 1, error_rate, by_site = TRUE))

  phy_lh <- colSums((phy_gts_search * site_lh_1) + ((1L - phy_gts_search) * site_lh_0))

  if (!is.null(perm_gts)) {
    perm_lh <- colSums((perm_gts * site_lh_1) + ((1L - perm_gts) * site_lh_0))
    p_val_perm <- map_dbl(phy_lh, ~ (sum(. < perm_lh, na.rm = T)) / (length(perm_lh)))
  } else {
    p_val_perm <- NA_real_
  }

  result <-
    tibble(phylotype = colnames(phy_gts_search),
           likelihood = phy_lh,
           p_val_perm = p_val_perm) %>%
    arrange(desc(likelihood), p_val_perm)

  return(result)
}

#' @importFrom purrr map map_dbl pmap_dbl map2_dbl map2
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter summarise
#' @importFrom tibble rownames_to_column
phy_mix_lh <- function(data, phy_gts_last, phy_gts_search, mix_prop_last, mix_prop_search,
                       error_rate = 0.005,
                       max_p_val_perm = 0.001,
                       max_p_val_wsrst = 0.001,
                       n_perm = 1000L,
                       score_best_only = TRUE) {

  # check args
  stopifnot(is.data.frame(data),
            is.matrix(phy_gts_last),
            is.matrix(phy_gts_search),
            nrow(data) == nrow(phy_gts_search),
            nrow(phy_gts_last) == nrow(phy_gts_search),
            !is.null(colnames(phy_gts_search)),
            is_proportion(mix_prop_last),
            is_scalar_proportion(mix_prop_search),
            is_scalar_proportion(max_p_val_perm),
            is_scalar_proportion(max_p_val_wsrst),
            is_integerish(n_perm))
  #
  site_model_last <- colSums(t(phy_gts_last) * mix_prop_last)
  site_lh_last <- with(data, binom_likelihood(bac, dp, site_model_last, error_rate, by_site = TRUE))
  lh_last <- sum(site_lh_last)
  n_sites <- length(site_lh_last)

  mix_props <-
    map_df(seq_along(mix_prop_last), function(i) {
      x <- combn(length(mix_prop_last), i)
      tibble(ix = map(seq_len(ncol(x)), function(j) x[,j]))
    }) %>%
    mutate(tot = map_dbl(ix, ~ sum(mix_prop_last[.]))) %>%
    filter(tot > mix_prop_search) %>%
    mutate(mp = map2(tot, ix, function(tot, ix) {
      mp <- c(mix_prop_last, mix_prop_search)
      mp[ix] <- (tot - mix_prop_search) / tot * mp[ix]
      mp
    })) %>%
    pull(mp) %>%
    do.call(cbind, .)

  site_lh <-
    map(setNames(c(0L, 1L), c('0', '1')), function(x) {
      vapply(seq_len(ncol(mix_props)),
             function(i) binom_likelihood(data$bac, data$dp,
                                          force_to_interval(colSums(t(cbind(phy_gts_last, x)) * mix_props[,i])),
                                          error_rate, by_site = TRUE),
             numeric(n_sites))
    })

  phy_lh <-
    vapply(seq_len(ncol(mix_props)),
           function(i) colSums((phy_gts_search * site_lh$`1`[,i]) + ((1L - phy_gts_search) * site_lh$`0`[,i])),
           numeric(ncol(phy_gts_search)))

  result <-
    set_colnames(phy_lh, seq_len(ncol(phy_lh))) %>%
    as.data.frame() %>%
    rownames_to_column('phylotype') %>%
    pivot_longer(-phylotype, names_to = 'mpi', values_to = 'likelihood') %>%
    mutate(mpi = as.integer(mpi)) %>%
    group_by(phylotype) %>%
    slice(which.max(likelihood)) %>%
    ungroup() %>%
    mutate(p_val_perm = NA_real_,
           p_val_wsrst = NA_real_) %>%
    arrange(desc(likelihood))

  site_states <- as.factor(unname(site_model_last))

  for (i in seq_along(result$phylotype)) {

    if (score_best_only & result$likelihood[i] < lh_last) { break }

    phylotype <- result$phylotype[i]
    mpi <- result$mpi[i]
    # calc p_val_perm
    perm_gt <- permute_similar(states = site_states,
                               new_gts = phy_gts_search[, phylotype],
                               n_perm = n_perm)
    perm_lh <- colSums((perm_gt * site_lh$`1`[,mpi]) + ((1L - perm_gt) * site_lh$`0`[,mpi]))
    result$p_val_perm[i] <- sum(result$likelihood[i] < perm_lh, na.rm = T) / n_perm
    # calc p_val_wsrst
    site_model <- colSums(t(cbind(phy_gts_last, phy_gts_search[, phylotype])) * mix_props[,mpi])
    sites_diff <- unname(which(site_model != site_model_last))
    site_lh_diff <- (phy_gts_search[sites_diff, phylotype] * site_lh$`1`[sites_diff, mpi]) +
      ((1L - phy_gts_search[sites_diff, phylotype]) * site_lh$`0`[sites_diff, mpi])
    result$p_val_wsrst[i] <-
      suppressWarnings(
        wilcox.test(site_lh_last[sites_diff], site_lh_diff,
                    paired = TRUE, alternative = 'less')$p.value)
    # check if passing
    if (score_best_only) {
      if (result$p_val_perm[i] > max_p_val_perm) { next }
      if (result$p_val_wsrst[i] > max_p_val_wsrst) { next }
      break
    }
  }

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
  top_fit[which.max(top_fit)] <- top_fit[which.max(top_fit)] + resolution - sum(top_fit)
  top_fit <- as.integer(top_fit)
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
  hs <- hash_set()
  n_iter <- 0L
  coeff <- 1
  converged <- FALSE
  # lh_path <- rep(NA_real_, max_iter + 1)
  # new_states <- rep(0L, max_iter + 1)
  # states_scored <- tibble()

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
      filter(!hs_contains(hs, fit_hash)) %>%
      group_by(fit_hash) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(lh = map_dbl(fit, function(fit_) {
        mix_likelihood(data_sub, mix_model(mgt_sub, fit_ / resolution), error_rate = error_rate, by_site = FALSE)
      }))
    # add hashes to search set so we don't search again
    # new_states[n_iter] <- nrow(search_2)
    hs <- hs_add(hs, search_2$fit_hash)

    # states_scored <- bind_rows(states_scored, select(search_2, fit, lh))

    search_2 <-search_2 %>% filter(lh > top_lh) %>% arrange(desc(lh))

    if (nrow(search_2) == 0) {
      # try again, reducing coeff by half
      coeff <- coeff / 2
    } else {
      # accept new best solution
      top_fit <- search_2$fit[[1]]
      top_lh <- search_2$lh[[1]]
    }
  }

  # states_scored %>%
  #   mutate(p1 = map_dbl(fit, ~ .[1]),
  #          p2 = map_dbl(fit, ~ .[2]),
  #          p3 = map_dbl(fit, ~ .[3])) %>%
  #   arrange(desc(lh)) %>%
  #   mutate(rank = row_number()) %>%
  #   ggplot(aes(p1, p3, col = rank)) +
  #   geom_point(size = 0.7)

  return(list(prop = top_fit / resolution,
              likelihood = top_lh + lh_same,
              n_iter = n_iter,
              converged = converged))
}
