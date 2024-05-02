
#' @export
#' @importFrom dplyr filter select mutate case_when if_else
#' @importFrom tidyr nest chop unnest expand_grid
#' @importFrom purrr map_lgl
#' @importFrom SeqArray seqOpen seqClose
assign_dr <- function(tbtype_results, gds,
                      dr_panel = TBtypeR::who_mut_cat_2021,
                      verbose = FALSE,
                      min_allele_count = 3,
                      min_allele_freq = 0.01,
                      min_depth = 10L,
                      conf_int = 0.95,
                      min_posterior = 0.90) {

  # TODO: check tbtype results
  assert_that(
    is.data.frame(tbtype_results),
    is_gds(gds),
    is.null(dr_panel) || check_panel(dr_panel, "dr"),
    is_scalar_integerish(min_allele_count) && min_allele_count > 0,
    is_scalar_proportion(min_allele_freq),
    is_scalar_proportion(conf_int),
    is_scalar_proportion(min_posterior)
  )

  drug_set <-
    dr_panel %>%
    select(drugs) %>%
    distinct() %>%
    separate_rows(drugs, sep = ';') %>%
    distinct() %>%
    pull()

  status_levels <- c('SENS', 'MIX-HET', 'HET-RES', 'HOM-RES')

  if (!is_open_gds(gds)) {
    message("Opening gds file")
    gds <- seqOpen(gds$filename, allow.duplicate = T)
    on.exit({
      seqClose(gds)
    })
  }

  SeqArray::seqSetFilter(gds, sample.id = unique(tbtype_results$sample_id), verbose = verbose)

  dr_res <-
    get_allele_counts(gds,
                      panel = dr_panel,
                      as_tibble = TRUE,
                      verbose = verbose,
                      max_ext_freq = 1
    ) %>%
    filter(
      alt_ac >= min_allele_count,
      baf >= min_allele_freq,
      depth >= min_depth
    ) %>%
    inner_join(panel_with_vid(dr_panel) %>% select(vid, drugs),
               by = "vid"
    ) %>%
    select(sample_id, vid, alt_ac, depth, baf, drugs) %>%
    nest(drug_ac = c(-sample_id)) %>%
    inner_join(tbtype_results %>%
                 filter(!is.na(error_rate)) %>%
                 unnest_mixtures(warn=FALSE) %>%
                 select(sample_id, n_phy, error_rate, mix_phylotype, mix_prop) %>%
                 nest(type_data = -c(sample_id, n_phy)),
               by = "sample_id"
    ) %>%
    mutate(drug_assign = map2(drug_ac, type_data, function(drug_ac, type_data) {
      assert_that(sum(type_data$mix_prop) == 1)
      err <- first(type_data$error_rate)
      mix_n <- nrow(type_data)
      x1 <-
        tibble(
        mix_phylotype = combinations(type_data$mix_phylotype),
        phy_set = seq_along(mix_phylotype)
      ) %>%
        unnest(mix_phylotype) %>%
        left_join(select(type_data, mix_phylotype, mix_prop),
                  by = "mix_phylotype"
        ) %>%
        chop(-phy_set) %>%
        mutate(
          mix_prop = map_dbl(mix_prop, sum, na.rm = T),
          p = mix_prop * (1 - err) + (1 - mix_prop) * (err)
        ) %>%
        expand_grid(drug_ac, .) %>%
        rowwise() %>%
        mutate(binom_prob = binom.test(alt_ac, depth, p)$p.value) %>%
        group_by(vid, alt_ac, depth, baf, drugs) %>%
        mutate(posterior = binom_prob / sum(binom_prob)) %>%
        slice(which.max(binom_prob)) %>%
        select(-phy_set, -mix_prop, -p) %>%
        mutate(status = case_when(
          binom_prob >= (1-conf_int) & posterior >= min_posterior ~ 'HOM-RES',
          mix_n == 1                                              ~ 'HET-RES',
          TRUE                                                    ~ 'MIX-HET'
        )) %>%
        # remove cases where error is most likely
        filter(!(status == 'HOM-RES' &&
                   lengths(mix_phylotype) == 1 &&
                   is.na(mix_phylotype[[1]]))) %>%
        ungroup() %>%
        mutate(mix_phylotype = if_else(status == "MIX-HET",
                                       list(type_data$mix_phylotype),
                                       as.list(mix_phylotype)),
               binom_prob = if_else(status == "MIX-HET", NA_real_, binom_prob),
               posterior  = if_else(status == "MIX-HET", NA_real_, posterior))
    }))

  if (nrow(dr_res)) {
    dr_res <- dr_res %>%
      select(sample_id, n_phy, drug_assign) %>%
      unnest(drug_assign) %>%
      unnest(mix_phylotype) %>%
      separate_rows(drugs, sep = ';') %>%
      rename(drug = drugs) %>%
      mutate(variant_status = status) %>%
      group_by(sample_id, n_phy, mix_phylotype, drug) %>%
      mutate(status = min(status)) %>%
      ungroup() %>%
      nest(variant_data = -c(sample_id, n_phy, mix_phylotype, drug, status)) %>%
      nest(mix_drug_res = -c(sample_id, n_phy, mix_phylotype))
  } else {
    dr_res <- tibble(sample_id = character(),
                     n_phy = double(),
                     mix_phylotype = character(),
                     mix_drug_res = list())
  }

  ret <-
    tbtype_results %>%
    unnest_mixtures(warn=FALSE) %>%
    left_join(dr_res, by = c('sample_id','n_phy', 'mix_phylotype')) %>%
    mutate(mix_drug_res = map(mix_drug_res, function(x) {
      `if`(is.null(x),
           expand_grid(drug = drug_set,
                       status = ordered('SENS', status_levels),
                       variant_data = list(tibble(
                         vid            = character(),
                         alt_ac         = integer(),
                         depth          = integer(),
                         baf            = double(),
                         binom_prob     = double(),
                         posterior      = double(),
                         variant_status = character(),
                       ))),
           complete(x, drug = drug_set) %>%
             mutate(status = if_else(is.na(status), 'SENS', status) %>%
                      ordered(status_levels)))
    })) %>%
    nest_mixtures(warn=FALSE) %>%
    mutate(drugs_resistant = map(mix_drug_res, function(x) {
      bind_rows(x) %>%
        filter(status != 'SENS') %>%
        pull(drug) %>%
        unique() %>%
        na.omit() %>%
        sort()
    }))

  return(ret)
}

