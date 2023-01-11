
#' @export
#' @importFrom dplyr filter select mutate case_when if_else
#' @importFrom tidyr nest chop unnest expand_grid
#' @importFrom purrr map_lgl
#' @importFrom SeqArray seqOpen seqClose
assign_dr <- function(tbtype_results, gds,
                      dr_panel = TBtyper::who_dr_panel,
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
          binom_prob >= (1-conf_int) & posterior >= min_posterior ~ 'homoresistant',
          mix_n == 1                                              ~ 'heteroresistant',
          TRUE                                                    ~ 'uncertain'
        )) %>%
        # remove cases where error is most likely
        filter(!(status == 'homoresistant' &&
                   lengths(mix_phylotype) == 1 &&
                   is.na(mix_phylotype[[1]]))) %>%
        ungroup() %>%
        mutate(mix_phylotype = if_else(status == "uncertain",
                                       list(type_data$mix_phylotype),
                                       as.list(mix_phylotype)),
               binom_prob = if_else(status == "uncertain", NA_real_, binom_prob),
               posterior  = if_else(status == "uncertain", NA_real_, posterior))
    })) %>%
    select(sample_id, n_phy, drug_assign) %>%
    unnest(drug_assign) %>%
    mutate(status = ordered(status, c('homoresistant', 'heteroresistant', 'uncertain'))) %>%
    unnest(mix_phylotype) %>%
    separate_rows(drugs, sep = ';') %>%
    rename(drug = drugs) %>%
    mutate(variant_status = status) %>%
    group_by(sample_id, n_phy, mix_phylotype, drug) %>%
    mutate(status = min(status)) %>%
    ungroup() %>%
    nest(variant_data = -c(sample_id, n_phy, mix_phylotype, drug, status)) %>%
    nest(mix_drug_res = -c(sample_id, n_phy, mix_phylotype))

  ret <-
    tbtype_results %>%
    unnest_mixtures(warn=FALSE) %>%
    left_join(dr_res, by = c('sample_id','n_phy', 'mix_phylotype')) %>%
    mutate(mix_drug_res = map(mix_drug_res, function(x) {
      `if`(is.null(x), tibble(drug = NA_character_, status = NA_character_), x)
    })) %>%
    nest_mixtures(warn=FALSE) %>%
    mutate(drugs_resistant = map(mix_drug_res, function(x) {
      bind_rows(x) %>% pull(drug) %>% unique() %>% na.omit() %>% sort()
    }))

  return(ret)
}

