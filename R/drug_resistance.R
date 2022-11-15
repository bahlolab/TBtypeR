
#' @export
#' @importFrom dplyr filter select mutate case_when if_else
#' @importFrom tidyr nest chop unnest expand_grid
#' @importFrom purrr map_lgl
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
    warning("Opening closed gds file")
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
                 unnest_mixtures() %>%
                 select(sample_id, error_rate, mix_phylotype, mix_prop) %>%
                 nest(type_data = c(-sample_id)),
               by = "sample_id"
    ) %>%
    mutate(drug_assign = map2(drug_ac, type_data, function(drug_ac, type_data) {
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
        # expand_grid(drug_ac %>% mutate(alt_ac = replace(alt_ac, 1, as.integer(252 * 0.30))), .) %>%
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
        mutate(mix_phylotype = as.list(mix_phylotype),
               binom_prob = if_else(status == "uncertain", NA_real_, binom_prob),
               posterior  = if_else(status == "uncertain", NA_real_, posterior))
    })) %>%
    select(sample_id, drug_assign) %>%
    unnest(drug_assign)


  ret <-
    tbtype_results %>%
    unnest_mixtures() %>%
    left_join(
      dr_res %>%
        filter(status != 'homoresistant') %>%
        select(-mix_phylotype) %>%
        nest(homo_res = -sample_id),
      by = 'sample_id'
    ) %>%
    left_join(
      dr_res %>%
        filter(status == 'homoresistant') %>%
        unnest(mix_phylotype) %>%
        nest(other_res = -c(sample_id, mix_phylotype)),
      by = c('sample_id', 'mix_phylotype')
    ) %>%
    mutate(mix_drug_res = map2(homo_res, other_res, bind_rows) %>%
             map(arrange_all)) %>%
    select(-homo_res, -other_res) %>%
    nest_mixtures() %>%
    mutate(drugs_resistant = map(mix_drug_res, function(x) {
      bind_rows(x, tibble(drugs = character())) %>%
        select(drugs) %>%
        separate_rows(drugs, sep = ';') %>%
        distinct() %>%
        pull(drugs) %>%
        sort()
    }))

}

