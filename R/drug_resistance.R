
#' @export
#' @importFrom dplyr filter select mutate case_when if_else
#' @importFrom tidyr nest chop unnest expand_grid
#' @importFrom purrr map_lgl
assign_dr <- function(tbtype_results, gds,
                      dr_panel = TBtyper::who_dr_panel,
                      verbose = FALSE,
                      min_allele_count = 3,
                      min_allele_freq = 0.01,
                      max_p_val = 0.01) {

  # TODO: check tbtype results
  assert_that(
    is.data.frame(tbtype_results),
    is_gds(gds),
    is.null(dr_panel) || check_panel(dr_panel, "dr"),
    is_scalar_integerish(min_allele_count) && min_allele_count > 0,
    is_scalar_proportion(min_allele_freq)
  )

  if (!is_open_gds(gds)) {
    warning("Opening closed gds file")
    gds <- seqOpen(gds$filename, allow.duplicate = T)
    on.exit({
      seqClose(gds)
    })
  }

  SeqArray::seqSetFilter(gds, sample.id = unique(tbtype_results$sample_id), verbose = verbose)

  dr_counts <-
    get_allele_counts(gds,
      panel = dr_panel,
      as_tibble = TRUE,
      verbose = verbose,
      max_ext_freq = 1
    ) %>%
    filter(
      alt_ac >= min_allele_count,
      baf >= min_allele_freq
    ) %>%
    inner_join(panel_with_vid(dr_panel) %>% select(vid, drugs),
      by = "vid"
    ) %>%
    select(sample_id, vid, alt_ac, depth, baf, drugs) %>%
    nest(drug_ac = c(-sample_id)) %>%
    inner_join(tbtype_results %>%
      filter(!is.na(error_rate)) %>%
      select(sample_id, error_rate, phylotype, mix_prop) %>%
      unnest(c(mix_prop, phylotype)) %>%
      nest(type_data = c(-sample_id)),
    by = "sample_id"
    ) %>%
    mutate(drug_assign = map2(drug_ac, type_data, function(drug_ac, type_data) {
      err <- first(type_data$error_rate)
      mix_n <- nrow(type_data)
      tibble(
        phylotype = combinations(type_data$phylotype),
        phy_set = seq_along(phylotype)
      ) %>%
        unnest(phylotype) %>%
        left_join(select(type_data, phylotype, mix_prop),
          by = "phylotype"
        ) %>%
        chop(-phy_set) %>%
        mutate(
          mix_prop = map_dbl(mix_prop, sum, na.rm = T),
          p = mix_prop * (1 - err) + (1 - mix_prop) * (err)
        ) %>%
        expand_grid(drug_ac, .) %>%
        rowwise() %>%
        mutate(p.value = binom.test(alt_ac, depth, p)$p.value) %>%
        group_by(vid) %>%
        slice(which.max(p.value)) %>%
        ungroup() %>%
        mutate(
          status = case_when(
            p.value < max_p_val ~ "uncertain",
            map_lgl(phylotype, ~ length(.) == 1 && is.na(.)) ~ "none",
            map_lgl(phylotype, ~ length(.) == mix_n) ~ "all",
            map_lgl(phylotype, ~ length(.) == 1) ~ "single",
            map_lgl(phylotype, ~ length(.) > 1) ~ "multiple"
          ),
          phylotype = as.list(phylotype),
          p = if_else(status == "uncertain", NA_real_, p),
          p.value = if_else(status == "uncertain", NA_real_, p.value)
        ) %>%
        select(vid, baf, alt_ac, depth, drugs, status, phylotype, p.exp = p, p.value)
    })) %>%
    select(sample_id, drug_assign, type_data) %>%
    unnest(drug_assign)
  # separate_rows(drugs, sep = ';')

  # Need reasonable prior-probabilities on strain-snp relationship
  # How to deal with hetero-resistance?
}
