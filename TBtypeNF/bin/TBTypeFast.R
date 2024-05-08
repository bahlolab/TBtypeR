#!/usr/bin/env Rscript

stopifnot(
  require(SeqArray),
  require(TBtypeR),
  require(future),
  require(docopt),
  require(tidyverse)
)
"
Usage:
  TBTypeFast.R <fastlin_output> <panel> [options]

Options:
  fastlin_output         fastlin output tsv file.
  panel                  Phylotype panel
  --sample-id=<s>        Sample Identifier [default: sample].
  --threads=<i>          Number of parallel threads to uses [default: 8].
  --max-mix=<i>          Maximum number of mixture components [default: 2].
  --min-mix-prop=<d>     Minimum mixture proportion detectable [default: 0.005]
  --args-json=<a>        Additional arguments in json format for TBtypeR::tbtype
" -> doc

opts <- docopt(doc)

stopifnot(
  file.exists(opts$fastlin_output),
  file.exists(opts$panel)
)

threads <- max(as.integer(opts$threads), 1)
max_mix <- as.integer(opts$max_mix)
min_mix_prop <- as.numeric(opts$min_mix_prop)

options(future.globals.maxSize = Inf)
plan(multicore, workers = threads)

panel <- TBtypeR::read_panel(opts$panel)
allele_counts <- TBtypeR:::fastlin_allele_counts(
  opts$fastlin_output,
  panel = panel,
  sample_id = opts$sample_id
)

args <- list(
  allele_counts = allele_counts,
  panel = panel,
  max_phylotypes = max_mix,
  min_mix_prop = min_mix_prop
)

if (!is.null(opts$args_json)) {
  args <- c(args, jsonlite::fromJSON(opts$args_json))
}

tbt_res <-
  do.call(TBtypeR::tbtype, args) %>%
  TBtypeR::filter_tbtype(
    max_phylotypes = max_mix,
    min_mix_prop = min_mix_prop
  ) %>%
  TBtypeR::unnest_mixtures()

saveRDS(tbt_res, str_c(opts$sample_id, ".tbt_res.rds"))

tbt_res %>%
  bind_rows(
    tibble( 
      sample_id = character(),
      n_phy = numeric(),
      likelihood = numeric(),
      error_rate = numeric(),
      p_val_perm = numeric(),
      p_val_wsrst = numeric(),
      abs_diff = numeric(),
      mix_phylotype = character(),
      mix_prop = numeric())
  ) %>% 
  select(
    "sample_id", "n_phy", "likelihood", "error_rate", "p_val_perm",
    "p_val_wsrst", "abs_diff", "mix_phylotype", "mix_prop", 
  ) %>%
  write_csv(str_c(opts$sample_id, ".tbt_calls.csv"))
