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
  TBType.R <gds_file> <panel> [options]

Options:
  gds_file               gds filename.
  panel                  Phylotype panel.
  --sample-meta=<f>      TSV file with columns sample, group and label
  --output=<o>           Prefix of output files [default: output].
  --dr-panel=<p>         Drug resistance panel.
  --threads=<i>          Number of parallel threads to uses [default: 8].
  --max-mix=<i>          Maximum number of mixture components [default: 2].
  --args-json=<a>        Additional arguments in json format for TBtypeR::tbtype
  --no-report            Disable report output
" -> doc

opts <- docopt(doc)

stopifnot(
  file.exists(opts$gds_file),
  file.exists(opts$panel)
)

threads <- max(as.integer(opts$threads), 1)
max_mix <- as.integer(opts$max_mix)

options(future.globals.maxSize = Inf)
plan(multicore, workers = threads)

gds <- seqOpen(opts$gds)
panel <- TBtypeR::read_panel(opts$panel)

sample_meta <- NULL
if (!is.null(opts$sample_meta)) {
  tmp <- read_tsv(opts$sample_meta)
  if ("sample" %in% colnames(tmp) && sum(c("group", "label") %in% colnames(tmp)) > 0) {
    sample_meta <- tmp
  }
}

args <- list(gds = gds, panel = panel)
if (!is.null(opts$args_json)) {
  args <- c(args, jsonlite::fromJSON(opts$args_json))
}

tbt_res <-
  do.call(TBtypeR::tbtype, args) %>%
  TBtypeR::filter_tbtype(max_phylotypes = max_mix) %>%
  TBtypeR::unnest_mixtures()

saveRDS(tbt_res, str_c(opts$output, ".tbt_res.rds"))

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
  write_csv(str_c(opts$output, ".tbt_calls.csv"))

snp_dist <- TBtypeR::snp_distance(gds, regions = NULL)

if (!is.null(opts$dr_panel)) {
  dr_panel <- TBtypeR::read_panel(opts$dr_panel, phylo = FALSE)
  tbt_res <-
    tbt_res %>%
    TBtypeR::assign_dr(
      gds = gds,
      dr_panel = dr_panel
    )
}

if (!opts$no_report) {
  TBtypeR::batch_report(
    tbt_res = tbt_res,
    snp_dist = snp_dist,
    report_title = opts$output,
    sample_meta = sample_meta,
    report_filename = str_c(opts$output, ".html"),
    browse = FALSE
  )
}

