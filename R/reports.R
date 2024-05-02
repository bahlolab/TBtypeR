
#' @export
batch_report <- function(tbt_res, snp_dist,
                         sample_meta = NULL,
                         report_title = 'TBtypeR Batch Report',
                         report_filename = NULL,
                         browse = TRUE) {

  assert_that(
    is.data.frame(tbt_res),
    inherits(snp_dist, 'dist'),
    (is.null(sample_meta) || is.data.frame(sample_meta)),
    is_scalar_character(report_title),
    is.null(report_filename) || is_scalar_character(report_filename)
  )

  if ('mix_drug_res' %in% colnames(tbt_res)) {
    rmd_file <- system.file(file.path("Rmd", "TBtypeR_batch_report_DR.Rmd"),
                            package = "TBtypeR", mustWork = TRUE)
  } else {
    rmd_file <- system.file(file.path("Rmd", "TBtypeR_batch_report.Rmd"),
                            package = "TBtypeR", mustWork = TRUE)
  }

  tmp_dir <- tempdir(check = T)

  file.copy(rmd_file, tmp_dir)

  rmd_file <- file.path(tmp_dir, basename(rmd_file))

  rmd_env <- list2env(
    list(tbt_res = tbt_res,
         snp_dist = snp_dist,
         sample_meta = sample_meta,
         report_title = report_title),
    envir = new.env())

  if (is.null(report_filename)) {
    report_filename <-
      stringr::str_replace_all(report_title, '\\s+', '_') %>%
      stringr::str_c('.html')
  }

  report_dir <- dirname(report_filename)
  report_filename <- basename(report_filename)

  suppressMessages(suppressWarnings(
    rmarkdown::render(output_file = report_filename,
                      input = rmd_file,
                      envir = rmd_env
    )))

  file.copy(file.path(tmp_dir, report_filename),
            report_dir,
            overwrite = T)

  unlink(tmp_dir, recursive = T)

  report_filename <- normalizePath(file.path(report_dir, report_filename))

  if (browse) {
    browseURL(path.expand(report_filename))
  }

  return(report_filename)
}

subchunkify <- function(g, name, fig_height=7, fig_width=8) {
  g_deparsed <- paste0(deparse(function() {g}), collapse = '')
  sub_chunk <- paste0("`", "``{r sub_chunk_", name, ", fig.height=",
                      fig_height, ", fig.width=", fig_width, ", dpi=200, echo=FALSE}",
                      "\n(",
                      g_deparsed
                      , ")()",
                      "\n`","``")
  knitr::knit(text = knitr::knit_expand(text = sub_chunk), output = 'subchunkify.txt', quiet = TRUE)
  scan('subchunkify.txt', what = character(), sep = '\n') %>%
    stringr::str_extract('((?<=\\().+(?=\\)))|((?<=<img src=")[^"]+(?="))')
}

data_table <- function(x) {
  DT::datatable(x,
                rownames = F,
                filter = 'top',
                extensions = 'Buttons',
                options = list(dom = 'rtipB',
                               buttons = c('csv', 'excel')))
}

