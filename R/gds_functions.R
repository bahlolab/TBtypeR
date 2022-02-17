
#' @export
#' @importFrom SeqArray seqGetData seqSetFilter seqNumAllele
#' @importFrom SeqVarTools variantInfo
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter arrange rename
#' @importFrom magrittr "%>%"
#' @importFrom tidyr spread
get_allele_counts_gds <- function(gds,
                                  var_info = NULL,
                                  verbose = FALSE,
                                  max_ext_freq = 0.25) {

  if (is.null(var_info)) {
    var_info <- varinf_v1
  }

  # check_args
  stopifnot(is_gds(gds),
            is.data.frame(var_info),
            setequal(c('variant_id', 'chr', 'pos', 'ref', 'alt'), colnames(var_info)),
            is_scalar_double(max_ext_freq) && max_ext_freq >= 0 && max_ext_freq <= 1)

  # find matching sites in gds
  var_id <- seqGetData(gds, 'variant.id')
  sam_id <- seqGetData(gds, 'sample.id')
  gr <- with(var_info, GenomicRanges::GRanges(chr, IRanges::IRanges(start = pos, width = 1L)))
  seqSetFilter(gds, gr, verbose = verbose)

  var_index <-
    variantInfo(gds, expand = TRUE) %>%
    as_tibble() %>%
    left_join(tibble(variant.id = seqGetData(gds, 'variant.id'),
                     num_allele = seqNumAllele(gds)),
              'variant.id') %>%
    inner_join(rename(var_info, target = alt), by = c('chr', 'pos', 'ref')) %>%
    (function(x) {
      bind_rows(
        filter(x, allele.index == 1) %>% mutate(allele.index = 0) %>% rename(allele = ref) %>% select(-alt),
        rename(x, allele = alt) %>% select(-ref)
      )
    }) %>%
    arrange(variant.id, allele.index) %>%
    filter(allele.index == 0 | allele == target) %>%
    select(-target) %>%
    group_by(variant_id) %>%
    filter(n() > 1) %>%
    ungroup() %>%
      (function(x) {
        inner_join(
          filter(x, allele.index == 0) %>%
            select(variant.id, num_allele, ref = allele) %>%
            distinct(),
          filter(x, allele.index > 0) %>%
            select(variant.id, alt = allele, alt.index = allele.index, variant_id) %>%
            chop(-variant.id),
          by = 'variant.id')
      }) %>%
    arrange(variant.id) %>%
    mutate(ref.index = 1L + cumsum(num_allele) - num_allele,
           alt.index = map2(ref.index, alt.index, ~ .x + .y)) %>%
    mutate(., ext.index = purrr::pmap(., function(num_allele, ref.index, alt.index, ...) {
      `if`(num_allele <= 2,
           integer(0),
           setdiff(seq.int(from = ref.index + 1, to = ref.index + num_allele -1), alt.index))
    })) %>%
    unnest(c(alt, alt.index, variant_id))

  # extract Allele Counts
  seqSetFilter(gds, variant.id = unique(var_index$variant.id), sample.id = sam_id, verbose = verbose)
  AD <- gds_get_AD_parallel(gds, verbose = verbose)$data

  ref_ac <- AD[, var_index$ref.index, drop = FALSE]
  alt_ac <- AD[, var_index$alt.index, drop = FALSE]

  ext_ac <- matrix(NA_integer_, nrow = nrow(ref_ac),ncol = ncol(ref_ac))
  for (i in seq_len(ncol(ref_ac))) {
    ext_ac[, i] <- rowSums(AD[, var_index$ext.index[[i]], drop = F], na.rm = T)
  }

  alt_ac[is.na(alt_ac)] <- 0L
  alt_ac[is.na(ref_ac)] <- NA_integer_
  flt <- which( ext_ac / (alt_ac + ref_ac + ext_ac) > max_ext_freq)
  ref_ac[flt] <- NA_integer_
  alt_ac[flt] <- NA_integer_

  allele_counts <-
    array(c(ref_ac, alt_ac),
          dim = c(length(sam_id), nrow(var_index), 2),
          dimnames = list(sample = sam_id,
                          variant = var_index$variant_id,
                          allele = c('Ref', 'Alt')))

  return(allele_counts)
}


#' @importFrom magrittr "%>%"
#' @importFrom SeqArray seqGetData seqSetFilter
gds_get_AD_parallel <- function(gds, verbose = FALSE) {

  fn = gds$filename
  var.id = seqGetData(gds, 'variant.id')
  sam.id = seqGetData(gds, 'sample.id')
  workers = future::nbrOfWorkers()

  parallel::splitIndices(length(var.id), workers) %>%
    map( ~ var.id [.] ) %>%
    { .[lengths(.) > 0 ] } %>%
    furrr::future_map( ~{
      gds <- SeqArray::seqOpen(fn, allow.duplicate = T)
      seqSetFilter(gds, variant.id = ., sample.id = sam.id, verbose = verbose)
      seqGetData(gds, 'annotation/format/AD')
    }) %>%
    purrr::reduce(function(x, y) {
      list(length = c(x$length, y$length),
           data = cbind(x$data, y$data))
    })
}








