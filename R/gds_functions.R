
#' @export
#' @importFrom SeqArray seqGetData seqSetFilter seqNumAllele
#' @importFrom SeqVarTools variantInfo
#' @importFrom dplyr distinct group_by ungroup as_tibble mutate filter arrange rename add_count
#' @importFrom dplyr inner_join anti_join semi_join
#' @importFrom tidyr pivot_wider unchop chop separate_rows
#' @importFrom magrittr "%>%"
get_allele_counts <- function(gds, panel,
                              as_tibble = FALSE,
                              verbose = FALSE,
                              max_ext_freq = 0.25) {

  # check_args
  assert_that(is_gds(gds),
              is.data.frame(panel),
              all(c('chrom', 'pos', 'ref', 'alt') %in% colnames(panel)),
              is_bool(verbose),
              is_scalar_double(max_ext_freq) && max_ext_freq >= 0 && max_ext_freq <= 1)

  # find matching sites in gds
  if (!is_open_gds(gds)) {
    warning('Opening closed gds file')
    gds <- seqOpen(gds$filename, allow.duplicate = T)
    on.exit({seqClose(gds)})
  }
  sam_id <- seqGetData(gds, 'sample.id')
  var_id <- seqGetData(gds, 'variant.id')
  gr <- with(panel, GenomicRanges::GRanges(chrom, IRanges::IRanges(start = pos, width = nchar(ref))))
  seqSetFilter(gds, gr, verbose = verbose)
  panel_var <- panel_with_vid(panel) %>% select(vid, chrom, pos, ref, alt)

  gds_panel_var <-
    variantInfo(gds) %>%
    as_tibble() %>%
    mutate(num_allele = seqNumAllele(gds)) %>%
    separate_rows(alt, sep = ',') %>%
    mutate(alt = if_else(nchar(alt)>0, alt, NA_character_)) %>%
    group_by(variant.id) %>%
    mutate(alt_index = if_else(is.na(alt), NA_integer_, cumsum(!is.na(alt)))) %>%
    ungroup() %>%
    rename(chrom = chr) %>%
    semi_join(panel_var, by = c('chrom', 'pos', 'ref')) %>%
    (function(x) {
      bind_rows(
        # ref allele
        inner_join(x, select(panel_var, -alt), by = c("chrom", "pos", "ref")) %>%
          left_join(mutate(panel_var, priority = TRUE), by = c("chrom", "pos", "ref", "alt", "vid")) %>%
          mutate(allele = 'ref', allele_index = 1L) %>%
          select(variant.id, vid, chrom, pos, allele, allele_index, num_allele, priority) %>%
          arrange(vid, priority) %>%
          group_by(vid) %>% slice(1) %>% ungroup() %>%
          select(-priority), # take first ref
        # alt allele
        inner_join(x, panel_var, by = c("chrom", "pos", "ref", "alt")) %>%
          mutate(allele = 'alt', allele_index = 1L + alt_index) %>%
          select(variant.id, vid, chrom, pos, allele, allele_index, num_allele) %>%
          group_by(vid) %>% slice(1) %>% ungroup(), # take first alt
        # ext allele
        anti_join(x, panel_var, by = c("chrom", "pos", "ref", "alt")) %>%
          inner_join(select(panel_var, -alt), by = c("chrom", "pos", "ref")) %>%
          group_by(vid, chrom, pos, ref, alt) %>% slice(1) %>% ungroup() %>%
          mutate(allele = 'ext', allele_index = 1L + alt_index) %>%
          select(variant.id, vid, chrom, pos, allele, allele_index, num_allele)
      )
    }) %>%
    mutate(allele = factor(allele, c('ref', 'alt', 'ext'))) %>%
    arrange_all() %>%
    chop(c(-variant.id, -num_allele)) %>%
    mutate(offset = cumsum(num_allele) - num_allele)

  seqSetFilter(gds, variant.id = gds_panel_var$variant.id, verbose = verbose)

  var_index <-
    gds_panel_var %>%
    select(-variant.id, -num_allele) %>%
    unchop(-offset) %>%
    mutate(allele_index = offset + allele_index) %>%
    select(vid, allele, allele_index) %>%
    (function(x) {
      filter(x, allele != 'ext') %>%
        pivot_wider(names_from = allele, values_from = allele_index, values_fill = NA) %>%
        left_join(
          filter(x, allele == 'ext') %>%
            rename(ext = allele_index) %>%
            select(vid, ext) %>%
            chop(ext),
          by = 'vid')
    })

  ext_index <-
    var_index %>%
    mutate(vidx = seq_along(vid),
           ii = map(ext, seq_along)) %>%
    select(vidx, ext, ii) %>%
    unnest(c(ext, ii)) %>%
    pivot_wider(names_from = ii, values_from = ext) %>%
    { left_join(tibble(vidx = seq_along(var_index$vid)), ., 'vidx') }


  AD <- seqGetData(gds, 'annotation/format/AD')$data
  # reset gds filter
  seqSetFilter(gds, variant.id = var_id, verbose=verbose)
  ref_ac <- AD[, var_index$ref, drop = F]
  alt_ac <-  AD[, var_index$alt, drop = F]
  alt_ac[is.na(alt_ac)] <- 0L
  alt_ac[is.na(ref_ac)] <- NA_integer_

  ext_ac <- matrix(0L, nrow(ref_ac), ncol(ref_ac))
  for (i in as.character(seq_len(ncol(ext_index)-1))) {
    ei <- na.omit(select(ext_index, vidx, adidx = all_of(i)))
    tmp <-  AD[, ei$adidx, drop=F]
    tmp[is.na(tmp)] <- 0L
    ext_ac[, ei$vidx] <- ext_ac[, ei$vidx, drop = F] + tmp
  }

  flt <- which(ext_ac / (ref_ac + alt_ac + ext_ac) > max_ext_freq)
  ref_ac[flt] <- NA_integer_
  alt_ac[flt] <- NA_integer_

  if (as_tibble) {
    allele_counts <-
      as_tibble(data.frame(sample_id = sam_id,
                           vid = rep(var_index$vid, each = length(sam_id)),
                           ref_ac = c(ref_ac),
                           alt_ac = c(alt_ac),
                           depth = c(ref_ac + alt_ac + ext_ac))) %>%
      arrange_all()

  } else {
    allele_counts <-
      array(c(ref_ac, alt_ac),
            dim = c(length(sam_id), nrow(var_index), 2),
            dimnames = list(sample = sam_id,
                            variant = var_index$vid,
                            allele = c('ref', 'alt')))
  }

  return(allele_counts)
}

is_gds <- function(x) {
  inherits(x, "SeqVarGDSClass")
}

is_open_gds <- function(x) {
  is_gds(x) && !xptr::is_null_xptr(x$ptr)
}

#' @importFrom magrittr "%>%"
#' @importFrom SeqArray seqGetData seqSetFilter seqOpen seqClose
# gds_get_AD_parallel <- function(gds, verbose = FALSE) {
#
#   fn = gds$filename
#   var.id = seqGetData(gds, 'variant.id')
#   sam.id = seqGetData(gds, 'sample.id')
#   workers = future::nbrOfWorkers()
#
#   parallel::splitIndices(length(var.id), workers) %>%
#     map( ~ var.id [.] ) %>%
#     { .[lengths(.) > 0 ] } %>%
#     furrr::future_map(function(variant.id) {
#       gds <- seqOpen(fn, allow.duplicate = T)
#       seqSetFilter(gds, variant.id = variant.id, sample.id = sam.id, verbose = verbose)
#       AD <- seqGetData(gds, 'annotation/format/AD')
#       seqClose(gds)
#       return(AD)
#     }) %>%
#     purrr::reduce(function(x, y) {
#       list(length = c(x$length, y$length),
#            data = cbind(x$data, y$data))
#     })
# }



