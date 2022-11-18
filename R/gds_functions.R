
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
  assert_that(
    is_gds(gds),
    is.data.frame(panel),
    all(c("chrom", "pos", "ref", "alt") %in% colnames(panel)),
    is_bool(verbose),
    is_scalar_double(max_ext_freq) && max_ext_freq >= 0 && max_ext_freq <= 1
  )
  # ensure gds is open
  if (!is_open_gds(gds)) {
    warning("Opening closed gds file")
    gds <- seqOpen(gds$filename, allow.duplicate = T)
    on.exit({
      seqClose(gds)
    })
  }
  # find matching sites in gds
  sam_id <- seqGetData(gds, "sample.id")
  var_id <- seqGetData(gds, "variant.id")
  gr <- with(panel, GenomicRanges::GRanges(chrom, IRanges::IRanges(start = pos, width = nchar(ref))))
  seqSetFilter(gds, gr, verbose = verbose)
  panel_var <- panel_with_vid(panel) %>% select(vid, chrom, pos, ref, alt)

  gds_panel_var <-
    variantInfo(gds) %>%
    as_tibble() %>%
    mutate(num_allele = seqNumAllele(gds)) %>%
    separate_rows(alt, sep = ",") %>%
    mutate(alt = if_else(nchar(alt) > 0, alt, NA_character_)) %>%
    group_by(variant.id) %>%
    mutate(alt_index = if_else(is.na(alt), NA_integer_, cumsum(!is.na(alt)))) %>%
    ungroup() %>%
    rename(chrom = chr) %>%
    semi_join(panel_var, by = c("chrom", "pos", "ref")) %>%
    (function(x) {
      bind_rows(
        # ref allele
        inner_join(x, select(panel_var, -alt), by = c("chrom", "pos", "ref")) %>%
          left_join(mutate(panel_var, priority = TRUE), by = c("chrom", "pos", "ref", "alt", "vid")) %>%
          mutate(allele = "ref", allele_index = 1L) %>%
          select(variant.id, vid, chrom, pos, allele, allele_index, num_allele, priority) %>%
          arrange(vid, priority) %>%
          group_by(vid) %>%
          slice(1) %>%
          ungroup() %>%
          select(-priority), # take first ref
        # alt allele
        inner_join(x, panel_var, by = c("chrom", "pos", "ref", "alt")) %>%
          mutate(allele = "alt", allele_index = 1L + alt_index) %>%
          select(variant.id, vid, chrom, pos, allele, allele_index, num_allele) %>%
          group_by(vid) %>%
          slice(1) %>%
          ungroup(), # take first alt
        # ext allele
        anti_join(x, panel_var, by = c("chrom", "pos", "ref", "alt")) %>%
          inner_join(select(panel_var, -alt), by = c("chrom", "pos", "ref")) %>%
          group_by(vid, chrom, pos, ref, alt) %>%
          slice(1) %>%
          ungroup() %>%
          mutate(allele = "ext", allele_index = 1L + alt_index) %>%
          select(variant.id, vid, chrom, pos, allele, allele_index, num_allele)
      )
    }) %>%
    mutate(allele = factor(allele, c("ref", "alt", "ext"))) %>%
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
      filter(x, allele != "ext") %>%
        pivot_wider(names_from = allele, values_from = allele_index, values_fill = NA) %>%
        left_join(
          filter(x, allele == "ext") %>%
            rename(ext = allele_index) %>%
            select(vid, ext) %>%
            chop(ext),
          by = "vid"
        )
    })

  ext_index <-
    var_index %>%
    mutate(
      vidx = seq_along(vid),
      ii = map(ext, seq_along)
    ) %>%
    select(vidx, ext, ii) %>%
    unnest(c(ext, ii)) %>%
    pivot_wider(names_from = ii, values_from = ext) %>%
    {
      left_join(tibble(vidx = seq_along(var_index$vid)), ., "vidx")
    }


  AD <- seqGetData(gds, "annotation/format/AD")$data
  # reset gds filter
  seqSetFilter(gds, variant.id = var_id, verbose = verbose)
  ref_ac <- AD[, var_index$ref, drop = F]
  alt_ac <- AD[, var_index$alt, drop = F]
  alt_ac[is.na(alt_ac)] <- 0L
  alt_ac[is.na(ref_ac)] <- NA_integer_

  ext_ac <- matrix(0L, nrow(ref_ac), ncol(ref_ac))
  for (i in as.character(seq_len(ncol(ext_index) - 1))) {
    ei <- na.omit(select(ext_index, vidx, adidx = all_of(i)))
    tmp <- AD[, ei$adidx, drop = F]
    tmp[is.na(tmp)] <- 0L
    ext_ac[, ei$vidx] <- ext_ac[, ei$vidx, drop = F] + tmp
  }

  flt <- which(ext_ac / (ref_ac + alt_ac + ext_ac) > max_ext_freq)
  ref_ac[flt] <- NA_integer_
  alt_ac[flt] <- NA_integer_

  if (as_tibble) {
    allele_counts <-
      as_tibble(data.frame(
        sample_id = sam_id,
        vid = rep(var_index$vid, each = length(sam_id)),
        ref_ac = c(ref_ac),
        alt_ac = c(alt_ac),
        depth = c(ref_ac + alt_ac + ext_ac)
      )) %>%
      mutate(baf = alt_ac / depth) %>%
      arrange_all()
  } else {
    allele_counts <-
      array(c(ref_ac, alt_ac),
        dim = c(length(sam_id), nrow(var_index), 2),
        dimnames = list(
          sample = sam_id,
          variant = var_index$vid,
          allele = c("ref", "alt")
        )
      )
  }

  return(allele_counts)
}

is_gds <- function(x) {
  inherits(x, "SeqVarGDSClass")
}

is_open_gds <- function(x) {
  is_gds(x) && !xptr::is_null_xptr(x$ptr)
}

#' @export
#' @importFrom SeqArray seqGetData seqSetFilter seqClose seqOpen
#' @importFrom magrittr "%>%"
snp_distance <- function(gds,
                         regions = get_core_regions(),
                         allelic_mode = 4L,
                         min_median_depth = 5L,
                         max_var_missing = 0.05,
                         max_sample_missing = 0.05,
                         threads = min(8L, parallel::detectCores())) {

  assert_that(
    is_gds(gds),
    is.null(regions) || inherits(regions, 'GenomicRanges'),
    is_scalar_integer(allelic_mode) && allelic_mode > 1L,
    is_scalar_integer(min_median_depth) && min_median_depth >= 0L,
    is_scalar_proportion(max_sample_missing),
    is_scalar_proportion(max_var_missing)
  )
  # ensure gds is open
  if (!is_open_gds(gds)) {
    warning("Opening closed gds file")
    gds <- seqOpen(gds$filename, allow.duplicate = T)
    on.exit({
      seqClose(gds)
    })
  }
  # filter gds to regions and allelic_mode
  seqSetFilter(gds, regions)
  seqSetFilter(gds, SeqVarTools::isSNV(gds, biallelic=FALSE), action = 'intersect')
  seqSetFilter(gds, seqNumAllele(gds) == allelic_mode, action = 'intersect')
  assert_that(length(seqGetData(gds, 'variant.id')) > 0)

  # assign GT based on allelic depth
  AD0 <- seqGetData(gds, 'annotation/format/AD')$data
  nsam <- nrow(AD0)
  nvar <- ncol(AD0) / allelic_mode

  GT <- matrix(0L, nrow=nsam, ncol=nvar) %>%
    set_rownames(seqGetData(gds, 'sample.id'))
  AD_MX <- AD0[, 1 + 4L * (seq_len(nvar) - 1L)]
  DP <- AD_MX
  GT[is.na(AD_MX)] <- NA_integer_
  for (i in 2:allelic_mode) {
    AD_i <- AD0[, i + 4L * (seq_len(nvar) - 1L)]
    DP <- DP + AD_i
    is_greater <- which(AD_i > AD_MX)
    AD_MX[is_greater] <- AD_i[is_greater]
    GT[is_greater] <- i-1L
  }
  sam_med_dp <- apply(DP, 1, median, na.rm = TRUE)
  rm(DP, AD0, AD_i, AD_MX)

  # filter samples on median depth
  low_dp <- which(sam_med_dp < min_median_depth)
  if (length(low_dp)) {
    message("Excluded ", length(low_dp), '/', nsam, ' samples due to low depth')
    GT <- GT[-low_dp, ]
    nsam <- nrow(GT)
  }

  # filter variants or samples first, maximizing remaining data
  var_miss_0 <- colSums(is.na(GT)) / nsam
  sam_miss_0 <- rowSums(is.na(GT)) / nvar
  var_set_0 <- which(var_miss_0 <= max_var_missing)
  sam_set_0 <- which(sam_miss_0 <= max_sample_missing)

  var_miss_1 <- colSums(is.na(GT[sam_set_0, ])) / length(sam_set_0)
  sam_miss_1 <- rowSums(is.na(GT[, var_set_0])) / length(var_set_0)
  var_set_1 <- which(var_miss_1 <= max_var_missing)
  sam_set_1 <- which(sam_miss_1 <= max_sample_missing)

  n_v0_s1 <- length(var_set_0) * length(sam_set_1)
  n_v1_s0 <- length(var_set_1) * length(sam_set_0)

  if (n_v0_s1 >= n_v1_s0) {
    var_set <-  var_set_0
    sam_set <- sam_set_1
  } else {
    var_set <-  var_set_1
    sam_set <- sam_set_0
  }
  if (length(sam_set) < nsam) {
    message("Excluded ", nsam - length(sam_set), '/', nsam, ' samples due to high missingness')
  }
  if (length(var_set) < nvar) {
    message("Excluded ", nvar - length(var_set), '/', nvar, ' variants due to high missingness')
  }

  GT <- GT[sam_set, var_set]

  pdist(GT, threads = threads)

}
