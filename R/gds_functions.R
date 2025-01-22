
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
                              max_ext_freq = 0.25,
                              rename_chrom = character()) {


  # check_args
  assert_that(
    is_gds(gds),
    is.data.frame(panel),
    all(c("chrom", "pos", "ref", "alt") %in% colnames(panel)),
    is_bool(verbose),
    is_scalar_double(max_ext_freq) && max_ext_freq >= 0 && max_ext_freq <= 1,
    is.character(rename_chrom),
    length(rename_chrom) == 0 | !is.null(names(rename_chrom))
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
  gr <-
    panel %>%
    select(chrom, pos, ref) %>%
    mutate(chrom = map(chrom, function(x) { c(x, names(rename_chrom)[rename_chrom == x])})) %>%
    unnest(chrom) %>%
    arrange_all() %>%
    with(., GenomicRanges::GRanges(chrom, IRanges::IRanges(start = pos, width = nchar(ref))))

  seqSetFilter(gds, gr, verbose = verbose)
  if (length(seqGetData(gds, "variant.id")) == 0) {
    stop("No matching panel variants found")
  }

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
    mutate(chrom =
             replace(
               chrom,
               chrom %in% names(rename_chrom),
               rename_chrom[chrom[chrom %in% names(rename_chrom)]]
             )
    ) %>%
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
  is_gds(x) && !identical(x$ptr, new("externalptr"))
}

check_valid_gds <- function(gds, filename = NULL) {

  fn <- `if`(is.null(filename), gds$filename, filename)

  if (!is_open_gds(gds)) {
    warning("Opening closed gds file")
    gds <- seqOpen(gds$filename, allow.duplicate = T)
    on.exit({
      seqClose(gds)
    })
  }

  sample_id <- tryCatch(seqGetData(gds, 'sample.id'), error = function(e) NULL)
  if (is.null(sample_id) | length(sample_id) == 0) {
    stop("No samples present in VCF/GDS file \"", fn, "\"")
  }

  ad <- tryCatch(seqGetData(gds, "annotation/format/AD")$data, error = function(e) NULL)
  if (is.null(ad)) {
    stop(
      "Missing required VCF format field \"AD\" (Allelic Depth) from \"", fn, "\". ",
      "It is recommended to create VCF with BCFtools or TBtypeNF pipeline."
    )
  }
  if (!( is.matrix(ad) & is.integer(ad))) {
    stop(
      "VCF format field \"AD\" has incorrect format or data type. " ,
      "It is recommended to create VCF with BCFtools or TBtypeNF pipeline."
    )
  }
}


#' @export
#' @importFrom SeqArray seqGetData seqSetFilter seqClose seqOpen
#' @importFrom magrittr "%>%"
gds_filter_gt <- function(gds,
                          regions = TBtypeR::h37rv_core_genome,
                          allelic_mode = 4L,
                          MIN_MEDIAN_DEPTH = 5L,
                          MAX_VAR_MISSING = 0.05,
                          MAX_SAMPLE_MISSING = 0.05,
                          MIN_VAR_DP_SD = 3,
                          DROP_INVARIANT = TRUE,
                          source = c('GT', 'AD'),
                          biallelic_only = FALSE,
                          MIN_C_FOLD = 5) {

  source <- match.arg(source)

  assert_that(
    is_gds(gds),
    is.null(regions) || inherits(regions, 'GenomicRanges'),
    is_scalar_integer(allelic_mode) && allelic_mode > 1L,
    is_scalar_integer(MIN_MEDIAN_DEPTH) && MIN_MEDIAN_DEPTH >= 0L,
    is_scalar_proportion(MAX_SAMPLE_MISSING),
    is_scalar_proportion(MAX_VAR_MISSING),
    is_scalar_character(source)
  )

  # ensure gds is open
  if (!is_open_gds(gds)) {
    warning("Opening closed gds file")
    gds <- seqOpen(gds$filename, allow.duplicate = T)
    on.exit({
      seqClose(gds)
    })
  }
  # restore filter on exit
  flt <- SeqArray::seqGetFilter(gds)
  on.exit({SeqArray::seqSetFilter(gds, sample.sel = flt$sample.sel, variant.sel = flt$variant.sel)})

  # filter gds to regions and allelic_mode
  if (!is.null(regions)) { seqSetFilter(gds, regions) }
  seqSetFilter(gds, SeqVarTools::isSNV(gds, biallelic=FALSE), action = 'intersect')
  seqSetFilter(gds, seqNumAllele(gds) == allelic_mode, action = 'intersect')
  assert_that(length(seqGetData(gds, 'variant.id')) > 0)

  # assign GT based on allelic depth, create DP
  AD0 <- seqGetData(gds, 'annotation/format/AD')$data
  assert_that(ncol(AD0) %% allelic_mode == 0)

  nsam <- nrow(AD0)
  nvar <- ncol(AD0) / allelic_mode

  GT <- matrix(0L, nrow=nsam, ncol=nvar) %>%
    set_rownames(seqGetData(gds, 'sample.id')) %>%
    set_colnames(seqGetData(gds, 'variant.id'))

  AD_MX <- AD0[, 1 + allelic_mode * (seq_len(nvar) - 1L)]
  DP <- AD_MX
  GT[is.na(AD_MX)] <- NA_integer_
  for (i in 2:allelic_mode) {
    AD_i <- AD0[, i + allelic_mode * (seq_len(nvar) - 1L)]
    DP <- DP + AD_i
    is_greater <- which(AD_i > AD_MX)
    AD_MX[is_greater] <- AD_i[is_greater]
    GT[is_greater] <- i-1L
  }

  rm(AD0, AD_i, AD_MX)

  if (source == 'GT') {
    # replace GT with called genotypes

    gt12 <- seqGetData(gds, 'genotype')
    GT <- gt12[1,,]
    GT[GT != gt12[2,,]] <- NA_integer_
    rownames(GT) <- seqGetData(gds, 'sample.id')
    colnames(GT) <- seqGetData(gds, 'variant.id')
    rm(gt12)

  }

  if (biallelic_only) {
    AAC <- matrix(0L, nrow = allelic_mode, ncol = ncol(GT))
    for (i in seq_len(allelic_mode)) {
      AAC[i, ] <- as.integer(colSums(GT == i-1L, na.rm = T))
    }
    AACO <- apply(AAC, 2, order, decreasing = T)
    AACR <- vapply(seq_len(ncol(AACO)), function(i) AAC[,i][AACO[,i]], integer(nrow(AACO)))
    CAC <- colSums(AACR[-(1:2), ], na.rm = T)
    biallelic <- which(AACR[2, ] > MIN_C_FOLD * CAC)

    # redact non-bialleic sites with -1L (make NA after missingness filers)
    tibble(i = intersect(biallelic, which(CAC > 0))) %>%
      mutate(
        a1 = AACO[1, i] -1L,
        a2 = AACO[2, i] -1L) %>%
      chop(i) %>%
      mutate(a12 = map2(a1, a2, ~ sort(c(.x, .y)))) %>%
      select(i, a12) %>%
      chop(i) %>%
      mutate(i = map(i, ~ sort(unique(unlist(.))))) %>%
      pwalk(function(a12, i) {
        GT[, i][!GT[, i] %in% c(a12, NA_integer_)] <<- -1L
      })

    message("Excluded ", ncol(GT) - length(biallelic), '/', nvar, ' variants due not being biallelic')
    GT <- GT[, biallelic, drop = F]
    DP <- DP[, biallelic, drop = F]
  }

  if (DROP_INVARIANT) {
    drop <- which(colSums(GT, na.rm = T) == colMeans(GT, na.rm = T))
    if (length(drop)) {
      message("Excluded ", length(drop), '/', nvar, ' variants due being invariant')
      GT <- GT[, -drop, drop = F]
      DP <- DP[, -drop, drop = F]
    }
  }

  samples_exclude <- integer()
  # filter samples on median depth
  sam_low_dp <-
    (apply(DP, 1, median, na.rm = TRUE) < MIN_MEDIAN_DEPTH) %>%
    replace_na(FALSE)

  if (sum(sam_low_dp)) {
    message("Excluded ", sum(sam_low_dp), '/', nsam, ' samples due to low depth')
    samples_exclude <- which(sam_low_dp)
  }

  # filter samples by missingness
  sam_high_miss <- (rowSums(is.na(GT)) / nvar > MAX_SAMPLE_MISSING) & (!sam_low_dp)

  if (sum(sam_high_miss)) {
    message("Excluded ", sum(sam_high_miss), '/', nsam, ' samples due to high missingness')
    samples_exclude <- union(samples_exclude, which(sam_high_miss))
  }

  if (length(samples_exclude)) {
    GT <- GT[-samples_exclude,]
    DP <- DP[-samples_exclude,]
    nsam <- nrow(GT)
  }

  variants_exclude <- integer()

  # filter variants by depth
  var_med_dp <- apply(DP, 2, median, na.rm=TRUE)
  min_var_med_dp <- mean(var_med_dp) - MIN_VAR_DP_SD*sd(var_med_dp)
  var_low_dp <- (var_med_dp < min_var_med_dp) %>% replace_na(FALSE)

  if (sum(var_low_dp)) {
    message("Excluded ", sum(var_low_dp), '/', nvar, ' variants due to low depth')
    variants_exclude <- which(var_low_dp)
  }

  # filter variants by missingness
  var_high_miss <- (colSums(is.na(GT)) / nsam > MAX_VAR_MISSING) & (!var_low_dp)

  if (sum(var_high_miss)) {
    message("Excluded ", sum(var_high_miss), '/', nvar, ' variants due to high missingness')
    variants_exclude <- union(variants_exclude, which(var_high_miss))
  }

  if (length(variants_exclude)) {
    GT <- GT[, -variants_exclude]
    # DP <- DP[, -variants_exclude]
    nvar <- ncol(GT)
  }

  if (biallelic_only) {
    GT[GT == -1L] <- NA_integer_
  }

  GT
}

#' @export
snp_distance <- function(
    gds,
    regions = TBtypeR::h37rv_core_genome,
    allelic_mode = 4L,
    MIN_MEDIAN_DEPTH = 5L,
    MAX_VAR_MISSING = 0.05,
    MAX_SAMPLE_MISSING = 0.05,
    MIN_VAR_DP_SD = 3,
    source = c('GT', 'AD'),
    threads = min(8L, parallel::detectCores())
)
{

  source <- match.arg(source)

  GT <- gds_filter_gt(
    gds = gds,
    regions = regions,
    allelic_mode = 4L,
    MIN_MEDIAN_DEPTH = MIN_MEDIAN_DEPTH,
    MAX_VAR_MISSING = MAX_VAR_MISSING,
    MAX_SAMPLE_MISSING = MAX_SAMPLE_MISSING,
    MIN_VAR_DP_SD = MIN_VAR_DP_SD,
    source = source
  )

  pdist(GT, threads = threads)
}

#' @importFrom SeqArray seqOpen seqVCF2GDS
#' @importFrom assertthat assert_that
#' @importFrom rlang is_scalar_character
#' @importFrom stringr str_replace
vcf_to_gds <- function(vcf_fn) {

  assert_that(is_scalar_character(vcf_fn),
              file.exists(vcf_fn))

  gds_fn <- file.path(tempdir(), str_c(basename(vcf_fn), '.gds'))
  if (! file.exists(gds_fn)) {
    message("create GDS file: ", gds_fn)
    seqVCF2GDS(vcf.fn = vcf_fn,
               out.fn = gds_fn,
               storage.option = 'ZIP_RA',
               ignore.chr.prefix = '',
               verbose = FALSE)
  }
  seqOpen(gds_fn, allow.duplicate = TRUE)
}
