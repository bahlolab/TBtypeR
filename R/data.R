#' TBtyper SNP phylotype panel.
#'
#' @format A data frame with 9,272 rows and 8 columns:
#' \describe{
#'   \item{chrom}{H37Rv chromosome}
#'   \item{pos}{H37Rv genomic position}
#'   \item{ref}{H37Rv reference base}
#'   \item{alt}{H37Rv alternate base}
#'   \item{genotype}{The genoptype associated with the phylotype, 0 for referene, 1 for alternate}
#'   \item{phylotype}{Name of MTBC phylotype}
#'   \item{parent_phylotype}{Name of parent MTBC phylotype}
#'   \item{referece}{Paper referencing the SNP marker}
#' }
#' @source Napier et al. 2020: \url{https://doi.org/10.1186/s13073-020-00817-3}
#' @source Thawornwattana et al. 2021: \url{https://doi.org/10.1099/mgen.0.000697}
#' @source Coscolla et al. 2021: \url{https://doi.org/10.1099/mgen.0.000477}
"tbt_panel"
