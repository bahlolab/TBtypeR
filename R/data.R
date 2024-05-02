#' TBtypeR SNP phylotype panel.
#'
#' @format A data frame with 10,904 rows and 8 columns:
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
#' @source Zwyer et al. 2021: \url{https://doi.org/10.12688/openreseurope.14029.2}
#' @source Shuaib et al. 2022: \url{https://doi.org/10.3390/genes13060990}
"tbt_panel"


#' Panel of TB drug resistance mutations derived from the WHO Mutation Catalogue 2021
#'
#' @format A data frame with 1,145 rows and 7 columns:
#' \describe{
#'   \item{chrom}{H37Rv chromosome}
#'   \item{pos}{H37Rv genomic position}
#'   \item{ref}{H37Rv reference base}
#'   \item{alt}{H37Rv alternate base}
#'   \item{drugs}{Three letter code for antibiotics to which polymorphism is associated with resistance, separated by ';'}
#'   \item{ppv}{Positive predictive value for each drug, separated by ';'}
#'   \item{referece}{reference}
#' }
#' @source Catalogue of mutations in Mycobacterium tuberculosis complex and their association with drug resistance,
#' WHO TEAM Global Tuberculosis Programme, 2021 \url{https://www.who.int/publications/i/item/9789240028173}
#'
#' Licence: Creative Commons Attribution-NonCommercial-ShareAlike 3.0 IGO (CC BY-NC-SA 3.0 IGO) \url{https://creativecommons.org/licenses/by-nc-sa/3.0/igo/deed.en}
"who_mut_cat_2021"

#' MTBC core genomic regions (H37RV coordinates)
#'
#' @format GRanges object with 1155 ranges and 0 metadata columns
#' @source Pan and Core Genome Analysis of 183 Mycobacterium tuberculosis Strains Revealed a High Inter-Species Diversity among the Human Adapted Strains \url{https://doi.org/10.3390/antibiotics10050500}
#'
"h37rv_core_genome"
