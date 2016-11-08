
#'
#' chillingham_genotype
#'
#' A dataset containing the 0/1/2 SNP genotypes of 17 Chillingham cattle on BTA1
#'
#' \itemize{
#'   \item{\code{FID}}: breeds (race)
#'   \item{\code{IID}}: id of animals (individual)
#'   \item{\code{PAT}}: paternal (who is father?)
#'   \item{\code{MAT}}: maternal
#'   \item{\code{SEX}}: sex (1/2)?
#'   \item{\code{PHENOTYPE}}: phenotype (or -9 if missing)
#'   \item{\code{BovineHD0100015049_0}}: sample name
#'   \item{\code{BovineHD0100015050_0}}: sample name
#'   \item{\code{BovineHD0100015051_0}}: sample name
#'   \item{\code{...}}
#' }
#' @examples
#' data(chillingam)
#' head(chillingham_genotype)

#'
#' @source Chillingham dataset (cite)
#' @format A data frame with 17 rows and 2641 variables
#' @keywords datasets
#'
"chillingham_genotype"

#'
#' chillingham_map
#'
#' A dataset containing the Map file for SNP on BTA1
#'
#' \describe{
#'   \item{\code{Chrom}}{chromosome name}
#'   \item{\code{SNP}}{SNP name}
#'   \item{\code{cM}}{position in cM}
#'   \item{\code{bps}}{position in bp}
#' }
#' @examples
#' data(chillingam)
#' head(chillingham_map)
#'
#' @source Chillingham dataset (cite)
#' @format A data frame with 2635 rows and 4 variables
#' @keywords datasets
#'
"chillingham_map"
