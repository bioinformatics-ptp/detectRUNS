#' Chillingham data to be used as example in detectRUNS.R
#'
#' @docType data
#' @usage data(beispiel)
#' @keywords datasets
#'
#' A dataset containing the 0/1/2 SNP genotypes of 17 Chillingham cattle on BTA1
#'
#'  dati The variables are as follows:
#'
#' @format A data frame with 17 rows and 2641 variables:
#' \itemize{
#'  \item IID, id of animals
#' }
#'
#'
#' A dataset containing the Map file for SNP on BTA1
#'
#'  mappa The variables are as follows:
#'
#' @format A data frame with 2635 rows and 4 variables
#' \itemize{
#'   \item V1: chromosome
#'   \item V2: SNP name
#'   \item V3: cM
#'   \item V4: bps
#' }
#'
#'
#' @examples
#' data(beispiel)
#' times <- attr(grav, "time")
#' head(mappa)
#' nrow(dati)
