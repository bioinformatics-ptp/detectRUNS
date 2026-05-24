
# Deprecated aliases — old ROH-naming kept for backward compatibility

#' Deprecated: use \code{\link{saveRUNS}} instead
#' @param ... passed to \code{saveRUNS}
#' @export
#' @keywords internal
saveROH <- function(...) {
  .Deprecated("saveRUNS")
  saveRUNS(...)
}

#' Deprecated: use \code{\link{loadRUNS}} instead
#' @param ... passed to \code{loadRUNS}
#' @export
#' @keywords internal
loadROH <- function(...) {
  .Deprecated("loadRUNS")
  loadRUNS(...)
}

#' Deprecated: use \code{\link{as_RUNS}} instead
#' @param ... passed to \code{as_RUNS}
#' @export
#' @keywords internal
as_ROH <- function(...) {
  .Deprecated("as_RUNS")
  as_RUNS(...)
}

#' Deprecated: use \code{\link{runsIslands}} instead
#' @param ... passed to \code{runsIslands}
#' @export
#' @keywords internal
rohIslands <- function(...) {
  .Deprecated("runsIslands")
  runsIslands(...)
}

# Suppress R CMD check NOTEs for column names used in data.table / ggplot2 NSE
utils::globalVariables(c(
  "CLASS", "group", "chrom", "freq", "id",
  "lengthBps", "nSNP", "MB",
  "CHROMOSOME", "CHR_LENGTH",
  "POPULATION", "IND", "COUNT", "START", "END",
  "SNP_NAME", "POSITION", "BREED", "PERCENTAGE",
  "value", "variable", "Froh_genome",
  "n_ROH", "total_length_bp"
))

.onAttach <- function(libname, pkgname) {
  version = packageVersion("detectRUNS")
  packageStartupMessage(paste("Using detectRUNS", version))
}

.onUnload <- function (libpath) {
  library.dynam.unload("detectRUNS", libpath)
}
