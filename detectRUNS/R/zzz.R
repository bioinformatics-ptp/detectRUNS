
# Deprecated aliases — old ROH-naming kept for backward compatibility
#' @export
saveROH <- function(...) {
  .Deprecated("saveRUNS")
  saveRUNS(...)
}
#' @export
loadROH <- function(...) {
  .Deprecated("loadRUNS")
  loadRUNS(...)
}
#' @export
as_ROH <- function(...) {
  .Deprecated("as_RUNS")
  as_RUNS(...)
}
#' @export
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
