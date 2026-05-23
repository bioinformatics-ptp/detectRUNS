###########################################################
### Binary save / load for ROH scan results
###########################################################


#' Save ROH scan results to a compact binary file
#'
#' Serialises the \code{$runs} table from a \code{scanRUNS()} result to a
#' self-contained binary file (format ROHB v1).  The file is smaller and
#' faster to read than CSV and can be loaded back with \code{loadROH()}.
#'
#' @param result Named list returned by \code{scanRUNS()}.  Must contain
#'   at least \code{$runs} and \code{$chrom_map}.
#' @param path   Output file path (e.g. \code{"results.roh"}).
#'
#' @return Invisibly returns \code{path}.
#'
#' @seealso \code{\link{loadROH}}, \code{\link{scanRUNS}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata", "Kijas2016_Sheep_subset.bed",
#'                         package = "detectRUNS")
#' roh <- scanRUNS(bedFile, method = "sliding", minSNP = 15,
#'                 minLengthBps = 100000)
#' saveROH(roh, tempfile(fileext = ".roh"))
#' }
saveROH <- function(result, path) {
    if (!is.list(result) || is.null(result$runs))
        stop("'result' must be an ROH object from scanRUNS() or as_ROH()")
    if (!is.character(path) || length(path) != 1L)
        stop("'path' must be a single file path string")

    chrom_map <- result$chrom_map
    if (is.null(chrom_map)) {
        # PED-engine results don't carry a chrom_map; synthesize a 0-based
        # integer index from the chromosome labels present in the runs.
        chroms    <- sort(unique(as.character(result$runs$chrom)))
        chrom_map <- if (length(chroms) > 0L)
            setNames(seq_along(chroms) - 1L, chroms)
        else
            setNames(integer(0L), character(0L))
    }

    C_save_roh(
        runs_df     = as.data.frame(result$runs),
        chrom_map_r = chrom_map,
        path        = path
    )
    invisible(path)
}


#' Load ROH scan results from a binary file
#'
#' Reads a binary file written by \code{saveROH()} and returns a named list
#' whose \code{$runs} element is a \code{data.table} in the same format as
#' \code{scanRUNS()$runs} (columns: group, id, chrom, nSNP, from, to,
#' lengthBps).
#'
#' @param path Path to a \code{.roh} binary file.
#'
#' @return A named list with element \code{$runs} (a \code{data.table}).
#'
#' @seealso \code{\link{saveROH}}, \code{\link{scanRUNS}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata", "Kijas2016_Sheep_subset.bed",
#'                         package = "detectRUNS")
#' roh     <- scanRUNS(bedFile, method = "sliding", minSNP = 15,
#'                     minLengthBps = 100000)
#' rohFile <- tempfile(fileext = ".roh")
#' saveROH(roh, rohFile)
#' roh2 <- loadROH(rohFile)
#' nrow(roh2$runs)
#' }
loadROH <- function(path) {
    if (!file.exists(path))
        stop(paste("File not found:", path))

    result       <- C_load_roh(path)
    result$runs  <- data.table::as.data.table(result$runs)
    result
}
