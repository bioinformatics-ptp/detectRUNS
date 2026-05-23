#####################
## ROH S3 CLASS
#####################

#' @keywords internal
new_ROH <- function(runs, summary, chrom_lengths, sample_info, snp_map,
                    method, type, snp_freq = NULL, chrom_map = NULL,
                    bed_path = NULL, scan_params = NULL) {
  structure(
    list(
      runs          = runs,
      summary       = summary,
      chrom_lengths = chrom_lengths,
      sample_info   = sample_info,
      snp_map       = snp_map,
      method        = method,
      type          = type,
      snp_freq      = snp_freq,
      chrom_map     = chrom_map,
      bed_path      = bed_path,
      scan_params   = scan_params
    ),
    class = "ROH"
  )
}


#' Print a summary of an ROH object
#'
#' @param x An \code{ROH} object returned by \code{\link{scanRUNS}} or
#'   \code{\link{as_ROH}}.
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}.
#' @export
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata", "Kijas2016_Sheep_subset.bed",
#'                         package = "detectRUNS")
#' roh <- scanRUNS(bedFile, method = "sliding", minSNP = 15,
#'                 minLengthBps = 100000)
#' print(roh)
#' }
print.ROH <- function(x, ...) {
  n_samp  <- nrow(x$sample_info)
  n_group <- length(unique(x$sample_info$group))
  n_roh   <- nrow(x$runs)
  n_chr   <- length(unique(x$runs$chrom))
  cat(sprintf("ROH object  [method: %s | type: %s]\n", x$method, x$type))
  cat(sprintf("  Samples  : %d  (%d group%s)\n", n_samp, n_group,
              if (n_group != 1L) "s" else ""))
  cat(sprintf("  Runs     : %d  across %d chromosome%s\n", n_roh, n_chr,
              if (n_chr != 1L) "s" else ""))
  if (n_roh > 0L)
    cat(sprintf("  Length   : mean %.0f bp  |  total %.2f Mbp\n",
                mean(x$runs$lengthBps),
                sum(as.numeric(x$runs$lengthBps)) / 1e6))
  invisible(x)
}


#' Extract runs as a plain data.frame from an ROH object
#'
#' @param x An \code{ROH} object returned by \code{\link{scanRUNS}} or
#'   \code{\link{as_ROH}}.
#' @param ... Ignored.
#'
#' @return A \code{data.frame} with columns \code{group}, \code{id},
#'   \code{chrom}, \code{nSNP}, \code{from}, \code{to}, \code{lengthBps}.
#' @export
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata", "Kijas2016_Sheep_subset.bed",
#'                         package = "detectRUNS")
#' roh <- scanRUNS(bedFile, method = "sliding", minSNP = 15,
#'                 minLengthBps = 100000)
#' df <- as.data.frame(roh)
#' head(df)
#' }
as.data.frame.ROH <- function(x, ...) as.data.frame(x$runs)


#' Build an ROH object from pre-existing run results
#'
#' Use this when you have run results from a previous analysis (e.g. loaded via
#' \code{\link{readExternalRuns}} or saved from an earlier session) and want to
#' use the clean statistics and plot API without passing file paths repeatedly.
#'
#' @param runs A data.frame with runs (columns: group, id, chrom, nSNP, from,
#'   to, lengthBps). Output of \code{\link{readExternalRuns}} or a previously
#'   saved result.
#' @param mapFile Path to the PLINK \code{.map} file used for the original
#'   analysis (PED-format input). Provide either \code{mapFile} or
#'   \code{bedFile}.
#' @param genotypeFile Path to the PLINK \code{.ped} file. Required only if
#'   you need \code{summaryRuns(snpInRuns=TRUE)}, \code{tableRuns()},
#'   \code{plot_SnpsInRuns()}, or \code{plot_manhattanRuns()}. When omitted,
#'   \code{sample_info} is derived from the runs (individuals with zero runs
#'   will be missing from percentage denominators).
#' @param bedFile Path to the PLINK \code{.bed} file (BED-format input). The
#'   matching \code{.bim} and \code{.fam} files are found automatically.
#'   Provide either \code{bedFile} or \code{mapFile}.
#' @param method Detection method: \code{"sliding"} or \code{"consecutive"}.
#' @param type Run type: \code{"ROHom"} or \code{"ROHet"}.
#'
#' @return An \code{ROH} object usable with all downstream statistics and plot
#'   functions without supplying file paths.
#' @export
#'
#' @examples
#' \dontrun{
#' runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.sliding.csv",
#'                          package = "detectRUNS")
#' mapFile  <- system.file("extdata", "Kijas2016_Sheep_subset.map",
#'                          package = "detectRUNS")
#' pedFile  <- system.file("extdata", "Kijas2016_Sheep_subset.ped",
#'                          package = "detectRUNS")
#'
#' runs <- readExternalRuns(runsFile, program = "detectRUNS")
#' roh  <- as_ROH(runs, mapFile = mapFile, genotypeFile = pedFile,
#'                method = "sliding", type = "ROHom")
#'
#' Froh_inbreeding(roh)
#' summaryRuns(roh, Class = 2)
#' }
as_ROH <- function(runs,
                   mapFile      = NULL,
                   genotypeFile = NULL,
                   bedFile      = NULL,
                   method       = c("sliding", "consecutive"),
                   type         = c("ROHom", "ROHet")) {
  method <- match.arg(method)
  type   <- match.arg(type)
  runs   <- as.data.frame(runs)

  # --- snp_map and chrom_lengths ---
  if (!is.null(mapFile)) {
    raw_map       <- as.data.frame(readMapFile(mapFile))
    snp_map       <- data.frame(CHR      = raw_map$CHR,
                                SNP_NAME = raw_map$SNP_NAME,
                                POSITION = raw_map$POSITION,
                                stringsAsFactors = FALSE)
    chrom_lengths <- .chrom_lengths_from_snp_map(snp_map)
  } else if (!is.null(bedFile)) {
    base    <- tools::file_path_sans_ext(bedFile)
    bim     <- as.data.frame(readBimFile(paste0(base, ".bim")))
    snp_map <- data.frame(CHR      = bim$chrom,
                          SNP_NAME = bim$snp_id,
                          POSITION = bim$bp_pos,
                          stringsAsFactors = FALSE)
    chrom_lengths <- .chrom_lengths_from_snp_map(snp_map)
  } else {
    stop("Provide mapFile= (for PED input) or bedFile= (for BED input).")
  }

  # --- sample_info ---
  if (!is.null(genotypeFile)) {
    pops        <- as.data.frame(readPOPCpp(genotypeFile))
    colnames(pops) <- c("group", "id")
    sample_info <- pops
  } else if (!is.null(bedFile)) {
    base        <- tools::file_path_sans_ext(bedFile)
    fam         <- as.data.frame(readFamFile(paste0(base, ".fam")))
    sample_info <- data.frame(group = fam$fid, id = fam$iid,
                              stringsAsFactors = FALSE)
  } else {
    sample_info <- unique(runs[, c("group", "id"), drop = FALSE])
    warning(paste(
      "Neither genotypeFile nor bedFile provided.",
      "sample_info is derived from the runs data:",
      "individuals with 0 runs will be missing from SNP percentage denominators."
    ))
  }

  summary_df <- .build_summary_from_runs(runs, sample_info)

  new_ROH(
    runs          = data.table::as.data.table(runs),
    summary       = data.table::as.data.table(summary_df),
    chrom_lengths = chrom_lengths,
    sample_info   = sample_info,
    snp_map       = snp_map,
    method        = method,
    type          = type
  )
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' @keywords internal
.chrom_lengths_from_snp_map <- function(snp_map) {
  snp_map <- snp_map[snp_map$POSITION > 0, , drop = FALSE]
  chr_lengths <- tapply(snp_map$POSITION, snp_map$CHR, max)
  data.frame(
    CHROMOSOME = names(chr_lengths),
    CHR_LENGTH = as.numeric(chr_lengths),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


#' @keywords internal
.build_summary_from_runs <- function(runs, sample_info) {
  if (nrow(runs) == 0L) {
    out                 <- sample_info
    out$n_ROH           <- 0L
    out$total_length_bp <- 0L
    out$mean_length     <- 0
    out$n_snps_in_roh   <- 0L
    return(out)
  }
  key   <- paste(runs$group, runs$id, sep = "\001")
  ukeys <- unique(key)
  n_roh   <- tapply(runs$lengthBps, key, length)[ukeys]
  tot_len <- tapply(runs$lengthBps, key, sum)[ukeys]
  mn_len  <- tapply(runs$lengthBps, key, mean)[ukeys]
  n_snps  <- tapply(runs$nSNP,      key, sum)[ukeys]
  parts   <- strsplit(ukeys, "\001", fixed = TRUE)
  summ    <- data.frame(
    group           = vapply(parts, `[[`, "", 1L),
    id              = vapply(parts, `[[`, "", 2L),
    n_ROH           = as.integer(n_roh),
    total_length_bp = as.integer(tot_len),
    mean_length     = as.double(mn_len),
    n_snps_in_roh   = as.integer(n_snps),
    stringsAsFactors = FALSE
  )
  result <- merge(sample_info, summ, by = c("group", "id"), all.x = TRUE)
  na_idx <- is.na(result$n_ROH)
  if (any(na_idx)) {
    result$n_ROH[na_idx]           <- 0L
    result$total_length_bp[na_idx] <- 0L
    result$mean_length[na_idx]     <- 0
    result$n_snps_in_roh[na_idx]   <- 0L
  }
  result
}


#' Extract the runs data.frame from an ROH object or return as-is
#' @keywords internal
.get_runs <- function(x) {
  if (inherits(x, "ROH")) as.data.frame(x$runs) else as.data.frame(x)
}


#' Get chromosome lengths from an ROH object or by reading mapFile
#' @keywords internal
.get_chrom_lengths <- function(x, mapFile = NULL) {
  if (inherits(x, "ROH")) return(x$chrom_lengths)
  if (!is.null(mapFile))  return(chromosomeLength(mapFile))
  stop("Provide an ROH object from scanRUNS() / as_ROH(), or mapFile=.")
}


#' Get sample info from an ROH object or by reading genotypeFile
#' @keywords internal
.get_sample_info <- function(x, genotypeFile = NULL) {
  if (inherits(x, "ROH")) return(x$sample_info)
  if (!is.null(genotypeFile)) {
    pops <- as.data.frame(readPOPCpp(genotypeFile))
    colnames(pops) <- c("group", "id")
    return(pops)
  }
  stop("Provide an ROH object from scanRUNS() / as_ROH(), or genotypeFile=.")
}


#' Get SNP map from an ROH object or by reading mapFile
#' @keywords internal
.get_snp_map <- function(x, mapFile = NULL) {
  if (inherits(x, "ROH")) return(x$snp_map)
  if (!is.null(mapFile)) {
    raw <- as.data.frame(readMapFile(mapFile))
    return(data.frame(CHR      = raw$CHR,
                      SNP_NAME = raw$SNP_NAME,
                      POSITION = raw$POSITION,
                      stringsAsFactors = FALSE))
  }
  stop("Provide an ROH object from scanRUNS() / as_ROH(), or mapFile=.")
}
