###########################################################
### Detect genomic runs from PLINK binary (BED/BIM/FAM) or text (PED/MAP)
###########################################################


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

#' @keywords internal
.detect_format <- function(genoFile) {
    ext  <- tolower(tools::file_ext(genoFile))
    base <- if (nchar(ext) > 0) tools::file_path_sans_ext(genoFile) else genoFile

    if (ext == "bed" || ext %in% c("bim", "fam")) {
        return(list(format = "bed", base = base))
    }
    if (ext == "ped" || ext == "map") {
        return(list(format = "ped", base = base))
    }
    # No or unknown extension — probe filesystem
    if (file.exists(paste0(base, ".bed"))) {
        return(list(format = "bed", base = base))
    }
    if (file.exists(paste0(base, ".ped"))) {
        return(list(format = "ped", base = base))
    }
    stop(sprintf(
        "Cannot determine file format for '%s'.\nSupply a path ending in .bed/.ped, or a base name when the files exist on disk.",
        genoFile
    ))
}


#' @keywords internal
.build_ped_summary <- function(runs_df, ped_file) {
    pop <- as.data.frame(readPOPCpp(ped_file))
    colnames(pop) <- c("group", "id")

    if (nrow(runs_df) == 0) {
        pop$n_ROH           <- 0L
        pop$total_length_bp <- 0L
        pop$mean_length     <- 0
        pop$n_snps_in_roh   <- 0L
        return(data.table::as.data.table(pop))
    }

    runs_df <- as.data.frame(runs_df)
    key     <- paste(runs_df$group, runs_df$id, sep = "\001")
    ukeys   <- unique(key)

    n_roh   <- tapply(runs_df$lengthBps, key, length)[ukeys]
    tot_len <- tapply(runs_df$lengthBps, key, sum)[ukeys]
    mn_len  <- tapply(runs_df$lengthBps, key, mean)[ukeys]
    n_snps  <- tapply(runs_df$nSNP,      key, sum)[ukeys]

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

    result <- merge(pop, summ, by = c("group", "id"), all.x = TRUE)
    na_idx <- is.na(result$n_ROH)
    if (any(na_idx)) {
        result$n_ROH[na_idx]           <- 0L
        result$total_length_bp[na_idx] <- 0L
        result$mean_length[na_idx]     <- 0
        result$n_snps_in_roh[na_idx]   <- 0L
    }
    data.table::as.data.table(result)
}


# ---------------------------------------------------------------------------
# Main public function
# ---------------------------------------------------------------------------

#' Detect runs of homozygosity or heterozygosity from PLINK files
#'
#' Unified entry point for ROH/ROHet detection.  Accepts either PLINK binary
#' (\code{.bed/.bim/.fam}) or PLINK text (\code{.ped/.map}) input.  The file
#' format is auto-detected from the extension of \code{genoFile}; you may also
#' pass the base name (without extension) and the function will probe the
#' filesystem (preferring BED over PED when both exist).
#'
#' Detection uses either the consecutive method (Marras et al. 2015) or the
#' sliding-window method (Bjelland et al. 2013).  The BED engine is fully
#' parallelised via OpenMP; the PED engine is single-threaded.
#'
#' @param genoFile Path to the genotype file.  Can be a \code{.bed}/\code{.ped}
#'   path or a base name without extension (e.g. \code{"mydata"} — BED is tried
#'   first, then PED).  When all three BED files are co-located with the same
#'   base name, passing only \code{genoFile} is sufficient.
#' @param bimFile Path to the \code{.bim} file (BED format only).  When
#'   \code{NULL} (default) the file is inferred from \code{genoFile}'s base
#'   name.  Specify explicitly to use a different \code{.bim} file.
#' @param famFile Path to the \code{.fam} file (BED format only).  When
#'   \code{NULL} (default) the file is inferred from \code{genoFile}'s base
#'   name.  Specify explicitly to use a different \code{.fam} file (useful for
#'   swapping population labels without re-running the scan).
#' @param mapFile Path to the \code{.map} file (PED format only).  When
#'   \code{NULL} (default) the file is inferred from \code{genoFile}'s base
#'   name.  Ignored for BED input.
#' @param method Detection method: \code{"consecutive"} (Marras 2015, default)
#'   or \code{"sliding"} (Bjelland 2013).
#' @param ROHet If \code{TRUE}, detect runs of heterozygosity (ROHet).
#'   If \code{FALSE} (default), detect runs of homozygosity (ROHom).
#' @param minSNP Minimum number of SNPs in a qualifying run.  Default 3.
#' @param maxOpp Maximum number of opposite-type genotypes allowed inside a run
#'   (consecutive) or a window (sliding).  Default 1.
#' @param maxMiss Maximum number of missing genotypes allowed inside a run
#'   (consecutive) or a window (sliding).  Default 1.
#' @param minLengthBps Minimum run length in base pairs.  Default 1000.
#' @param maxGap Maximum gap between consecutive SNPs (bp); a gap >= this value
#'   breaks a run.  Default 1e6.
#' @param windowSize Sliding-window width in SNPs (\code{method = "sliding"}
#'   only).  Default 15.
#' @param threshold Bjelland coverage-ratio threshold (strictly >); a SNP is
#'   called in a run when the fraction of overlapping passing windows exceeds
#'   this value (\code{method = "sliding"} only).  Default 0.05.
#' @param nThreads Number of parallel threads/cores.  For the BED engine this
#'   sets the OpenMP thread count; for the PED engine it sets the number of
#'   \code{mclapply} workers (Mac/Linux only — Windows always uses 1).
#'   \code{NULL} (default) or \code{0} auto-detects the number of physical
#'   cores via \code{parallel::detectCores(logical = FALSE)}.
#' @param verbose If \code{TRUE} (default), print a progress bar and a summary
#'   on completion.
#'
#' @return A named list with four elements:
#' \describe{
#'   \item{runs}{A \code{data.table} with one row per detected run and columns
#'     \code{group}, \code{id}, \code{chrom}, \code{nSNP}, \code{from},
#'     \code{to}, \code{lengthBps}.}
#'   \item{summary}{A \code{data.table} with per-individual summaries:
#'     \code{group}, \code{id}, \code{n_ROH}, \code{total_length_bp},
#'     \code{mean_length}, \code{n_snps_in_roh}.}
#'   \item{snp_freq}{Named integer vector counting how many individuals have a
#'     ROH covering each SNP (BED engine only; \code{NULL} for PED input).}
#'   \item{chrom_map}{Named integer vector mapping chromosome names to internal
#'     indices used by \code{saveRUNS()} / \code{loadRUNS()} (BED engine only;
#'     \code{NULL} for PED input).}
#' }
#'
#' @export
#'
#' @importFrom stats median
#'
#' @examples
#' \dontrun{
#' # Old-style explicit paths (bed/bim/fam) — still supported
#' res <- scanRUNS(bed, bim, fam, method = "consecutive",
#'                 minSNP = 15, maxOpp = 0, maxMiss = 0)
#'
#' # Auto-detect from .bed path (bim/fam inferred from same base name)
#' res <- scanRUNS("mydata.bed", method = "sliding",
#'                 minSNP = 15, maxOpp = 1, maxMiss = 1, minLengthBps = 100000)
#'
#' # Swap in a different fam file (e.g. updated population labels)
#' res <- scanRUNS("mydata.bed", famFile = "newlabels.fam", method = "consecutive")
#'
#' # PED input
#' res <- scanRUNS("mydata.ped", method = "consecutive",
#'                 minSNP = 15, maxOpp = 0, maxMiss = 0)
#'
#' # Use result with existing plot / statistics functions
#' plot_Runs(res$runs)
#' summaryRuns(res$runs, mapFile, genotypeFile)
#' }
scanRUNS <- function(
    genoFile,
    bimFile      = NULL,   # explicit .bim path (BED only); auto-found when NULL
    famFile      = NULL,   # explicit .fam path (BED only); auto-found when NULL
    mapFile      = NULL,   # explicit .map path (PED only); auto-found when NULL
    method       = c("consecutive", "sliding"),
    ROHet        = FALSE,
    minSNP       = 3,
    maxOpp       = 1,
    maxMiss      = 1,
    minLengthBps = 1000,
    maxGap       = 1e6,
    windowSize   = 15,
    threshold    = 0.05,
    nThreads     = NULL,
    verbose      = TRUE
) {
    method <- match.arg(method)

    # --- Resolve nThreads ---
    # NULL  → auto-detect physical cores
    # 0     → same as NULL (convenience)
    # N > 0 → use exactly N
    if (is.null(nThreads) || identical(nThreads, 0L) || identical(nThreads, 0)) {
        nThreads <- max(1L, parallel::detectCores(logical = FALSE), na.rm = TRUE)
    }
    if (!is.na(Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)))
        nThreads <- min(nThreads, 2L)
    nThreads <- as.integer(nThreads)

    # --- Parameter validation ---
    if (!is.character(genoFile) || length(genoFile) != 1L || nchar(genoFile) == 0L)
        stop("genoFile must be a non-empty file path string.\n",
             "Pass a .bed/.ped path or a base name (e.g. \"mydata\"); ",
             "bim/fam/map files are found automatically.")
    if (!is.logical(ROHet) || length(ROHet) != 1L)
        stop("ROHet must be TRUE or FALSE")
    if (!is.numeric(minSNP) || minSNP < 1L)
        stop("minSNP must be a positive integer")
    if (!is.numeric(maxOpp) || maxOpp < 0L)
        stop("maxOpp must be a non-negative integer")
    if (!is.numeric(maxMiss) || maxMiss < 0L)
        stop("maxMiss must be a non-negative integer")
    if (!is.numeric(minLengthBps) || minLengthBps < 0L)
        stop("minLengthBps must be >= 0")
    if (!is.numeric(maxGap) || maxGap < 1L)
        stop("maxGap must be >= 1")
    if (!is.numeric(windowSize) || length(windowSize) != 1L)
        stop("windowSize must be a single positive integer (n. of SNPs)")
    if (!is.numeric(threshold) || length(threshold) != 1L)
        stop("threshold must be a single number between 0 and 1")
    if (method == "sliding") {
        if (windowSize < 1L)
            stop("windowSize must be >= 1")
        if (threshold < 0 || threshold > 1)
            stop("threshold must be between 0 and 1")
    }

    # --- Auto-detect format ---
    detected <- .detect_format(genoFile)
    fmt  <- detected$format
    base <- detected$base
    ext  <- tolower(tools::file_ext(genoFile))

    if (verbose)
        message(sprintf("Scanning | format: %s | method: %s | type: %s | threads: %d",
                        toupper(fmt), method,
                        if (ROHet) "ROHet" else "ROHom",
                        nThreads))

    # --- Dispatch ---
    t0 <- proc.time()

    if (fmt == "bed") {
        bed_path <- if (ext == "bed") genoFile else paste0(base, ".bed")
        bim_path <- if (!is.null(bimFile)) bimFile else paste0(base, ".bim")
        fam_path <- if (!is.null(famFile)) famFile else paste0(base, ".fam")

        for (f in c(bed_path, bim_path, fam_path)) {
            if (!file.exists(f))
                stop(paste("File not found:", f))
        }

        result <- C_scan_roh_bed(
            bed_path      = bed_path,
            bim_path      = bim_path,
            fam_path      = fam_path,
            method        = if (method == "consecutive") 0L else 1L,
            roh_type      = if (ROHet) 1L else 0L,
            min_snps      = as.integer(minSNP),
            max_opposite  = as.integer(maxOpp),
            max_missing   = as.integer(maxMiss),
            min_length_bp = as.integer(minLengthBps),
            max_gap       = as.integer(maxGap),
            window_size   = as.integer(windowSize),
            threshold     = as.double(threshold),
            n_threads     = as.integer(nThreads),
            verbose       = isTRUE(verbose)
        )

        result$runs    <- data.table::as.data.table(result$runs)
        result$summary <- data.table::as.data.table(result$summary)

        if (verbose)
            .print_scan_summary(result$runs, result$summary, method, ROHet,
                                fmt = "BED", elapsed = (proc.time() - t0)[["elapsed"]])

        bim         <- as.data.frame(readBimFile(bim_path))
        snp_map_bed <- data.frame(CHR      = bim$chrom,
                                  SNP_NAME = bim$snp_id,
                                  POSITION = bim$bp_pos,
                                  stringsAsFactors = FALSE)
        chrom_len_bed <- .chrom_lengths_from_snp_map(snp_map_bed)
        sample_info_bed <- data.frame(
            group = result$summary$group,
            id    = result$summary$id,
            stringsAsFactors = FALSE
        )

        return(new_RUNS(
            runs          = result$runs,
            summary       = result$summary,
            chrom_lengths = chrom_len_bed,
            sample_info   = sample_info_bed,
            snp_map       = snp_map_bed,
            method        = method,
            type          = if (ROHet) "ROHet" else "ROHom",
            snp_freq      = result$snp_freq,
            chrom_map     = result$chrom_map,
            bed_path      = bed_path,
            scan_params   = list(
                input_format = "bed",
                genoFile     = normalizePath(bed_path, mustWork = FALSE),
                bimFile      = normalizePath(bim_path, mustWork = FALSE),
                famFile      = normalizePath(fam_path, mustWork = FALSE),
                mapFile      = NA_character_,
                minSNP       = as.integer(minSNP),
                maxOpp       = as.integer(maxOpp),
                maxMiss      = as.integer(maxMiss),
                minLengthBps = as.integer(minLengthBps),
                maxGap       = as.integer(maxGap),
                windowSize   = as.integer(windowSize),
                threshold    = as.double(threshold),
                ROHet        = ROHet,
                nThreads     = nThreads
            ),
            meta = list(
                pkg_version = as.character(utils::packageVersion("detectRUNS")),
                r_version   = paste(R.version$major, R.version$minor, sep = "."),
                timestamp   = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                platform    = R.version$platform
            )
        ))

    } else {
        # PED path
        ped_path <- paste0(base, ".ped")
        map_path <- if (!is.null(mapFile)) mapFile else paste0(base, ".map")

        if (!file.exists(ped_path))
            stop(paste("File not found:", ped_path))
        if (!file.exists(map_path))
            stop(sprintf(
                "Map file not found: %s\nSpecify mapFile= explicitly or place it alongside the .ped file.",
                map_path
            ))

        map_df <- readMapFile(map_path)
        colnames(map_df) <- c("Chrom", "SNP", "bps")

        if (method == "consecutive") {
            runs_df <- .consecutive_ped(
                ped_file     = ped_path,
                map_df       = map_df,
                ROHet        = ROHet,
                maxOppRun    = maxOpp,
                maxMissRun   = maxMiss,
                minSNP       = minSNP,
                minLengthBps = minLengthBps,
                maxGap       = maxGap,
                nCores       = nThreads,
                verbose      = isTRUE(verbose)
            )
        } else {
            gaps <- diff(map_df$bps)
            runs_df <- .sliding_ped(
                ped_file      = ped_path,
                map_df        = map_df,
                gaps          = gaps,
                windowSize    = windowSize,
                threshold     = threshold,
                minSNP        = minSNP,
                ROHet         = ROHet,
                maxOppWindow  = maxOpp,
                maxMissWindow = maxMiss,
                maxGap        = maxGap,
                minLengthBps  = minLengthBps,
                minDensity    = 1/1000,
                maxOppRun     = NULL,
                maxMissRun    = NULL,
                nCores        = nThreads,
                verbose       = isTRUE(verbose)
            )
        }

        runs_dt  <- data.table::as.data.table(runs_df)
        summ_dt  <- .build_ped_summary(runs_df, ped_path)

        if (verbose)
            .print_scan_summary(runs_dt, summ_dt, method, ROHet,
                                fmt = "PED", elapsed = (proc.time() - t0)[["elapsed"]])

        snp_map_ped <- data.frame(
            CHR      = map_df$Chrom,
            SNP_NAME = map_df$SNP,
            POSITION = map_df$bps,
            stringsAsFactors = FALSE
        )
        chrom_len_ped <- .chrom_lengths_from_snp_map(snp_map_ped)
        sample_info_ped <- data.frame(
            group = summ_dt$group,
            id    = summ_dt$id,
            stringsAsFactors = FALSE
        )

        return(new_RUNS(
            runs          = runs_dt,
            summary       = summ_dt,
            chrom_lengths = chrom_len_ped,
            sample_info   = sample_info_ped,
            snp_map       = snp_map_ped,
            method        = method,
            type          = if (ROHet) "ROHet" else "ROHom",
            snp_freq      = NULL,
            chrom_map     = NULL,
            scan_params   = list(
                input_format = "ped",
                genoFile     = normalizePath(ped_path, mustWork = FALSE),
                bimFile      = NA_character_,
                famFile      = NA_character_,
                mapFile      = normalizePath(map_path, mustWork = FALSE),
                minSNP       = as.integer(minSNP),
                maxOpp       = as.integer(maxOpp),
                maxMiss      = as.integer(maxMiss),
                minLengthBps = as.integer(minLengthBps),
                maxGap       = as.integer(maxGap),
                windowSize   = as.integer(windowSize),
                threshold    = as.double(threshold),
                ROHet        = ROHet,
                nThreads     = nThreads
            ),
            meta = list(
                pkg_version = as.character(utils::packageVersion("detectRUNS")),
                r_version   = paste(R.version$major, R.version$minor, sep = "."),
                timestamp   = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                platform    = R.version$platform
            )
        ))
    }
}


#' @keywords internal
.print_scan_summary <- function(runs_dt, summ_dt, method, ROHet, fmt, elapsed = NULL) {
    n_roh   <- nrow(runs_dt)
    n_indiv <- nrow(summ_dt)
    n_with  <- sum(summ_dt$n_ROH > 0)
    label   <- if (ROHet) "ROHet" else "ROHom"

    cat(sprintf("\n=== %s scan complete (%s) [%s] ===\n",
                toupper(method), label, fmt))
    cat(sprintf("  Total ROH detected:       %s\n",
                format(n_roh, big.mark = ",")))
    cat(sprintf("  Individuals with ROH:     %d / %d (%.1f%%)\n",
                n_with, n_indiv,
                if (n_indiv > 0) 100 * n_with / n_indiv else 0))

    if (n_roh > 0) {
        cat(sprintf("  Mean ROH / individual:    %.1f\n",
                    mean(summ_dt$n_ROH)))
        cat(sprintf("  Mean ROH length:          %s bp\n",
                    format(round(mean(runs_dt$lengthBps)), big.mark = ",")))
        cat(sprintf("  Median ROH length:        %s bp\n",
                    format(stats::median(runs_dt$lengthBps), big.mark = ",")))
        cat(sprintf("  Total genome in ROH:      %.2f Mbp\n",
                    sum(as.numeric(summ_dt$total_length_bp)) / 1e6))
    }
    if (!is.null(elapsed))
        cat(sprintf("  Elapsed time:             %.1f sec\n", elapsed))
    cat("\n")
}
