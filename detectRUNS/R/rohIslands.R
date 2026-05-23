###########################################################
### Permutation-based ROH island detection
###########################################################


#' Detect ROH islands via permutation testing
#'
#' Implements the permutation-based ROH island detection method of
#' Falchi et al. (2026, BMC Genomics).  For each chromosome, a null
#' distribution of SNPROH values (the number of individuals that have a given
#' SNP inside a ROH) is built by randomly permuting sample identity
#' \code{n_perm} times and re-running the same ROH scan on the permuted data.
#' The observed SNPROH of each SNP is then compared against the
#' chromosome-specific threshold derived as the \code{percentile}-th quantile
#' of the pooled null distribution.  SNPs whose real SNPROH exceeds the
#' threshold are declared ROH islands.
#'
#' Permutation is performed in C++ (OpenMP-parallel over permutations) so that
#' even 1000 permutations finish in a reasonable time.  All scan parameters
#' are taken directly from the \code{ROH} object so the permuted scans use
#' exactly the same settings as the original scan.
#'
#' @param roh        An \code{ROH} object returned by \code{\link{scanRUNS}}
#'   with BED-format input.  Must contain \code{$snp_freq}, \code{$bed_path},
#'   and \code{$scan_params} (all present when the object was created with the
#'   current package version).
#' @param bed_path   Path to the \code{.bed} file.  Only needed when
#'   \code{roh$bed_path} is \code{NULL} (objects created before this field
#'   was added); otherwise inferred automatically.
#' @param n_perm     Number of permutations.  Default 1000.  Use a smaller
#'   value (e.g. 100) for exploratory runs.
#' @param percentile Quantile used to derive the chromosome-specific threshold
#'   from the null distribution.  Default 0.99 (99th percentile, as in Falchi
#'   et al. 2026).
#' @param nThreads   Number of OpenMP threads for the permutation loop.
#'   \code{NULL} (default) auto-detects physical cores.
#' @param seed       Integer RNG seed for reproducibility.  Default 0 (random).
#' @param verbose    If \code{TRUE} (default), print per-chromosome progress.
#'
#' @return An object of class \code{"ROHIslands"}, a named list with:
#' \describe{
#'   \item{islands}{A \code{data.table} of SNPs flagged as ROH islands with
#'     columns \code{SNP_NAME}, \code{CHR}, \code{POSITION}, \code{snp_freq}
#'     (real SNPROH count), \code{threshold} (chromosome-specific threshold),
#'     \code{pct_animals} (SNPROH as a percentage of total individuals).}
#'   \item{snp_table}{A \code{data.table} with the same columns for \emph{all}
#'     SNPs (not just islands), suitable for Manhattan-style plots.}
#'   \item{thresholds}{Named numeric vector: chromosome name to threshold.}
#'   \item{n_samples}{Total number of individuals in the scan.}
#'   \item{n_perm}{Number of permutations used.}
#'   \item{percentile}{Percentile used.}
#' }
#'
#' @references
#' Falchi L, Cesarani A, Brito LF, Mastrangelo S, Pauciullo A, Macciotta NPP,
#' Gaspa G (2026). Runs of homozygosity in Italian Holstein bulls: a permutation
#' approach and time-based mapping of the genomic regions potentially under
#' selection. \emph{BMC Genomics}, 27, 203.
#' \doi{10.1186/s12864-026-12564-7}
#'
#' @seealso \code{\link{scanRUNS}}, \code{\link{plot_manhattanRuns}}
#'
#' @importFrom stats setNames
#' @importFrom tools file_path_sans_ext
#' @export
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata", "Kijas2016_Sheep_subset.bed",
#'                         package = "detectRUNS")
#' roh <- scanRUNS(bedFile, method = "consecutive",
#'                 minSNP = 15, maxOpp = 1, maxMiss = 1,
#'                 minLengthBps = 100000)
#' # Use a small n_perm for speed; increase to 1000 for publication-quality results
#' islands <- rohIslands(roh, n_perm = 100, seed = 42)
#' print(islands)
#' head(islands$islands)
#' }
rohIslands <- function(
    roh,
    bed_path   = NULL,
    n_perm     = 1000L,
    percentile = 0.99,
    nThreads   = NULL,
    seed       = 0L,
    verbose    = TRUE
) {
    # --- Input validation ---
    if (!inherits(roh, "ROH"))
        stop("'roh' must be an ROH object from scanRUNS()")
    if (is.null(roh$snp_freq))
        stop(paste(
            "'roh$snp_freq' is NULL.",
            "rohIslands() requires BED-format input (not PED).",
            "Re-run scanRUNS() with a .bed file."))
    if (is.null(roh$scan_params))
        stop(paste(
            "'roh$scan_params' is NULL.",
            "Re-run scanRUNS() with the current package version",
            "to store scan parameters in the ROH object."))

    bp <- if (!is.null(bed_path)) bed_path else roh$bed_path
    if (is.null(bp))
        stop(paste(
            "bed_path is required: either pass it explicitly as bed_path=",
            "or re-run scanRUNS() with the current package version",
            "(which stores it in roh$bed_path)."))
    if (!file.exists(bp))
        stop(paste("BED file not found:", bp))

    base     <- tools::file_path_sans_ext(bp)
    bim_path <- paste0(base, ".bim")
    fam_path <- paste0(base, ".fam")
    for (f in c(bim_path, fam_path))
        if (!file.exists(f))
            stop(paste("Required file not found:", f))

    if (!is.numeric(n_perm) || length(n_perm) != 1L || n_perm < 1L)
        stop("n_perm must be a single positive integer")
    if (!is.numeric(percentile) || length(percentile) != 1L ||
        percentile <= 0 || percentile >= 1)
        stop("percentile must be a single number strictly between 0 and 1")
    if (!is.numeric(seed) || length(seed) != 1L)
        stop("seed must be a single integer (0 = random)")

    # --- Resolve nThreads ---
    if (is.null(nThreads) || identical(nThreads, 0L) || identical(nThreads, 0))
        nThreads <- max(1L, parallel::detectCores(logical = FALSE), na.rm = TRUE)
    if (!is.na(Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)))
        nThreads <- min(nThreads, 2L)
    nThreads <- as.integer(nThreads)

    # --- Decode method and scan params from ROH object ---
    method_int <- if (roh$method == "consecutive") 0L else 1L
    roh_type   <- if (roh$type   == "ROHet")       1L else 0L
    sp         <- roh$scan_params

    n_samples <- nrow(roh$sample_info)

    if (verbose)
        message(sprintf(
            "rohIslands | method: %s | type: %s | n_perm: %d | percentile: %.3f | threads: %d",
            roh$method, roh$type, as.integer(n_perm), percentile, nThreads))

    # --- Call C++ permutation engine ---
    perm_res <- C_perm_roh_islands(
        bed_path     = bp,
        bim_path     = bim_path,
        fam_path     = fam_path,
        snp_freq_r   = as.integer(roh$snp_freq),
        method       = method_int,
        roh_type     = roh_type,
        min_snps     = sp$minSNP,
        max_opposite = sp$maxOpp,
        max_missing  = sp$maxMiss,
        min_length_bp = sp$minLengthBps,
        max_gap      = sp$maxGap,
        window_size  = sp$windowSize,
        threshold    = sp$threshold,
        n_threads    = nThreads,
        n_perm       = as.integer(n_perm),
        percentile   = as.double(percentile),
        seed         = as.integer(seed)
    )

    # --- Assemble result ---
    thresholds <- perm_res$thresholds   # named by chrom name
    is_island  <- as.logical(perm_res$is_island)

    snp_map <- roh$snp_map
    chr_col <- as.character(snp_map$CHR)
    snp_df  <- data.frame(
        SNP_NAME    = snp_map$SNP_NAME,
        CHR         = snp_map$CHR,
        POSITION    = snp_map$POSITION,
        snp_freq    = as.integer(roh$snp_freq),
        threshold   = thresholds[chr_col],
        is_island   = is_island,
        pct_animals = as.numeric(roh$snp_freq) / n_samples * 100,
        stringsAsFactors = FALSE,
        row.names = NULL
    )

    snp_table <- data.table::as.data.table(snp_df)
    islands   <- snp_table[snp_table$is_island == TRUE, ]

    if (verbose)
        message(sprintf("rohIslands | %d ROH island SNPs identified across %d chromosome(s)",
                        nrow(islands), length(unique(islands$CHR))))

    structure(
        list(
            islands    = islands,
            snp_table  = snp_table,
            thresholds = thresholds,
            n_samples  = n_samples,
            n_perm     = as.integer(n_perm),
            percentile = percentile
        ),
        class = "ROHIslands"
    )
}


#' Print a summary of an ROHIslands object
#'
#' @param x An \code{ROHIslands} object returned by \code{\link{rohIslands}}.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.ROHIslands <- function(x, ...) {
    n_isl  <- nrow(x$islands)
    n_chr  <- if (n_isl > 0L) length(unique(x$islands$CHR)) else 0L
    cat(sprintf("ROHIslands  [n_perm: %d | percentile: %.3f | individuals: %d]\n",
                x$n_perm, x$percentile, x$n_samples))
    cat(sprintf("  Island SNPs : %d  across %d chromosome%s\n",
                n_isl, n_chr, if (n_chr != 1L) "s" else ""))
    if (length(x$thresholds) > 0L) {
        cat("  Chromosome thresholds (SNPROH count):\n")
        for (nm in names(x$thresholds))
            cat(sprintf("    chr %-6s : %.0f\n", nm, x$thresholds[[nm]]))
    }
    invisible(x)
}
