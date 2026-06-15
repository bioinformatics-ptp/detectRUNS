###########################################################
### Permutation-based runs island detection
###########################################################


#' Detect runs islands via permutation testing
#'
#' Implements the permutation-based runs island detection method of
#' Falchi et al. (2026, BMC Genomics).  For each chromosome, a null
#' distribution of SNPROH values (the number of individuals that have a given
#' SNP inside a ROH) is built by randomly permuting sample identity
#' \code{n_perm} times and re-running the same runs scan on the permuted data.
#' The observed SNPROH of each SNP is then compared against the
#' chromosome-specific threshold derived as the \code{percentile}-th quantile
#' of the pooled null distribution.  SNPs whose real SNPROH exceeds the
#' threshold are declared runs islands.
#'
#' Permutation is performed in C++ (OpenMP-parallel over permutations) so that
#' even 1000 permutations finish in a reasonable time.  All scan parameters
#' are taken directly from the \code{RUNS} object so the permuted scans use
#' exactly the same settings as the original scan.
#'
#' @param roh        A \code{RUNS} object returned by \code{\link{scanRUNS}}
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
#' @return An object of class \code{"RunsIslands"}, a named list with:
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
#' islands <- runsIslands(roh, n_perm = 100, seed = 42)
#' print(islands)
#' head(islands$islands)
#' }
runsIslands <- function(
    roh,
    bed_path   = NULL,
    n_perm     = 1000L,
    percentile = 0.99,
    nThreads   = NULL,
    seed       = 0L,
    verbose    = TRUE
) {
    # --- Input validation ---
    if (!inherits(roh, "RUNS"))
        stop("'roh' must be a RUNS object from scanRUNS()")
    if (is.null(roh$snp_freq))
        stop(paste(
            "'roh$snp_freq' is NULL.",
            "runsIslands() requires BED-format input (not PED).",
            "Re-run scanRUNS() with a .bed file."))
    if (is.null(roh$scan_params))
        stop(paste(
            "'roh$scan_params' is NULL.",
            "Re-run scanRUNS() with the current package version",
            "to store scan parameters in the RUNS object."))

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

    # --- Decode method and scan params from RUNS object ---
    method_int <- if (roh$method == "consecutive") 0L else 1L
    roh_type   <- if (roh$type   == "ROHet")       1L else 0L
    sp         <- roh$scan_params

    n_samples <- nrow(roh$sample_info)

    if (verbose)
        message(sprintf(
            "runsIslands | method: %s | type: %s | n_perm: %d | percentile: %.3f | threads: %d",
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
        message(sprintf("runsIslands | %d runs island SNPs identified across %d chromosome(s)",
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
        class = "RunsIslands"
    )
}


#' Print a summary of an RunsIslands object
#'
#' @param x An \code{RunsIslands} object returned by \code{\link{runsIslands}}.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.RunsIslands <- function(x, ...) {
    n_isl  <- nrow(x$islands)
    n_chr  <- if (n_isl > 0L) length(unique(x$islands$CHR)) else 0L
    cat(sprintf("RunsIslands  [n_perm: %d | percentile: %.3f | individuals: %d]\n",
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


#' Summarise runs islands as contiguous genomic regions
#'
#' Groups consecutive island SNPs (adjacent in BIM order, on the same
#' chromosome) into contiguous regions and returns one row per region with
#' start/end coordinates, SNP count, peak SNPROH percentage, and region width.
#'
#' @param object An \code{RunsIslands} object returned by
#'   \code{\link{runsIslands}}.
#' @param ...    Ignored.
#'
#' @return A \code{data.table} with columns:
#' \describe{
#'   \item{CHR}{Chromosome name.}
#'   \item{start_bp}{Start position of the island region (bp).}
#'   \item{end_bp}{End position of the island region (bp).}
#'   \item{n_snps}{Number of island SNPs in the region.}
#'   \item{peak_pct}{Highest SNPROH percentage observed in the region.}
#'   \item{width_mb}{Region width in megabases (\code{end_bp - start_bp}).}
#' }
#'
#' @seealso \code{\link{runsIslands}}, \code{\link{plot.RunsIslands}}
#' @export
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata", "Kijas2016_Sheep_subset.bed",
#'                         package = "detectRUNS")
#' roh     <- scanRUNS(bedFile, method = "consecutive", minSNP = 15,
#'                     maxOpp = 1, maxMiss = 1, minLengthBps = 100000)
#' islands <- runsIslands(roh, n_perm = 100, seed = 42)
#' summary(islands)
#' }
summary.RunsIslands <- function(object, ...) {
    df <- as.data.frame(object$snp_table)   # BIM order preserved

    empty <- data.table::data.table(
        CHR      = character(),
        start_bp = integer(),
        end_bp   = integer(),
        n_snps   = integer(),
        peak_pct = numeric(),
        width_mb = numeric()
    )

    if (nrow(df) == 0L || sum(df$is_island) == 0L) return(empty)

    # Assign a region ID: increment when transitioning FALSE→TRUE or chr changes.
    chr_vec     <- as.character(df$CHR)
    isl_vec     <- df$is_island
    chr_change  <- c(TRUE, chr_vec[-1] != chr_vec[-length(chr_vec)])
    prev_false  <- c(TRUE, !isl_vec[-length(isl_vec)])
    new_region  <- isl_vec & (prev_false | chr_change)
    df$region_id <- cumsum(new_region)
    df$region_id[!isl_vec] <- NA_integer_

    island_df <- df[!is.na(df$region_id), ]

    regions <- lapply(split(island_df, island_df$region_id), function(g) {
        data.frame(
            CHR      = as.character(g$CHR[1]),
            start_bp = min(g$POSITION),
            end_bp   = max(g$POSITION),
            n_snps   = nrow(g),
            peak_pct = round(max(g$pct_animals), 2),
            width_mb = round((max(g$POSITION) - min(g$POSITION)) / 1e6, 3),
            stringsAsFactors = FALSE
        )
    })

    out <- do.call(rbind, regions)
    out <- out[order(out$CHR, out$start_bp), ]
    row.names(out) <- NULL
    data.table::as.data.table(out)
}


#' Manhattan plot of runs island detection results
#'
#' Plots SNPROH frequency (percentage of individuals with a given SNP inside a
#' runs) across the genome, with chromosome-specific permutation thresholds
#' shown as dashed lines and island SNPs highlighted in a distinct colour.
#' The plot mirrors Fig. 1 of Falchi et al. (2026, \emph{BMC Genomics}).
#'
#' @param x          An \code{RunsIslands} object returned by
#'   \code{\link{runsIslands}}.
#' @param col_island Colour for island SNPs.  Default \code{"firebrick"}.
#' @param col_snp    Two-element character vector of alternating colours for
#'   non-island SNPs (one per chromosome, alternating).
#'   Default \code{c("grey60", "grey80")}.
#' @param col_threshold Colour for the per-chromosome threshold lines.
#'   Default \code{"steelblue"}.
#' @param title      Plot title.  Default \code{"Runs Island Detection"}.
#' @param pt_size    Point size for SNPs.  Default \code{0.6}.
#' @param pt_alpha   Point transparency.  Default \code{0.8}.
#' @param ...        Ignored.
#'
#' @return A \code{ggplot2} object (invisible).  The plot is printed as a
#'   side-effect.
#'
#' @seealso \code{\link{runsIslands}}, \code{\link{summary.RunsIslands}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_segment scale_x_continuous scale_y_continuous expansion labs theme_bw theme element_blank element_text
#' @export
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata", "Kijas2016_Sheep_subset.bed",
#'                         package = "detectRUNS")
#' roh     <- scanRUNS(bedFile, method = "consecutive", minSNP = 15,
#'                     maxOpp = 1, maxMiss = 1, minLengthBps = 100000)
#' islands <- runsIslands(roh, n_perm = 100, seed = 42)
#' plot(islands)
#' }
plot.RunsIslands <- function(
    x,
    col_island     = "firebrick",
    col_snp        = c("grey60", "grey80"),
    col_threshold  = "steelblue",
    title          = "Runs Island Detection",
    pt_size        = 0.6,
    pt_alpha       = 0.8,
    ...)
{
    df   <- as.data.frame(x$snp_table)
    chrs <- unique(as.character(df$CHR))   # BIM order

    # Compute cumulative x-axis offsets (small gap of 2% max_pos between chroms)
    chr_maxpos <- sapply(chrs, function(ch)
        max(df$POSITION[as.character(df$CHR) == ch]))
    gap     <- sum(chr_maxpos) * 0.02 / max(length(chrs) - 1L, 1L)
    offsets <- c(0, cumsum(chr_maxpos[-length(chr_maxpos)] + gap))
    names(offsets) <- chrs

    df$cum_pos  <- df$POSITION + offsets[as.character(df$CHR)]
    chr_idx     <- match(as.character(df$CHR), chrs)
    df$pt_color <- ifelse(df$is_island,
                          col_island,
                          col_snp[(chr_idx %% 2L) + 1L])

    # Chromosome midpoints for x labels
    chr_mids <- sapply(chrs, function(ch) {
        pos <- df$cum_pos[as.character(df$CHR) == ch]
        (min(pos) + max(pos)) / 2
    })

    # Per-chromosome threshold segments (convert count → % animals)
    thr_segs <- do.call(rbind, lapply(chrs, function(ch) {
        pos <- df$cum_pos[as.character(df$CHR) == ch]
        thr <- if (!is.null(x$thresholds[ch]) && !is.na(x$thresholds[ch]))
            x$thresholds[[ch]] / x$n_samples * 100
        else NA_real_
        data.frame(xmin = min(pos), xmax = max(pos),
                   thr = thr, stringsAsFactors = FALSE)
    }))
    thr_segs <- thr_segs[!is.na(thr_segs$thr), ]

    p <- ggplot2::ggplot(df,
             ggplot2::aes(x = .data[["cum_pos"]] / 1e6,
                          y = .data[["pct_animals"]])) +
        ggplot2::geom_point(color = df$pt_color,
                            size  = pt_size,
                            alpha = pt_alpha) +
        ggplot2::geom_segment(
            data = thr_segs,
            ggplot2::aes(x    = .data[["xmin"]] / 1e6, xend = .data[["xmax"]] / 1e6,
                         y    = .data[["thr"]],         yend = .data[["thr"]]),
            color     = col_threshold,
            linewidth = 0.8,
            linetype  = "dashed",
            inherit.aes = FALSE) +
        ggplot2::scale_x_continuous(
            breaks = chr_mids / 1e6,
            labels = chrs,
            expand = c(0.01, 0)) +
        ggplot2::scale_y_continuous(
            limits = c(0, NA),
            expand = ggplot2::expansion(mult = c(0, 0.05))) +
        ggplot2::labs(
            title   = title,
            x       = "Chromosome",
            y       = "% individuals with SNP in ROH",
            caption = sprintf("n_perm = %d  |  threshold percentile = %.2f  |  n = %d individuals",
                              x$n_perm, x$percentile, x$n_samples)) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
            panel.grid.minor   = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_blank(),
            plot.title         = ggplot2::element_text(hjust = 0.5),
            plot.caption       = ggplot2::element_text(size = 8, colour = "grey50")
        )

    print(p)
    invisible(p)
}
