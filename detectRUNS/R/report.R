###############################################################################
## reportRUNS — comprehensive run-analysis report
###############################################################################

#' Generate a comprehensive report for a RUNS object
#'
#' Produces a summary report (Markdown, HTML, or PDF) from an \code{ROH} object
#' returned by \code{\link{scanRUNS}}.  The report includes an executive
#' summary, dataset overview, scan parameters, run statistics, per-chromosome
#' coverage, inbreeding coefficients, individual outlier flags, top chromosomes
#' and SNPs in runs, common ROH regions, SNP-frequency Manhattan plot, and
#' (optionally) ROH island detection results.  Works for both ROHom and ROHet
#' objects, with ROHet-specific content when appropriate.
#'
#' @param runs An \code{ROH} object returned by \code{\link{scanRUNS}} or
#'   \code{\link{as_ROH}}.
#' @param output_dir Directory where the report and plot sub-folder are written.
#'   The directory must already exist.  Default \code{"."}.
#' @param prefix Base name for all output files.  When \code{NULL} (default),
#'   auto-generated as \code{"detectRUNS_report_<YYYYMMDD_HHMMSS>"}.
#' @param islands Optional \code{ROHIslands} object from \code{\link{rohIslands}}.
#'   When supplied, an islands section and Manhattan plot are added.
#' @param format Output format: \code{"markdown"} (default), \code{"html"}, or
#'   \code{"pdf"}.  PDF requires the \pkg{rmarkdown} package and a working
#'   pandoc/LaTeX installation; if unavailable the function stops with a clear
#'   message.
#' @param Class ROH length class interval in Mb passed to
#'   \code{\link{summaryRuns}} and \code{\link{Froh_inbreedingClass}}.
#'   Default \code{2}.
#' @param snp_table If \code{TRUE} (default), compute per-SNP run frequency via
#'   \code{summaryRuns(snpInRuns = TRUE)} and add a top-10 chromosomes and
#'   top-10 SNPs table.  Set to \code{FALSE} to skip this step on large
#'   datasets.
#' @param table_threshold Frequency threshold (0–1) passed to
#'   \code{\link{tableRuns}}.  Default \code{0.50}.
#' @param froh_class If \code{TRUE} (default), include an
#'   \code{\link{Froh_inbreedingClass}} table.  Ignored when
#'   \code{runs$type == "ROHet"}.
#' @param include_plots If \code{TRUE} (default), generate and embed plots:
#'   violin (total run length), class distribution, Froh box plot, SNP-frequency
#'   Manhattan plot, and (if \code{islands} is supplied) a ROH-island Manhattan
#'   plot.
#' @param plot_width Plot width in inches.  Default \code{8}.
#' @param plot_height Plot height in inches.  Default \code{5}.
#' @param outlier_sd Number of standard deviations from the group mean used to
#'   flag individual outliers.  Default \code{2}.
#' @param overwrite If \code{FALSE} (default), stop when the output file already
#'   exists.  Set to \code{TRUE} to overwrite silently.
#' @param open If \code{TRUE}, open the report after writing
#'   (\code{\link[utils]{browseURL}} for HTML; file path printed otherwise).
#'   Silently ignored in non-interactive sessions.  Default \code{FALSE}.
#' @param verbose Print progress messages.  Default \code{TRUE}.
#'
#' @return Invisibly returns a named list:
#' \describe{
#'   \item{report_file}{Absolute path to the written report file.}
#'   \item{plot_dir}{Absolute path to the plot sub-directory, or \code{NA}.}
#'   \item{plots}{Named character vector: plot tag → absolute file path.}
#'   \item{summary}{The full \code{\link{summaryRuns}} result list.}
#' }
#'
#' @export
#' @importFrom grDevices png dev.off dev.cur
#' @importFrom utils browseURL packageVersion
#' @importFrom stats aggregate median sd
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata", "Kijas2016_Sheep_subset.bed",
#'                         package = "detectRUNS")
#' runs <- scanRUNS(bedFile, method = "sliding", minSNP = 15,
#'                  minLengthBps = 100000)
#'
#' # Markdown report in temp directory
#' out <- reportRUNS(runs, output_dir = tempdir())
#' cat("Report at:", out$report_file, "\n")
#'
#' # HTML report with ROH islands
#' islands <- rohIslands(runs, n_perm = 100, seed = 42)
#' reportRUNS(runs, islands = islands, format = "html",
#'            output_dir = tempdir(), prefix = "my_report")
#' }
reportRUNS <- function(
    runs,
    output_dir      = ".",
    prefix          = NULL,
    islands         = NULL,
    format          = c("markdown", "html", "pdf"),
    Class           = 2,
    snp_table       = TRUE,
    table_threshold = 0.50,
    froh_class      = TRUE,
    include_plots   = TRUE,
    plot_width      = 8,
    plot_height     = 5,
    outlier_sd      = 2,
    overwrite       = FALSE,
    open            = FALSE,
    verbose         = TRUE
) {
  # --------------------------------------------------------------------------
  # Validate
  # --------------------------------------------------------------------------
  if (!inherits(runs, "ROH"))
    stop("'runs' must be an ROH object from scanRUNS() or as_ROH().")
  if (!is.null(islands) && !inherits(islands, "ROHIslands"))
    stop("'islands' must be a ROHIslands object from rohIslands(), or NULL.")
  format <- match.arg(format)
  if (format == "pdf" && !requireNamespace("rmarkdown", quietly = TRUE))
    stop("format='pdf' requires the rmarkdown package.\n",
         "Install it with install.packages('rmarkdown').")
  if (!dir.exists(output_dir))
    stop(sprintf("output_dir '%s' does not exist. Create it first.", output_dir))

  is_rohet <- identical(runs$type, "ROHet")

  # --------------------------------------------------------------------------
  # Output paths
  # --------------------------------------------------------------------------
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  if (is.null(prefix))
    prefix <- paste0("detectRUNS_report_", ts)

  md_file     <- file.path(output_dir, paste0(prefix, ".md"))
  report_file <- switch(format,
    markdown = md_file,
    html     = file.path(output_dir, paste0(prefix, ".html")),
    pdf      = file.path(output_dir, paste0(prefix, ".pdf"))
  )

  if (!overwrite && file.exists(report_file))
    stop(sprintf("Output file '%s' already exists. Use overwrite = TRUE to replace it.",
                 report_file))

  plot_dir    <- NA_character_
  saved_plots <- character()

  if (include_plots) {
    plot_dir <- file.path(output_dir, paste0(prefix, "_plots"))
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  }

  .msg <- function(...) if (verbose) message(...)

  # --------------------------------------------------------------------------
  # Helper: save one plot via png device
  # --------------------------------------------------------------------------
  .save_plot <- function(tag, fn, ...) {
    if (is.na(plot_dir)) return(NULL)
    path <- file.path(plot_dir, paste0(tag, ".png"))
    tryCatch({
      grDevices::png(path,
                     width  = round(plot_width  * 120),
                     height = round(plot_height * 120),
                     res    = 120)
      tryCatch(fn(...), finally = grDevices::dev.off())
      rel <- file.path(paste0(prefix, "_plots"), paste0(tag, ".png"))
      saved_plots <<- c(saved_plots, stats::setNames(normalizePath(path), tag))
      rel
    }, error = function(e) {
      if (grDevices::dev.cur() > 1L) grDevices::dev.off()
      .msg("  Plot '", tag, "' failed: ", conditionMessage(e))
      NULL
    })
  }

  # --------------------------------------------------------------------------
  # Run analyses
  # --------------------------------------------------------------------------
  .msg("reportRUNS: running summaryRuns()")
  summ <- summaryRuns(runs, Class = Class, snpInRuns = snp_table)

  froh_gw <- NULL
  froh_cls <- NULL
  froh_breed <- data.frame(Group = character(), Mean_Froh = numeric(),
                           SD_Froh = numeric(), stringsAsFactors = FALSE)
  top_inbred <- data.frame()
  outliers   <- data.frame()

  if (!is_rohet) {
    .msg("reportRUNS: running Froh_inbreeding()")
    froh_gw <- Froh_inbreeding(runs, genome_wide = TRUE)

    if (froh_class) {
      .msg("reportRUNS: running Froh_inbreedingClass()")
      froh_cls <- Froh_inbreedingClass(runs, Class = Class)
    }

    froh_means <- aggregate(Froh_genome ~ group, data = froh_gw, FUN = mean)
    froh_sds   <- aggregate(Froh_genome ~ group, data = froh_gw, FUN = sd)
    froh_breed <- data.frame(
      Group     = froh_means$group,
      Mean_Froh = round(froh_means$Froh_genome, 4),
      SD_Froh   = round(froh_sds$Froh_genome,   4),
      stringsAsFactors = FALSE
    )
    froh_breed <- froh_breed[order(froh_breed$Mean_Froh, decreasing = TRUE), , drop = FALSE]

    top_inbred <- head(froh_gw[order(froh_gw$Froh_genome, decreasing = TRUE), ], 10)
    top_inbred$Froh_genome <- round(top_inbred$Froh_genome, 4)

    # -- Outlier detection: per-group |Froh - mean| > outlier_sd * sd
    froh_gw$z <- NA_real_
    for (grp in unique(froh_gw$group)) {
      idx <- froh_gw$group == grp
      gm  <- mean(froh_gw$Froh_genome[idx], na.rm = TRUE)
      gs  <- sd(froh_gw$Froh_genome[idx],   na.rm = TRUE)
      if (!is.na(gs) && gs > 0)
        froh_gw$z[idx] <- (froh_gw$Froh_genome[idx] - gm) / gs
    }
    out_idx  <- !is.na(froh_gw$z) & abs(froh_gw$z) > outlier_sd
    outliers <- froh_gw[out_idx, , drop = FALSE]
    outliers$Direction <- ifelse(outliers$z > 0, "High", "Low")
    outliers$Froh_genome <- round(outliers$Froh_genome, 4)
    outliers$z           <- round(outliers$z,           2)
    outliers <- outliers[order(abs(outliers$z), decreasing = TRUE), , drop = FALSE]
  }

  .msg("reportRUNS: running tableRuns()")
  tbl_runs <- tableRuns(runs, threshold = table_threshold)

  # -- Per-chromosome coverage
  r        <- as.data.frame(runs$runs)
  cl_df    <- runs$chrom_lengths
  chr_runs <- aggregate(list(n_runs = r$lengthBps,
                             total_bp = r$lengthBps),
                        by = list(chrom = r$chrom),
                        FUN = function(x) c(length(x), sum(as.numeric(x))))
  chr_n   <- aggregate(lengthBps ~ chrom, data = r, FUN = length)
  chr_tot <- aggregate(lengthBps ~ chrom, data = r,
                       FUN = function(x) sum(as.numeric(x)))
  chr_mn  <- aggregate(lengthBps ~ chrom, data = r,
                       FUN = function(x) round(mean(x) / 1e6, 3))
  chr_cov <- Reduce(function(a, b) merge(a, b, by = "chrom"), list(
    data.frame(chrom = chr_n$chrom,  N_runs    = chr_n$lengthBps,   stringsAsFactors = FALSE),
    data.frame(chrom = chr_tot$chrom, Total_Mbp = round(chr_tot$lengthBps / 1e6, 2), stringsAsFactors = FALSE),
    data.frame(chrom = chr_mn$chrom,  Mean_Mbp  = chr_mn$lengthBps, stringsAsFactors = FALSE)
  ))
  if (!is.null(cl_df) && nrow(cl_df) > 0) {
    cl_df$CHROMOSOME <- as.character(cl_df$CHROMOSOME)
    chr_cov$chrom    <- as.character(chr_cov$chrom)
    chr_cov <- merge(chr_cov,
                     data.frame(chrom = cl_df$CHROMOSOME,
                                Chr_Mbp = round(cl_df$CHR_LENGTH / 1e6, 2),
                                stringsAsFactors = FALSE),
                     by = "chrom", all.x = TRUE)
    chr_cov$Coverage_pct <- round(chr_cov$Total_Mbp / chr_cov$Chr_Mbp * 100, 1)
  }
  chr_cov$N_runs <- format(chr_cov$N_runs, big.mark = ",")
  chr_ord <- tryCatch(
    chr_cov[order(as.integer(chr_cov$chrom)), , drop = FALSE],
    error = function(e) chr_cov[order(chr_cov$chrom), , drop = FALSE]
  )

  # --------------------------------------------------------------------------
  # Generate plots
  # --------------------------------------------------------------------------
  plot_refs <- list()
  if (include_plots) {
    .msg("reportRUNS: generating plots")
    plot_refs$violin <- .save_plot("violin_sum",
                                    plot_ViolinRuns, runs, method = "sum",
                                    savePlots = FALSE)
    plot_refs$dist   <- .save_plot("dist_class",
                                    plot_DistributionRuns, runs,
                                    style = "MeanClass", Class = Class,
                                    savePlots = FALSE)
    if (!is_rohet)
      plot_refs$froh <- .save_plot("froh_boxplot",
                                    plot_InbreedingChr, runs,
                                    style = "FrohBoxPlot", savePlots = FALSE)

    plot_refs$manhattan <- .save_plot("manhattan_snp_freq",
                                       plot_manhattanRuns, runs,
                                       savePlots = FALSE)
    if (!is.null(islands))
      plot_refs$islands_manhattan <- .save_plot("roh_islands_manhattan",
                                                 plot, islands)
  }

  # --------------------------------------------------------------------------
  # Build markdown
  # --------------------------------------------------------------------------
  md <- character()

  .h1  <- function(t)    { md <<- c(md, paste0("# ",  t), "") }
  .h2  <- function(t)    { md <<- c(md, paste0("## ", t), "") }
  .h3  <- function(t)    { md <<- c(md, paste0("### ",t), "") }
  .p   <- function(...)  { md <<- c(md, paste0(...),      "") }
  .hr  <- function()     { md <<- c(md, "---",            "") }
  .bq  <- function(...)  { md <<- c(md, paste0("> ", ...), "") }
  .img <- function(rel, cap = "")
    if (!is.null(rel)) md <<- c(md, paste0("![", cap, "](", rel, ")"), "")
  .tbl <- function(df) {
    df  <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
    hdr <- paste("|", paste(names(df), collapse = " | "), "|")
    sep <- paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|")
    rows <- apply(df, 1L, function(rw) paste("|", paste(rw, collapse = " | "), "|"))
    md  <<- c(md, hdr, sep, rows, "")
  }

  run_type <- if (!is.null(runs$type)) runs$type else "RUNS"
  sec      <- 0L
  .sec <- function(title) { sec <<- sec + 1L; .h2(sprintf("%d. %s", sec, title)) }

  n_samp   <- nrow(runs$sample_info)
  n_groups <- length(unique(runs$sample_info$group))
  n_snps   <- if (!is.null(runs$snp_map)) nrow(runs$snp_map) else NA_integer_
  n_chr    <- length(unique(r$chrom))
  mean_froh_overall <- if (!is_rohet && nrow(froh_breed) > 0)
    round(mean(froh_gw$Froh_genome, na.rm = TRUE), 4) else NA

  # -- Header
  .h1(sprintf("detectRUNS — %s Analysis Report", run_type))
  .p("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  .p("**Package:** detectRUNS ", as.character(utils::packageVersion("detectRUNS")))
  if (!is.null(runs$meta))
    .p("**Scanned:** ", runs$meta$timestamp,
       "  [R ", runs$meta$r_version, " | ", runs$meta$platform, "]")
  .hr()

  # --------------------------------------------------------------------------
  # Executive summary
  # --------------------------------------------------------------------------
  .h2("Executive Summary")
  exec_rows <- list(
    c("Run type",    run_type),
    c("Method",      runs$method),
    c("Samples",     format(n_samp, big.mark = ",")),
    c("Groups",      as.character(n_groups)),
    c("Chromosomes", as.character(n_chr)),
    c("SNPs",        if (!is.na(n_snps)) format(n_snps, big.mark = ",") else "N/A"),
    c("Total runs",  format(nrow(r), big.mark = ",")),
    c("Mean run length", sprintf("%.1f kbp", mean(r$lengthBps) / 1000)),
    c("Total genome in runs", sprintf("%.3f Gbp", sum(as.numeric(r$lengthBps)) / 1e9))
  )
  if (!is_rohet)
    exec_rows <- c(exec_rows, list(
      c("Mean F_ROH (all)",    as.character(mean_froh_overall)),
      c("Outlier individuals", as.character(nrow(outliers)))
    ))
  exec_df <- data.frame(
    Metric = vapply(exec_rows, `[[`, "", 1L),
    Value  = vapply(exec_rows, `[[`, "", 2L),
    stringsAsFactors = FALSE
  )
  .tbl(exec_df)

  if (is_rohet)
    .bq("**Note:** This is a Runs of Heterozygosity (ROHet) analysis. ",
        "F_ROH inbreeding coefficients are not computed — they are only ",
        "meaningful for runs of homozygosity.")
  .hr()

  # -- 1. Scan Parameters
  .sec("Scan Parameters")
  if (!is.null(runs$scan_params)) {
    sp      <- runs$scan_params
    sliding <- identical(runs$method, "sliding")
    param_names <- c("Method", "Type", "minSNP", "maxOpp", "maxMiss",
                     "minLengthBps", "maxGap",
                     if (sliding) c("windowSize", "threshold"),
                     "ROHet", "Input format", "Threads")
    param_vals  <- c(runs$method, run_type,
                     sp$minSNP, sp$maxOpp, sp$maxMiss,
                     format(sp$minLengthBps, big.mark = ","),
                     format(sp$maxGap,       big.mark = ","),
                     if (sliding) c(sp$windowSize, sp$threshold),
                     sp$ROHet,
                     sp$input_format,
                     if (!is.null(sp$nThreads)) sp$nThreads else "N/A")
    .tbl(data.frame(Parameter = param_names, Value = as.character(param_vals),
                    stringsAsFactors = FALSE))
    if (!is.na(sp$genoFile))
      .p("> **Input file:** `", basename(sp$genoFile), "`")
  } else {
    .bq("*Scan parameters not available.",
        " Re-run `scanRUNS()` with detectRUNS >= 1.1.0 to capture them.*")
  }
  .hr()

  # -- 2. Dataset
  .sec("Dataset Overview")
  .tbl(data.frame(
    Property = c("Samples", "Groups", "Chromosomes", "SNPs"),
    Value    = c(format(n_samp, big.mark = ","), n_groups, n_chr,
                 if (!is.na(n_snps)) format(n_snps, big.mark = ",") else "N/A"),
    stringsAsFactors = FALSE
  ))
  .h3("Group composition")
  grp_cnt <- sort(table(runs$sample_info$group), decreasing = TRUE)
  .tbl(data.frame(
    Group = names(grp_cnt),
    N     = as.integer(grp_cnt),
    Pct   = sprintf("%.1f%%", as.integer(grp_cnt) / n_samp * 100),
    stringsAsFactors = FALSE
  ))
  .hr()

  # -- 3. RUNS overview
  .sec("RUNS Overview")
  .tbl(data.frame(
    Metric = c("Total runs", "Mean length", "Median length",
               "Min length", "Max length", "Total genome in runs"),
    Value  = c(format(nrow(r), big.mark = ","),
               sprintf("%.1f kbp", mean(r$lengthBps) / 1000),
               sprintf("%.1f kbp", stats::median(r$lengthBps) / 1000),
               sprintf("%.1f kbp", min(r$lengthBps)  / 1000),
               sprintf("%.1f kbp", max(r$lengthBps)  / 1000),
               sprintf("%.3f Gbp", sum(as.numeric(r$lengthBps)) / 1e9)),
    stringsAsFactors = FALSE
  ))
  .h3("Runs per group")
  grp_n    <- aggregate(lengthBps ~ group, data = r, FUN = length)
  grp_mean <- aggregate(lengthBps ~ group, data = r,
                         FUN = function(x) round(mean(x) / 1000, 1))
  grp_tot  <- aggregate(lengthBps ~ group, data = r,
                         FUN = function(x) round(sum(as.numeric(x)) / 1e6, 2))
  grp_tbl  <- Reduce(function(a, b) merge(a, b, by = "Group"), list(
    data.frame(Group = grp_n$group,    N_runs   = format(grp_n$lengthBps, big.mark=","), stringsAsFactors=FALSE),
    data.frame(Group = grp_mean$group, Mean_kbp = grp_mean$lengthBps, stringsAsFactors=FALSE),
    data.frame(Group = grp_tot$group,  Total_Mbp = grp_tot$lengthBps, stringsAsFactors=FALSE)
  ))
  if (!is_rohet && nrow(froh_breed) > 0)
    grp_tbl <- merge(grp_tbl, froh_breed[, c("Group", "Mean_Froh")],
                     by = "Group", all.x = TRUE)
  grp_tbl <- grp_tbl[order(if (!is_rohet && "Mean_Froh" %in% names(grp_tbl))
    grp_tbl$Mean_Froh else grp_tbl$Total_Mbp, decreasing = TRUE), , drop = FALSE]
  .tbl(grp_tbl)
  .img(plot_refs$violin, "Total run length per individual by group")
  .hr()

  # -- 4. Per-chromosome coverage
  .sec("Per-Chromosome Coverage")
  .tbl(chr_ord)
  .hr()

  # -- 5. Class distribution
  .sec(sprintf("Run Length Class Distribution (%d Mb bins)", Class))
  cnt <- summ$summary_ROH_count
  .tbl(cbind(Class_Mb = rownames(cnt), as.data.frame(cnt)))
  .img(plot_refs$dist, sprintf("Run length distribution by %d Mb class", Class))
  .hr()

  # -- 6. Top chromosomes and SNPs
  if (snp_table && !is.null(summ$SNPinRun)) {
    .sec("Top Chromosomes and SNPs in RUNS")
    .h3("Top 10 chromosomes by mean run length")
    mchr <- summ$summary_ROH_mean_chr
    breed_cols <- setdiff(names(mchr), "chrom")
    mchr$Mean_all_Mb <- round(rowMeans(mchr[, breed_cols, drop = FALSE],
                                        na.rm = TRUE), 2)
    top_chr <- head(mchr[order(mchr$Mean_all_Mb, decreasing = TRUE), ], 10)
    for (bc in breed_cols) top_chr[[bc]] <- round(top_chr[[bc]], 2)
    .tbl(top_chr)

    .h3("Top 10 SNPs in RUNS by frequency")
    snp_in   <- summ$SNPinRun
    top_snps <- head(snp_in[order(snp_in$PERCENTAGE, decreasing = TRUE), ], 10)
    top_snps$PERCENTAGE <- round(top_snps$PERCENTAGE, 2)
    .tbl(top_snps)

    .h3("SNP frequency across genome")
    .img(plot_refs$manhattan, "SNP-in-RUNS frequency across chromosomes")
    .hr()
  } else if (include_plots && !is.null(plot_refs$manhattan)) {
    .sec("SNP Frequency Across Genome")
    .img(plot_refs$manhattan, "SNP-in-RUNS frequency across chromosomes")
    .hr()
  }

  # -- 7. Common regions
  .sec(sprintf("Common RUNS Regions (>= %.0f%% of individuals)", table_threshold * 100))
  if (nrow(tbl_runs) > 0) {
    .tbl(tbl_runs)
  } else {
    .bq(sprintf("No regions found in >= %.0f%% of individuals.", table_threshold * 100))
  }
  .hr()

  # -- 8. Inbreeding / heterozygosity metrics
  if (is_rohet) {
    .sec("Heterozygosity Metrics")
    .bq("F_ROH inbreeding coefficients are not applicable to ROHet objects. ",
        "Use ROHet runs to characterise regions of excess heterozygosity ",
        "rather than inbreeding. Per-group run counts and lengths are reported ",
        "in Section 3 (RUNS Overview).")
    .hr()
  } else {
    .sec("Inbreeding Coefficients (F_ROH)")

    .h3("Genome-wide F_ROH by group")
    .tbl(froh_breed)

    .h3("Top 10 most inbred individuals")
    top_cols  <- intersect(c("id", "group", "Froh_genome"), names(top_inbred))
    top_print <- top_inbred[, top_cols, drop = FALSE]
    names(top_print)[names(top_print) == "Froh_genome"] <- "Froh"
    .tbl(top_print)

    if (!is.null(froh_cls)) {
      .h3(sprintf("F_ROH by %d Mb run-length class (group means)", Class))
      drop_cols <- intersect(c("id", "sum"), names(froh_cls))
      cls_data  <- froh_cls[, setdiff(names(froh_cls), drop_cols), drop = FALSE]
      cls_means <- aggregate(. ~ group, data = cls_data, FUN = mean)
      cls_means[, -1L] <- round(cls_means[, -1L, drop = FALSE], 4)
      .tbl(cls_means)
    }

    .img(plot_refs$froh, "F_ROH distribution per individual")

    # -- Outliers
    .h3(sprintf("Individual Outliers (|z| > %g SD from group mean)", outlier_sd))
    if (nrow(outliers) > 0) {
      out_cols  <- intersect(c("id", "group", "Froh_genome", "z", "Direction"),
                             names(outliers))
      out_print <- outliers[, out_cols, drop = FALSE]
      names(out_print)[names(out_print) == "Froh_genome"] <- "Froh"
      .tbl(out_print)
    } else {
      .bq(sprintf("No outliers detected at |z| > %g SD.", outlier_sd))
    }
    .hr()
  }

  # -- 9. Islands (optional)
  if (!is.null(islands)) {
    .sec("ROH Islands")
    .tbl(data.frame(
      Property = c("Method", "n_perm", "Percentile",
                   "Island SNPs", "Chromosomes affected"),
      Value    = c(runs$method,
                   islands$n_perm,
                   islands$percentile,
                   sum(islands$snp_table$is_island),
                   length(unique(islands$islands$CHR))),
      stringsAsFactors = FALSE
    ))
    .h3("Island regions")
    isl_reg <- as.data.frame(summary(islands))
    if (nrow(isl_reg) > 0) .tbl(isl_reg)
    .h3("Per-chromosome thresholds")
    thr_df <- data.frame(
      CHR       = names(islands$thresholds),
      Threshold = unname(islands$thresholds),
      Pct       = round(unname(islands$thresholds) / islands$n_samples * 100, 2),
      stringsAsFactors = FALSE
    )
    .tbl(thr_df[order(as.integer(thr_df$CHR)), , drop = FALSE])
    .img(plot_refs$islands_manhattan, "ROH island Manhattan plot")
    .hr()
  }

  # -- Methods
  .sec("Methods")
  .p("Runs of ",
     if (is_rohet) "heterozygosity" else "homozygosity",
     " were detected using the detectRUNS package with the ",
     runs$method, "-window method",
     if (!is.null(runs$scan_params))
       sprintf(" (minSNP = %d, minLen = %s bp)",
               runs$scan_params$minSNP,
               format(runs$scan_params$minLengthBps, big.mark = ","))
     else "",
     ".")
  if (!is_rohet)
    .p("Inbreeding coefficients (F_ROH) represent the proportion of the ",
       "autosomal genome covered by runs exceeding the minimum length ",
       "threshold.  Individual outliers were identified as samples whose ",
       sprintf("F_ROH deviated by more than %g SD from their group mean.", outlier_sd))
  if (!is.null(islands))
    .p("ROH islands were identified by permutation testing using ",
       islands$n_perm, " permutations at the ",
       islands$percentile * 100, "th percentile.")
  .hr()

  md <- c(md, sprintf("*Report generated by detectRUNS %s — %s*",
                      as.character(utils::packageVersion("detectRUNS")),
                      format(Sys.time(), "%Y-%m-%d")))

  # --------------------------------------------------------------------------
  # Write output
  # --------------------------------------------------------------------------
  writeLines(md, md_file)

  if (format == "html") {
    html_lines <- .md_to_html_report(
      md, title = sprintf("detectRUNS %s Report", run_type))
    writeLines(html_lines, report_file)
    .msg(sprintf("reportRUNS: HTML written to %s", report_file))

  } else if (format == "pdf") {
    rmarkdown::render(md_file, output_file = report_file,
                      output_format = "pdf_document", quiet = !verbose)
    .msg(sprintf("reportRUNS: PDF written to %s", report_file))

  } else {
    report_file <- md_file
    .msg(sprintf("reportRUNS: Markdown written to %s", report_file))
  }

  if (format != "markdown" && file.exists(md_file))
    file.remove(md_file)

  if (open && interactive())
    utils::browseURL(normalizePath(report_file))

  invisible(list(
    report_file = normalizePath(report_file),
    plot_dir    = if (!is.na(plot_dir)) normalizePath(plot_dir) else NA_character_,
    plots       = saved_plots,
    summary     = summ
  ))
}


# ---------------------------------------------------------------------------
# Internal: Markdown → HTML
# ---------------------------------------------------------------------------

#' @keywords internal
.md_to_html_report <- function(md_lines, title = "detectRUNS Report") {
  css <- paste0(
    "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;",
    "max-width:1100px;margin:40px auto;padding:0 24px;color:#333;line-height:1.6}",
    "h1,h2,h3{color:#2c3e50;margin-top:1.8em}",
    "h1{border-bottom:2px solid #3498db;padding-bottom:.4em}",
    "h2{border-bottom:1px solid #ddd;padding-bottom:.3em}",
    "table{border-collapse:collapse;width:100%;margin:1em 0;font-size:.92em}",
    "th,td{text-align:left;padding:7px 12px;border:1px solid #ddd}",
    "th{background:#f0f4f8;font-weight:600}",
    "tr:nth-child(even){background:#fafafa}",
    "img{max-width:100%;height:auto;display:block;margin:1.5em 0;",
    "border:1px solid #eee;border-radius:4px}",
    "figure{margin:1.5em 0}",
    "figcaption{font-size:.85em;color:#666;text-align:center;margin-top:.3em}",
    "blockquote{border-left:4px solid #3498db;margin:1em 0;padding:.5em 1em;",
    "background:#f8f9fa;color:#555}",
    "hr{border:none;border-top:1px solid #eee;margin:2em 0}",
    "code{background:#f5f5f5;padding:2px 5px;border-radius:3px;",
    "font-family:monospace;font-size:.9em}",
    "p{margin:.8em 0}",
    ".exec-box{background:#eaf4fb;border-left:4px solid #3498db;",
    "padding:1em 1.5em;border-radius:4px;margin:1em 0}"
  )
  header <- c(
    "<!DOCTYPE html>", "<html lang=\"en\">", "<head>",
    "<meta charset=\"UTF-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">",
    paste0("<title>", title, "</title>"),
    paste0("<style>", css, "</style>"),
    "</head>", "<body>"
  )
  footer <- c("</body>", "</html>")
  c(header, .convert_md_lines(md_lines), footer)
}


#' @keywords internal
.convert_md_lines <- function(lines) {
  result  <- character()
  i       <- 1L
  n       <- length(lines)
  in_tbl  <- FALSE
  tbl_buf <- character()

  flush_table <- function() {
    if (length(tbl_buf) >= 2L)
      result <<- c(result, .md_tbl_to_html(tbl_buf))
    tbl_buf <<- character()
    in_tbl  <<- FALSE
  }

  while (i <= n) {
    line       <- lines[[i]]
    is_tbl_row <- grepl("^\\s*\\|", line)

    if (is_tbl_row) {
      in_tbl  <- TRUE
      tbl_buf <- c(tbl_buf, line)
      i <- i + 1L
      next
    }
    if (in_tbl) flush_table()

    if (grepl("^# ",   line)) {
      result <- c(result, paste0("<h1>", .md_inline(sub("^#+ ", "", line)), "</h1>"))
      i <- i + 1L; next
    }
    if (grepl("^## ",  line)) {
      result <- c(result, paste0("<h2>", .md_inline(sub("^#+ ", "", line)), "</h2>"))
      i <- i + 1L; next
    }
    if (grepl("^### ", line)) {
      result <- c(result, paste0("<h3>", .md_inline(sub("^#+ ", "", line)), "</h3>"))
      i <- i + 1L; next
    }
    if (grepl("^---\\s*$", line)) {
      result <- c(result, "<hr>"); i <- i + 1L; next
    }
    if (grepl("^!\\[", line)) {
      m_alt <- regmatches(line, regexpr("(?<=!\\[)[^]]*", line, perl = TRUE))
      m_src <- regmatches(line, regexpr("(?<=\\()([^)]+)(?=\\))", line, perl = TRUE))
      if (length(m_src) > 0)
        result <- c(result, sprintf(
          "<figure><img src=\"%s\" alt=\"%s\"><figcaption>%s</figcaption></figure>",
          m_src, m_alt, m_alt))
      i <- i + 1L; next
    }
    if (grepl("^> ", line)) {
      result <- c(result,
                  paste0("<blockquote><p>",
                         .md_inline(sub("^> +", "", line)),
                         "</p></blockquote>"))
      i <- i + 1L; next
    }
    if (trimws(line) == "") {
      result <- c(result, ""); i <- i + 1L; next
    }
    result <- c(result, paste0("<p>", .md_inline(line), "</p>"))
    i <- i + 1L
  }

  if (in_tbl) flush_table()
  result
}


#' @keywords internal
.md_inline <- function(text) {
  text <- gsub("\\*\\*([^*]+)\\*\\*", "<strong>\\1</strong>", text)
  text <- gsub("\\*([^*]+)\\*",       "<em>\\1</em>",         text)
  text <- gsub("`([^`]+)`",           "<code>\\1</code>",     text)
  text <- gsub("~([^~]+)~",           "<sub>\\1</sub>",       text)
  text
}


#' @keywords internal
.md_tbl_to_html <- function(rows) {
  if (length(rows) < 2L) return(paste(rows, collapse = "\n"))
  parse_row <- function(r) {
    r <- gsub("^\\s*\\|\\s*|\\s*\\|\\s*$", "", trimws(r))
    strsplit(r, "\\s*\\|\\s*")[[1L]]
  }
  hdr  <- parse_row(rows[[1L]])
  data <- if (length(rows) > 2L) rows[-(1L:2L)] else character()

  out <- c("<table>", "<thead><tr>",
           paste0("<th>", vapply(hdr, .md_inline, ""), "</th>"),
           "</tr></thead>", "<tbody>")
  for (row in data) {
    cells <- parse_row(row)
    out   <- c(out, "<tr>",
               paste0("<td>", vapply(cells, .md_inline, ""), "</td>"),
               "</tr>")
  }
  paste(c(out, "</tbody>", "</table>"), collapse = "\n")
}
