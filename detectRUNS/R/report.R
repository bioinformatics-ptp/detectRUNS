###############################################################################
## reportRUNS — comprehensive run-analysis report
###############################################################################

#' Generate a comprehensive report for a RUNS object
#'
#' Produces a summary report (Markdown, HTML, or PDF) from an \code{ROH} object
#' returned by \code{\link{scanRUNS}}.  The report includes dataset overview,
#' scan parameters, run statistics, inbreeding coefficients, top chromosomes and
#' SNPs in runs, common ROH regions, and (optionally) ROH island detection
#' results.  Works for both ROHom and ROHet objects.
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
#'   top-10 SNPs table.  Set to \code{FALSE} to skip this step on large datasets.
#' @param table_threshold Frequency threshold (0–1) passed to
#'   \code{\link{tableRuns}}.  Default \code{0.50}.
#' @param froh_class If \code{TRUE} (default), include an
#'   \code{\link{Froh_inbreedingClass}} table.
#' @param include_plots If \code{TRUE} (default), generate and embed plots:
#'   violin (total run length), class distribution, Froh box plot, and (if
#'   \code{islands} is supplied) a Manhattan plot.
#' @param plot_width Plot width in inches.  Default \code{8}.
#' @param plot_height Plot height in inches.  Default \code{5}.
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

  .msg("reportRUNS: running Froh_inbreeding()")
  froh_gw <- Froh_inbreeding(runs, genome_wide = TRUE)

  froh_cls <- NULL
  if (froh_class) {
    .msg("reportRUNS: running Froh_inbreedingClass()")
    froh_cls <- Froh_inbreedingClass(runs, Class = Class)
  }

  .msg("reportRUNS: running tableRuns()")
  tbl_runs <- tableRuns(runs, threshold = table_threshold)

  # Breed-level Froh summary
  froh_means           <- aggregate(Froh_genome ~ group, data = froh_gw, FUN = mean)
  froh_sds             <- aggregate(Froh_genome ~ group, data = froh_gw, FUN = sd)
  froh_breed           <- data.frame(
    Group      = froh_means$group,
    Mean_Froh  = round(froh_means$Froh_genome, 4),
    SD_Froh    = round(froh_sds$Froh_genome,   4),
    stringsAsFactors = FALSE
  )
  froh_breed <- froh_breed[order(froh_breed$Mean_Froh, decreasing = TRUE), , drop = FALSE]

  top_inbred <- head(froh_gw[order(froh_gw$Froh_genome, decreasing = TRUE), ], 10)
  top_inbred$Froh_genome <- round(top_inbred$Froh_genome, 4)

  # --------------------------------------------------------------------------
  # Generate plots
  # --------------------------------------------------------------------------
  plot_refs <- list()
  if (include_plots) {
    .msg("reportRUNS: generating plots")
    plot_refs$violin <- .save_plot("violin_sum",
                                    plot_ViolinRuns, runs, method = "sum", savePlots = FALSE)
    plot_refs$dist   <- .save_plot("dist_class",
                                    plot_DistributionRuns, runs,
                                    style = "MeanClass", Class = Class, savePlots = FALSE)
    plot_refs$froh   <- .save_plot("froh_boxplot",
                                    plot_InbreedingChr, runs,
                                    style = "FrohBoxPlot", savePlots = FALSE)
    if (!is.null(islands))
      plot_refs$manhattan <- .save_plot("roh_islands_manhattan",
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
  .img <- function(rel, cap = "")
    if (!is.null(rel)) md <<- c(md, paste0("![", cap, "](", rel, ")"), "")
  .tbl <- function(df) {
    df  <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
    hdr <- paste("|", paste(names(df), collapse = " | "), "|")
    sep <- paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|")
    rows <- apply(df, 1L, function(r) paste("|", paste(r, collapse = " | "), "|"))
    md  <<- c(md, hdr, sep, rows, "")
  }

  run_type <- if (!is.null(runs$type)) runs$type else "RUNS"
  sec      <- 0L
  .sec <- function(title) { sec <<- sec + 1L; .h2(sprintf("%d. %s", sec, title)) }

  # -- Header
  .h1(sprintf("detectRUNS — %s Analysis Report", run_type))
  .p("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  .p("**Package:** detectRUNS ", as.character(utils::packageVersion("detectRUNS")))
  if (!is.null(runs$meta))
    .p("**Scanned:** ", runs$meta$timestamp,
       "  [R ", runs$meta$r_version, " | ", runs$meta$platform, "]")
  .hr()

  # -- 1. Scan Parameters
  .sec("Scan Parameters")
  if (!is.null(runs$scan_params)) {
    sp <- runs$scan_params
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
    .p("> *Scan parameters not available.",
       " Re-run `scanRUNS()` with detectRUNS >= 1.1.0 to capture them.*")
  }
  .hr()

  # -- 2. Dataset
  .sec("Dataset Overview")
  n_samp   <- nrow(runs$sample_info)
  n_groups <- length(unique(runs$sample_info$group))
  n_snps   <- if (!is.null(runs$snp_map)) nrow(runs$snp_map) else NA_integer_
  n_chr    <- length(unique(runs$runs$chrom))

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
  r <- as.data.frame(runs$runs)
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
  grp_tbl  <- merge(merge(
    data.frame(Group = grp_n$group,    N_runs    = grp_n$lengthBps,    stringsAsFactors = FALSE),
    data.frame(Group = grp_mean$group, Mean_kbp  = grp_mean$lengthBps, stringsAsFactors = FALSE),
    by = "Group"), data.frame(Group = grp_tot$group, Total_Mbp = grp_tot$lengthBps,
                              stringsAsFactors = FALSE), by = "Group")
  grp_tbl  <- merge(grp_tbl, froh_breed[, c("Group", "Mean_Froh")],
                    by = "Group", all.x = TRUE)
  grp_tbl  <- grp_tbl[order(grp_tbl$Mean_Froh, decreasing = TRUE), , drop = FALSE]
  grp_tbl$N_runs <- format(grp_tbl$N_runs, big.mark = ",")
  .tbl(grp_tbl)

  .img(plot_refs$violin, "Total run length per individual by group")
  .hr()

  # -- 4. Class distribution
  .sec(sprintf("Run Length Class Distribution (%d Mb bins)", Class))
  cnt <- summ$summary_ROH_count
  .tbl(cbind(Class_Mb = rownames(cnt), as.data.frame(cnt)))
  .img(plot_refs$dist, sprintf("Run length distribution by %d Mb class", Class))
  .hr()

  # -- 5. Top chromosomes and SNPs
  if (snp_table && !is.null(summ$SNPinRun)) {
    .sec("Top Chromosomes and SNPs in RUNS")

    .h3("Top 10 chromosomes by mean run length")
    mchr <- summ$summary_ROH_mean_chr
    breed_cols <- setdiff(names(mchr), "chrom")
    mchr$Mean_all_Mb <- round(rowMeans(mchr[, breed_cols, drop = FALSE], na.rm = TRUE), 2)
    top_chr <- head(mchr[order(mchr$Mean_all_Mb, decreasing = TRUE), ], 10)
    for (bc in breed_cols) top_chr[[bc]] <- round(top_chr[[bc]], 2)
    .tbl(top_chr)

    .h3("Top 10 SNPs in RUNS by frequency")
    snp_in <- summ$SNPinRun
    top_snps <- head(snp_in[order(snp_in$PERCENTAGE, decreasing = TRUE), ], 10)
    top_snps$PERCENTAGE <- round(top_snps$PERCENTAGE, 2)
    .tbl(top_snps)
    .hr()
  }

  # -- 6. Common regions
  .sec(sprintf("Common RUNS Regions (>= %.0f%% of individuals)", table_threshold * 100))
  if (nrow(tbl_runs) > 0) {
    .tbl(tbl_runs)
  } else {
    .p(sprintf("> No regions found in >= %.0f%% of individuals.", table_threshold * 100))
  }
  .hr()

  # -- 7. Froh inbreeding
  .sec("Inbreeding Coefficients (F_ROH)")

  .h3("Genome-wide F_ROH by group")
  .tbl(froh_breed)

  .h3("Top 10 most inbred individuals")
  top_cols <- intersect(c("id", "group", "Froh_genome"), names(top_inbred))
  top_print <- top_inbred[, top_cols, drop = FALSE]
  names(top_print)[names(top_print) == "Froh_genome"] <- "Froh"
  .tbl(top_print)

  if (!is.null(froh_cls)) {
    .h3(sprintf("F_ROH by %d Mb run-length class (group means)", Class))
    drop_cols  <- intersect(c("id", "sum"), names(froh_cls))
    cls_data   <- froh_cls[, setdiff(names(froh_cls), drop_cols), drop = FALSE]
    cls_means  <- aggregate(. ~ group, data = cls_data, FUN = mean)
    cls_means[, -1L] <- round(cls_means[, -1L, drop = FALSE], 4)
    .tbl(cls_means)
  }

  .img(plot_refs$froh, "F_ROH distribution per individual")
  .hr()

  # -- 8. Islands (optional)
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

    .img(plot_refs$manhattan, "ROH island Manhattan plot")
    .hr()
  }

  # -- Methods note
  .sec("Methods")
  .p("Runs of ", if (identical(run_type, "ROHet")) "heterozygosity" else "homozygosity",
     " were detected using the detectRUNS package with the ",
     runs$method, "-window method",
     if (!is.null(runs$scan_params))
       sprintf(" (minSNP = %d, minLen = %s bp)",
               runs$scan_params$minSNP,
               format(runs$scan_params$minLengthBps, big.mark = ","))
     else "",
     ".  Inbreeding coefficients (F_ROH) represent the proportion of the ",
     "autosomal genome covered by runs exceeding the minimum length threshold.")
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
    "p{margin:.8em 0}"
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
