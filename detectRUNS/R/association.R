################################
## RUNS-PHENOTYPE ASSOCIATION
################################

#' Test association between run presence/absence and a quantitative phenotype
#'
#' For each run region present in at least \code{minFreq} fraction of individuals,
#' fits a linear regression of the phenotype on run presence (1) vs absence (0)
#' and returns the effect size and p-value. Works with both runs of homozygosity
#' (ROHom) and runs of heterozygosity (ROHet).
#'
#' @param runs an S3 object of class \code{ROH} (from \code{\link{scanRUNS}} or
#'   \code{\link{as_ROH}}) or a plain data.frame with runs results.
#'   Must contain columns: \code{id}, \code{chrom}, \code{from}, \code{to}.
#' @param pheno data.frame with at least two columns: \code{id} (matching the
#'   individual IDs in \code{runs}, no duplicates allowed) and a numeric
#'   phenotype column. Individuals with \code{NA} in the phenotype column are
#'   excluded from regression but still count toward frequency calculations.
#' @param minFreq minimum fraction of individuals (0, 1] that must carry a run
#'   region for it to be tested. Frequency is computed over all individuals in
#'   \code{pheno}; animals absent from \code{runs} are treated as non-carriers.
#'   Default: 0.10 (10\%).
#' @param phenoCol name of the phenotype column in \code{pheno}. Default:
#'   \code{"pheno"}.
#'
#' @details
#' Run regions are labelled as \code{chr<chrom>_<from>_<to>}. Regions carried
#' by fewer than \code{minFreq * n_individuals} animals are silently dropped
#' before testing.
#'
#' The model fitted for each region is:
#' \deqn{phenotype \sim run_{presence}}
#' where \eqn{run_{presence}} is a numeric 0/1 indicator. The reported
#' \code{beta} is the slope (mean phenotype difference between carriers and
#' non-carriers) and \code{pvalue} is the two-sided Wald p-value.
#'
#' If a region has no variation among phenotyped individuals (all carriers or
#' all non-carriers), \code{beta} and \code{pvalue} are returned as \code{NA}.
#'
#' \strong{Multiple testing:} hundreds or thousands of regions are typically
#' tested simultaneously. Raw p-values should not be used directly for
#' inference; use \code{pvalue_bonferroni} (conservative, controls family-wise
#' error rate) or \code{pvalue_fdr} (less conservative, controls false
#' discovery rate via Benjamini-Hochberg).
#'
#' \strong{Independence assumption:} the linear model assumes independent
#' observations. Related animals (family structure, population stratification)
#' share runs by descent, which can inflate type I error. For structured
#' populations, consider mixed-model approaches that account for the genomic
#' relationship matrix.
#'
#' @return A data.frame with one row per tested run region and columns:
#' \describe{
#'   \item{region}{run region label (\code{chr<chrom>_<from>_<to>})}
#'   \item{n_animals}{number of phenotyped individuals carrying the region}
#'   \item{beta}{regression slope (phenotype difference carriers vs non-carriers)}
#'   \item{pvalue}{raw two-sided p-value from the linear model}
#'   \item{pvalue_bonferroni}{Bonferroni-corrected p-value}
#'   \item{pvalue_fdr}{Benjamini-Hochberg FDR-adjusted p-value}
#' }
#' Returns an empty data.frame (with the same columns) if no region passes
#' \code{minFreq}.
#'
#' @importFrom stats lm coef p.adjust
#' @export
#'
#' @examples
#' # Minimal synthetic example
#' runs <- data.frame(
#'   group     = "PopA",
#'   id        = c(paste0("ind", 1:5), paste0("ind", 1:5)),
#'   chrom     = 1,
#'   nSNP      = 20,
#'   from      = c(rep(1000L, 5), rep(6000L, 5)),
#'   to        = c(rep(5000L, 5), rep(9000L, 5)),
#'   lengthBps = c(rep(4000L, 5), rep(3000L, 5)),
#'   stringsAsFactors = FALSE
#' )
#' pheno <- data.frame(
#'   id    = paste0("ind", 1:8),
#'   pheno = c(2, 2, 2, 2, 2, 1, 1, 1),
#'   stringsAsFactors = FALSE
#' )
#' runsAssociation(runs, pheno, minFreq = 0.30)
#'
runsAssociation <- function(runs, pheno, minFreq = 0.10, phenoCol = "pheno") {

  if (!is.data.frame(pheno))
    stop("'pheno' must be a data.frame")
  if (!"id" %in% names(pheno))
    stop("'pheno' must have a column named 'id'")
  if (anyDuplicated(pheno$id))
    stop("'pheno' contains duplicate 'id' values")
  if (!phenoCol %in% names(pheno))
    stop(sprintf("Column '%s' not found in 'pheno'", phenoCol))
  if (!is.numeric(minFreq) || length(minFreq) != 1L || minFreq <= 0 || minFreq > 1)
    stop("'minFreq' must be a single numeric value in (0, 1]")

  runs_df <- .get_runs(runs)

  required <- c("id", "chrom", "from", "to")
  missing_cols <- setdiff(required, names(runs_df))
  if (length(missing_cols))
    stop(sprintf("'runs' is missing required columns: %s",
                 paste(missing_cols, collapse = ", ")))

  runs_df$region <- paste0("chr", runs_df$chrom, "_", runs_df$from, "_", runs_df$to)

  # Use pheno individuals as the reference population so that animals with
  # zero runs (absent from runs) are correctly treated as non-carriers.
  pheno_ids  <- unique(pheno$id)
  n_id       <- length(pheno_ids)
  runs_ids   <- unique(runs_df$id)
  n_matched  <- length(intersect(pheno_ids, runs_ids))
  n_runs_only  <- length(setdiff(runs_ids,  pheno_ids))
  n_pheno_only <- length(setdiff(pheno_ids, runs_ids))
  n_regions_total <- length(unique(runs_df$region))

  uniq_pairs <- unique(runs_df[, c("id", "region")])
  # count frequency only among phenotyped individuals
  carriers_in_pheno <- uniq_pairs[uniq_pairs$id %in% pheno_ids, , drop = FALSE]
  region_freq       <- table(carriers_in_pheno$region) / n_id
  valid_regions     <- names(region_freq)[region_freq >= minFreq]

  if (length(valid_regions) == 0L) {
    message("No run regions pass the minFreq threshold.")
    return(data.frame(region            = character(),
                      n_animals         = integer(),
                      beta              = numeric(),
                      pvalue            = numeric(),
                      pvalue_bonferroni = numeric(),
                      pvalue_fdr        = numeric(),
                      stringsAsFactors  = FALSE))
  }

  # Build presence/absence matrix with ALL pheno individuals as rows so that
  # non-carriers (not appearing in runs) receive 0 for every region.
  mat <- matrix(0L,
                nrow     = n_id,
                ncol     = length(valid_regions),
                dimnames = list(pheno_ids, valid_regions))

  idx_row <- match(uniq_pairs$id,     pheno_ids)
  idx_col <- match(uniq_pairs$region, valid_regions)
  keep    <- !is.na(idx_row) & !is.na(idx_col)
  mat[cbind(idx_row[keep], idx_col[keep])] <- 1L

  runs_wide      <- as.data.frame(mat, stringsAsFactors = FALSE)
  runs_wide$id   <- rownames(runs_wide)
  rownames(runs_wide) <- NULL

  ok <- merge(pheno[, c("id", phenoCol), drop = FALSE], runs_wide, by = "id")

  if (nrow(ok) == 0L)
    stop("No individuals match between 'pheno' and 'runs' on column 'id'")

  results <- lapply(valid_regions, function(col) {
    x          <- ok[[col]]
    y          <- ok[[phenoCol]]
    valid_rows <- !is.na(y)
    n_valid    <- sum(valid_rows)
    n_present  <- sum(x[valid_rows] == 1L)

    if (n_present == 0L || n_present == n_valid)
      return(data.frame(region    = col,
                        n_animals = n_present,
                        beta      = NA_real_,
                        pvalue    = NA_real_,
                        stringsAsFactors = FALSE))

    # suppressWarnings: lm warns on perfect separation (all carriers have the
    # same phenotype value); the estimates are still correct in that case.
    fit   <- suppressWarnings(lm(y ~ x))
    coefs <- suppressWarnings(coef(summary(fit)))
    data.frame(region    = col,
               n_animals = n_present,
               beta      = coefs[2L, "Estimate"],
               pvalue    = coefs[2L, "Pr(>|t|)"],
               stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out$pvalue_bonferroni <- p.adjust(out$pvalue, method = "bonferroni")
  out$pvalue_fdr        <- p.adjust(out$pvalue, method = "BH")

  n_sig_bonf <- sum(out$pvalue_bonferroni < 0.05, na.rm = TRUE)
  n_sig_fdr  <- sum(out$pvalue_fdr        < 0.05, na.rm = TRUE)

  message("--- runsAssociation summary ---")
  message(sprintf("  Animals in runs       : %d", length(runs_ids)))
  message(sprintf("  Animals in pheno      : %d", n_id))
  message(sprintf("  Animals matched       : %d", n_matched))
  if (n_runs_only  > 0L) message(sprintf("  In runs only (ignored): %d", n_runs_only))
  if (n_pheno_only > 0L) message(sprintf("  In pheno only (non-carriers): %d", n_pheno_only))
  message(sprintf("  Run regions detected  : %d", n_regions_total))
  message(sprintf("  Run regions tested    : %d  (minFreq >= %.2f)", length(valid_regions), minFreq))
  message(sprintf("  Significant (Bonferroni p<0.05): %d / %d", n_sig_bonf, length(valid_regions)))
  message(sprintf("  Significant (FDR p<0.05)       : %d / %d", n_sig_fdr,  length(valid_regions)))

  out
}
