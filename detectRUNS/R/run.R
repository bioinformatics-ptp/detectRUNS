###########################################################
### Compute Genomic Runs in R (homozygosity/heterozygosity)
###########################################################


# ---------------------------------------------------------------------------
# Private PED scanning helpers (called by scanRUNS)
# ---------------------------------------------------------------------------

#' @keywords internal
.consecutive_ped <- function(ped_file, map_df, ROHet,
                              maxOppRun, maxMissRun,
                              minSNP, minLengthBps, maxGap,
                              nCores = 1L, verbose = FALSE) {

    # map_df must be a plain data.frame — no data.table, safe to fork
    map_df <- as.data.frame(map_df)
    lines  <- readLines(ped_file)
    n_snp  <- nrow(map_df)
    N      <- length(lines)

    .process_one <- function(oneLine) {
        geno <- as.character(strsplit(oneLine, " ")[[1]])
        if (length(geno) - 6 != n_snp * 2)
            stop("Number of markers differ in mapFile and genotype: are those the same dataset?")
        animal <- list(FID = geno[1], IID = geno[2])
        geno   <- pedConvertCpp(geno[7:length(geno)])
        consecutiveRunsCpp(
            geno, animal,
            mapFile             = map_df,
            ROHet               = ROHet,
            minSNP              = minSNP,
            maxOppositeGenotype = maxOppRun,
            maxMiss             = maxMissRun,
            minLengthBps        = minLengthBps,
            maxGap              = maxGap
        )
    }

    # mclapply (fork) is efficient here: map_df is already a plain data.frame,
    # no data.table pointers, copy-on-write — no serialisation overhead.
    # Falls back to progress-bar lapply on Windows (no fork support).
    # Parallel path: process in chunks so the progress bar advances between batches.
    if (verbose)
        pb <- utils::txtProgressBar(min = 0, max = N, style = 3,
                                    char = "#", width = 50)
    results <- vector("list", N)

    if (nCores > 1L && .Platform$OS.type == "unix") {
        chunk_size <- max(nCores * 4L, 20L)
        chunks     <- split(seq_len(N), ceiling(seq_len(N) / chunk_size))
        done       <- 0L
        for (ch in chunks) {
            results[ch] <- parallel::mclapply(lines[ch], .process_one,
                                              mc.cores = nCores)
            done <- done + length(ch)
            if (verbose) utils::setTxtProgressBar(pb, done)
        }
    } else {
        for (i in seq_len(N)) {
            results[[i]] <- .process_one(lines[[i]])
            if (verbose) utils::setTxtProgressBar(pb, i)
        }
    }
    if (verbose) { close(pb); cat("\n") }

    RUNs <- do.call(rbind, results)
    if (is.null(RUNs))
        RUNs <- data.frame(group = character(), id = character(),
                           chrom = character(), nSNP = integer(),
                           from  = integer(),   to   = integer(),
                           lengthBps = integer())
    row.names(RUNs) <- NULL
    RUNs
}


#' @keywords internal
.sliding_ped <- function(ped_file, map_df, gaps,
                          windowSize, threshold,
                          minSNP, ROHet,
                          maxOppWindow, maxMissWindow,
                          maxGap, minLengthBps, minDensity,
                          maxOppRun, maxMissRun,
                          nCores = 1L, verbose = FALSE) {

    parameters <- list(
        windowSize    = windowSize,
        threshold     = threshold,
        minSNP        = minSNP,
        ROHet         = ROHet,
        maxOppWindow  = maxOppWindow,
        maxMissWindow = maxMissWindow,
        maxGap        = maxGap,
        minLengthBps  = minLengthBps,
        minDensity    = minDensity,
        maxOppRun     = maxOppRun,
        maxMissRun    = maxMissRun
    )

    # map_df must be a plain data.frame — no data.table, safe to fork
    map_df <- as.data.frame(map_df)
    lines  <- readLines(ped_file)
    n_snp  <- nrow(map_df)
    N      <- length(lines)

    # Validate marker count on the first line before spawning parallel workers
    # so the error propagates directly rather than being wrapped by mclapply.
    first_geno <- as.character(strsplit(lines[1L], " ")[[1]])
    if (length(first_geno) - 6L != n_snp * 2L)
        stop("Number of markers differ in mapFile and genotype: are those the same dataset?")

    .process_one <- function(oneLine) {
        geno <- as.character(strsplit(oneLine, " ")[[1]])
        animal <- list(FID = geno[1], IID = geno[2])
        geno   <- pedConvertCpp(geno[7:length(geno)])
        slidingRuns(geno, animal, map_df, gaps, parameters)
    }

    # mclapply (fork) is efficient here: map_df is already a plain data.frame,
    # no data.table pointers, copy-on-write — no serialisation overhead.
    # Falls back to progress-bar lapply on Windows (no fork support).
    # Parallel path: process in chunks so the progress bar advances between batches.
    if (verbose)
        pb <- utils::txtProgressBar(min = 0, max = N, style = 3,
                                    char = "#", width = 50)
    results <- vector("list", N)

    if (nCores > 1L && .Platform$OS.type == "unix") {
        chunk_size <- max(nCores * 4L, 20L)
        chunks     <- split(seq_len(N), ceiling(seq_len(N) / chunk_size))
        done       <- 0L
        for (ch in chunks) {
            results[ch] <- parallel::mclapply(lines[ch], .process_one,
                                              mc.cores = nCores)
            done <- done + length(ch)
            if (verbose) utils::setTxtProgressBar(pb, done)
        }
    } else {
        for (i in seq_len(N)) {
            results[[i]] <- .process_one(lines[[i]])
            if (verbose) utils::setTxtProgressBar(pb, i)
        }
    }
    if (verbose) { close(pb); cat("\n") }

    RUNs <- do.call(rbind, results)
    if (is.null(RUNs))
        RUNs <- data.frame(group = character(), id = character(),
                           chrom = character(), nSNP = integer(),
                           from  = integer(),   to   = integer(),
                           lengthBps = integer())
    row.names(RUNs) <- NULL
    RUNs
}
