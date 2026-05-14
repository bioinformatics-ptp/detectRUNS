## Merge per-dataset summary CSVs into global files under Ext_Data/results/

args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0L)
    setwd(dirname(dirname(normalizePath(script_path))))

RES_DIR  <- "Ext_Data/results"
DATASETS <- c("pigData", "SELMOL", "Innovagen_HD", "ADAPTmap")

read_if <- function(path) {
    if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else NULL
}

merge_csv <- function(filename, out_name) {
    parts <- Filter(Negate(is.null),
                    lapply(DATASETS, function(ds) read_if(file.path(RES_DIR, ds, filename))))
    if (length(parts) == 0L) { cat("  No files found for:", filename, "\n"); return(invisible(NULL)) }
    merged <- do.call(rbind, parts)
    out    <- file.path(RES_DIR, out_name)
    write.csv(merged, out, row.names = FALSE)
    cat(sprintf("  %-30s  %d rows  ->  %s\n", filename, nrow(merged), out))
    invisible(merged)
}

cat("Merging summary CSVs...\n")
master <- merge_csv("master_summary.csv", "master_summary.csv")
checks <- merge_csv("sanity_checks.csv",  "sanity_checks.csv")
timing <- merge_csv("step_timing.csv",    "step_timing.csv")

if (!is.null(master)) {
    n_ok   <- sum(master$status == "OK",         na.rm = TRUE)
    n_fail <- sum(master$status != "OK",         na.rm = TRUE)
    cat(sprintf("\nTotal scans: %d  |  OK: %d  |  non-OK: %d\n",
                nrow(master), n_ok, n_fail))
    cat(sprintf("Total time:  %.1f minutes\n", sum(master$total_s, na.rm=TRUE) / 60))
}
if (!is.null(checks)) {
    cat(sprintf("Sanity checks: %d PASS  |  %d FAIL  |  %d ERROR\n",
                sum(checks$result == "PASS"),
                sum(checks$result == "FAIL"),
                sum(checks$result == "ERROR")))
}
cat("Done.\n")
