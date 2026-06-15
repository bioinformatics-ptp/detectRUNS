## Tests for plotting functions
## All plots are rendered to a tempfile pdf so no Rplots.pdf is created.

library(testthat)
library(detectRUNS)

context("Plotting functions")

# Fixtures (local test files — fast, CRAN-safe)
ped_file  <- "test.ped"
map_file  <- "test.map"
runs_file <- "test.ROHet.sliding.csv"

.runs_df <- suppressMessages(
    readExternalRuns(inputFile = runs_file, program = "detectRUNS")
)

.roh_plot <- suppressWarnings(suppressMessages(
    scanRUNS(ped_file, method = "consecutive",
             minSNP = 10, maxOpp = 1, maxMiss = 1,
             minLengthBps = 50000, verbose = FALSE)
))

# Helper: redirect any stray graphics output to a tempfile
.with_null_dev <- function(expr) {
    tmp <- tempfile(fileext = ".pdf")
    grDevices::pdf(tmp)
    on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
    force(expr)
}

# ---------------------------------------------------------------------------
# plot_manhattanRuns
# ---------------------------------------------------------------------------
test_that("plot_manhattanRuns produces a file in tempdir without error", {
    out <- tempfile()
    .with_null_dev(
        expect_no_error(
            suppressWarnings(suppressMessages(
                plot_manhattanRuns(.runs_df, ped_file, map_file,
                                   savePlots = TRUE, outputName = out)
            ))
        )
    )
})

# ---------------------------------------------------------------------------
# plot_Runs — prints to device, does not return a value
# ---------------------------------------------------------------------------
test_that("plot_Runs runs without error for non-empty results", {
    skip_if(nrow(.roh_plot$runs) == 0L, "no runs in test dataset")
    .with_null_dev(
        expect_no_error(suppressWarnings(suppressMessages(
            plot_Runs(.roh_plot)
        )))
    )
})

# ---------------------------------------------------------------------------
# plot_DistributionRuns
# ---------------------------------------------------------------------------
test_that("plot_DistributionRuns runs without error", {
    skip_if(nrow(.roh_plot$runs) == 0L, "no runs in test dataset")
    .with_null_dev(
        expect_no_error(suppressWarnings(suppressMessages(
            plot_DistributionRuns(.roh_plot)
        )))
    )
})

# ---------------------------------------------------------------------------
# plot_InbreedingChr
# ---------------------------------------------------------------------------
test_that("plot_InbreedingChr runs without error", {
    .with_null_dev(
        expect_no_error(suppressWarnings(suppressMessages(
            plot_InbreedingChr(.roh_plot)
        )))
    )
})
