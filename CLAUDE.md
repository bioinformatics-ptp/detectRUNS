# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

detectRUNS is an R package for detecting runs of homozygosity (ROH/ROHom) and runs of heterozygosity (ROHet, "heterozygosity-rich regions") in diploid genomes. The package uses two statistical approaches:
1. **Sliding window method** (similar to PLINK, Purcell et al. 2007)
2. **Consecutive runs method** (Marras et al. 2015)

Beyond detection, it provides functions for visualization, summary statistics, and inbreeding calculations (Froh).

## Build, Test, and Development Commands

### Installing Dependencies

Before building or testing for the first time, install all declared dependencies (Imports + Suggests) with:

```R
# In R (run once after cloning the repo):
pak::local_install_deps("detectRUNS", dependencies = TRUE)
```

This installs: `plyr`, `iterators`, `itertools`, `ggplot2`, `reshape2`, `Rcpp`, `gridExtra`, `data.table` (Imports) and `testthat`, `knitr`, `rmarkdown`, `prettydoc` (Suggests). A C++ compiler must already be available for `Rcpp` to compile.

If `pak` itself is not yet installed: `install.packages("pak")`.

### Building the Package
```R
# The package source is in the detectRUNS/ subdirectory.
# All devtools commands require the path argument when run from the repo root.
# In R:
devtools::load_all("detectRUNS")              # Load package for development
devtools::build("detectRUNS")                 # Build tarball
devtools::document("detectRUNS")              # Generate NAMESPACE and documentation from roxygen2 comments
Rcpp::compileAttributes("detectRUNS")         # Regenerate R/RcppExports.R and src/RcppExports.cpp from C++ code
```

### Running Tests
```R
# Run all tests
devtools::test("detectRUNS")

# Run specific test file (path relative to package root)
devtools::test_file("detectRUNS/tests/testthat/test_run.R")
devtools::test_file("detectRUNS/tests/testthat/test_functions.R")

# In R directly
testthat::test_dir("detectRUNS/tests/testthat/")
```

### Checking Package
```R
# Full R CMD check (equivalent to GitHub Actions CI)
devtools::check("detectRUNS")

# Coverage report
covr::codecov(path="detectRUNS/")
```

### Vignette Building
```R
# Build vignettes (required before rebuilding package to include them)
devtools::build_vignettes("detectRUNS")
```

## Code Architecture

### Directory Structure
The package code is located in `detectRUNS/` subdirectory:
- **R/**: Pure R code
  - `run.R`: Two main entry functions: `slidingRUNS.run()` and `consecutiveRUNS.run()`
  - `funktionen.R`: Core utility functions for window testing, SNP classification, consecutive run detection
  - `plots.R`: All visualization functions (8+ plotting functions)
  - `Stats.R`: Statistical functions (Froh inbreeding calculation, summaries)
  - `RcppExports.R`: Auto-generated R wrappers for C++ functions (DO NOT EDIT MANUALLY)
  - `zzz.R`: Package initialization hooks

- **src/**: C++ code (compiled with Rcpp for performance)
  - `functions.cpp`: Performance-critical algorithms for genotype conversion, window sliding, consecutive run detection
  - `RcppExports.cpp`: Auto-generated C++ wrapper code (DO NOT EDIT MANUALLY)

- **tests/testthat/**: Unit tests
  - `test_run.R`: Tests for main detection functions (sliding vs consecutive)
  - `test_functions.R`: Tests for utility functions (genoConvert, readMapFile, pedConvert)
  - `test_plots.R`, `test_Stats.R`: Tests for visualization and statistics
  - `test.ped`, `test.map`, `test.raw`: Small sample genotype/map files for testing

- **vignettes/**: Documentation
  - `detectRUNS.vignette.Rmd`: Complete tutorial with examples

- **man/**: Auto-generated documentation from roxygen2

### Core Architecture

#### Data Flow
1. **Input**: Plink format files (`.ped` genotype files and `.map` SNP position files)
2. **Processing**: Two parallel detection methods:
   - **Sliding Window**: `slidingRUNS.run()` → `slidingWindow()` (R) or `slidingWindowCpp()` (C++) → `snpInRun()` → run classification
   - **Consecutive**: `consecutiveRUNS.run()` → `consecutiveRunsCpp()` (C++) → direct run detection per individual
3. **Output**: Data frame with columns: group, id, chrom, nSNP, from, to, lengthBps

#### Key Functions

**Main Detection Functions** (exported):
- `slidingRUNS.run(genotypeFile, mapFile, windowSize=15, threshold=0.05, minSNP=3, ROHet=FALSE, ...)`
  - Parameters control window size, SNP calling threshold, run length/density constraints
  - Returns data frame of detected runs
- `consecutiveRUNS.run(genotypeFile, mapFile, minSNP=15, ROHet=FALSE, ...)`
  - Window-free direct SNP-by-SNP scanning
  - Simpler parameter set than sliding window

**Utility Functions** (internal):
- `genoConvert()` / `genoConvertCpp()`: Convert 0/1/2 genotypes to 0/1 (homozygous/heterozygous)
- `pedConvertCpp()`: Convert Plink ped format allele pairs to 0/1 genotypes
- `homoZygotTest()` / `homoZygotTestCpp()`: Check if window meets homozygosity criteria
- `heteroZygotTest()`: Check if window meets heterozygosity criteria
- `slidingWindow()` / `slidingWindowCpp()`: Slide window across genome
- `snpInRun()` / `snpInRunCpp()`: Determine if SNP is in a run based on window overlap
- `readMapFile()`: Parse Plink map files
- `readPOPCpp()`: Extract population/ID info from ped files
- `createRUNdf()`: Build run data frame from consecutive method

**Visualization Functions** (exported):
- `plot_Runs()`: Plot runs per individual by chromosome
- `plot_StackedRuns()`: Stacked (piled) runs per population/chromosome
- `plot_SnpsInRuns()`: % of times each SNP is in a run
- `plot_manhattanRuns()`: Manhattan-style plot of SNP frequencies in runs
- `plot_PatternRuns()`: Number of runs vs run length (scatter/violin plots)
- `plot_ViolinRuns()`: Distribution of run lengths
- `plot_InbreedingChr()`: Inbreeding coefficient by chromosome
- `plot_DistributionRuns()`: Run length distribution by size classes

**Statistical Functions** (exported):
- `Froh_inbreeding(runs, mapFile, genome_wide=TRUE)`: Calculate genome/chromosome-wide inbreeding coefficients
- `Froh_inbreedingClass(runs, mapFile, Class=2)`: Inbreeding by run size classes (Mb intervals)
- `summaryRuns()`: Comprehensive summary statistics with optional SNP-in-runs tabulation
- `tableRuns()`: Most common runs in population (threshold-based filtering)
- `snpInsideRunsCpp()`: Count SNP occurrences in runs per population

#### Rcpp Integration

C++ functions are performance-optimized versions of R equivalents. The build system uses:
1. `Rcpp::compileAttributes()` to auto-generate R/RcppExports.R and src/RcppExports.cpp from roxygen2 directives in src/functions.cpp
2. `useDynLib(detectRUNS)` in NAMESPACE to link compiled code
3. `.Call()` mechanism in R wrappers to invoke C++ functions

**Never manually edit** RcppExports files; always edit source files and regenerate.

### Parameter Key Concepts

**Window-based Parameters**:
- `windowSize`: Number of SNPs in sliding window (typical: 10-20)
- `threshold`: Proportion of overlapping windows required to call SNP in run (0-1, typical: 0.05-0.1)
- `maxOppWindow`: Max opposite genotypes allowed in a window (e.g., heterozygotes in ROHom detection)
- `maxMissWindow`: Max missing genotypes in a window

**Run-based Parameters**:
- `minSNP`: Minimum SNPs required to define a run (typical: 3-20)
- `minLengthBps`: Minimum run length in base pairs (typical: 100kb-1Mb)
- `minDensity`: Minimum SNP density (SNPs per kbp, typical: 1/10000)
- `maxGap`: Maximum gap between consecutive SNPs in a run (typical: 1Mb)
- `maxOppRun`: Max opposite genotypes in entire run
- `maxMissRun`: Max missing genotypes in entire run
- `ROHet`: Boolean; TRUE for heterozygosity runs, FALSE for homozygosity runs (default)

### Testing Strategy

Tests use testthat framework with reference data files. Key test categories:
1. **Conversion tests**: Verify genotype conversions (0/1/2 → 0/1, ped format → 0/1)
2. **File parsing**: Map file reading, ped file parsing
3. **Detection accuracy**: Compare sliding vs consecutive results against reference CSV files
4. **Error handling**: Mismatched file sizes, missing files
5. **Utility functions**: Individual function correctness (window testing, SNP classification)

Reference files are committed (test.ROHet.sliding.csv, test.ROHet.consecutive.csv) to ensure reproducible results.

## GitHub Actions CI/CD

Two workflows defined in `.github/workflows/`:

1. **R-CMD-check.yaml**: Runs `devtools::check()` on multiple platforms
   - macOS (release)
   - Windows (release)
   - Ubuntu (devel, release, oldrel-1)
   - Uploads check artifacts on failure

2. **test-coverage.yaml**: Runs `covr::codecov()` on Ubuntu
   - Generates coverage reports for codecov.io

## Key Dependencies

**Imports** (required):
- `plyr`: Data manipulation
- `iterators`, `itertools`: Iterator utilities for window sliding
- `ggplot2`: All plotting
- `reshape2`: Data reshaping for plots
- `Rcpp`: C++ compilation and integration
- `gridExtra`: Multi-plot arrangements
- `data.table`: Fast file reading (fread)

**Suggests** (optional, for development):
- `testthat`: Unit testing
- `knitr`, `rmarkdown`: Vignette building
- `prettydoc`: Vignette HTML theme

**System Requirements**:
- R >= 3.0.0
- C++ compiler (for Rcpp compilation)

## Important Notes

1. **Always use roxygen2 comments** (lines starting with `#'`) for documentation. The NAMESPACE and man/ files are auto-generated; never edit them directly.

2. **Rcpp workflow**: When modifying `src/functions.cpp`, always run `Rcpp::compileAttributes("detectRUNS")` to regenerate wrappers.

3. **Vignette updates**: Vignettes are not rebuilt automatically. Edit `vignettes/detectRUNS.vignette.Rmd` and run `devtools::build_vignettes()` before rebuilding the package.

4. **Data formats**: Package expects standard Plink ped/map format:
   - `.ped`: Space-separated, first 6 columns are FID, IID, PID, MID, SEX, PHENOTYPE; remainder are genotypes (pairs per SNP)
   - `.map`: Space-separated, columns are CHR, SNP_NAME, cM (ignored), POSITION

5. **TODO items** (in TODO file):
   - Calculate by chromosome (cut mapfile by chrom, subset genotype)
   - Handle compressed data
   - Update Google Docs documentation
   - Write GitHub WIKI

6. **CRAN Status**: Package is on CRAN. Check cran-comments.md before submitting updates; note that installed size is 7.1Mb (4.0Mb compiled C++ libs, 2.1Mb example data) due to Rcpp and vignette data.

## Quick Reference: Common Development Tasks

| Task | Command |
|------|---------|
| Install all dependencies (first time) | `pak::local_install_deps("detectRUNS", dependencies = TRUE)` |
| Load package for testing | `devtools::load_all("detectRUNS")` |
| Run all tests | `devtools::test("detectRUNS")` |
| Run one test file | `devtools::test_file("detectRUNS/tests/testthat/test_run.R")` |
| Update documentation | `devtools::document("detectRUNS")` |
| Regenerate Rcpp wrappers | `Rcpp::compileAttributes("detectRUNS")` |
| Build vignettes | `devtools::build_vignettes("detectRUNS")` |
| Full package check | `devtools::check("detectRUNS")` |
| Coverage report | `covr::codecov(path="detectRUNS/")` |
