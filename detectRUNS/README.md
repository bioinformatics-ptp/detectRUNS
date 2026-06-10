<!-- README.md is generated from README.Rmd. Please edit that file -->

# detectRUNS

detectRUNS is a R package for the detection of runs of homozygosity
(ROH/ROHom) and of heterozygosity (ROHet, a.k.a. “heterozygosity-rich
regions”) in diploid genomes. Besides runs detection, it implements
several functions to summarize and plot results.

## Installation

detectRUNS is installed as a standard R package. Some core functions are
written in C++ to increase efficieny of calculations: this makes use of
the R library Rcpp. detectRUNS uses other R packages for data
manipulation and plots. These packages are set as *Imports*, and
detectRUNS will try to install any missing packages upon installation.

## Dependencies

detectRUNS imports: plyr, iterators, itertools, ggplot2, reshape2, Rcpp,
gridExtra, data.table detectRUNS suggests: testthat, knitr, rmarkdown,
prettydoc

## Documentation

Please see the package vignette for a complete tutorial. What follows is
a minimal working example to give the gist of the tool.

## Example

This is a basic example which shows you how to detect runs of
homozygosity (ROH):

``` r
library(detectRUNS)

# Input files (bundled example data)
genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")

# 1) Detect runs with the sliding-window approach
runs <- slidingRUNS.run(
  genotypeFile, mapFile,
  windowSize = 15, threshold = 0.1,
  minSNP = 15, ROHet = FALSE,
  maxOppWindow = 1, maxMissWindow = 1,
  maxGap = 10^6, minLengthBps = 100000, minDensity = 1/10000
)

# 2) Summarise results
summaryList <- summaryRuns(
  runs = runs, mapFile = mapFile, genotypeFile = genotypeFile,
  Class = 6, snpInRuns = TRUE
)

# 3) Plot results
plot_Runs(runs = runs)
```

## Development

This repository is configured for AI-assisted development with [Claude
Code](https://claude.ai/code). A `CLAUDE.md` file at the repo root
documents the full project structure, build commands, and development
conventions for the AI assistant.

### Setting up a development environment

``` r
# 1. Install pak if not already available
install.packages("pak")

# 2. Install all package dependencies (Imports + Suggests)
pak::local_install_deps("detectRUNS", dependencies = TRUE)

# 3. Load the package for interactive development
devtools::load_all("detectRUNS")
```

### Common development tasks

``` r
devtools::document("detectRUNS")          # Regenerate documentation
devtools::test("detectRUNS")              # Run test suite
devtools::check("detectRUNS")             # Full R CMD check
Rcpp::compileAttributes("detectRUNS")     # Rebuild C++ wrappers after editing src/functions.cpp
```

### AI-assisted workflows

With Claude Code installed, the following slash commands are available
from the repo root:

- `/cran-check` — verifies version, date, NEWS.md, and runs
  `R CMD check --as-cran`
- `/code-review` — reviews the current diff for correctness and
  simplification opportunities
