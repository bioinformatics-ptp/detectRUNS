<!-- README.md is generated from README.Rmd. Please edit that file -->

# detectRUNS

detectRUNS is an R package for detecting **Runs of Homozygosity** (ROH/ROHom)
and **Runs of Heterozygosity** (ROHet) in diploid genomes. It implements two
detection methods (sliding-window and consecutive) and provides functions to
summarise, plot, and compute inbreeding coefficients from detected runs.

## Installation

```r
install.packages("detectRUNS")
```

Or install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("bioinformatics-ptp/detectRUNS/detectRUNS")
```

## Dependencies

`Imports`: ggplot2, Rcpp, gridExtra, data.table  
`Suggests`: testthat, knitr, rmarkdown  
`SystemRequirements`: OpenMP (optional, for parallel BED scanning)

## Quick start

```r
library(detectRUNS)

# Scan a PLINK BED file — auto-detects format, runs in parallel
bedFile <- system.file("extdata", "Kijas2016_Sheep_subset.bed",
                        package = "detectRUNS")
roh <- scanRUNS(bedFile,
                method       = "sliding",
                minSNP       = 15,
                minLengthBps = 100000)
print(roh)
#> ROH object  [method: sliding | type: ROHom]
#>   Samples  : 100  (2 groups)
#>   Runs     : 2388  across 2 chromosomes
#>   Length   : mean 2711583 bp  |  total 6475.26 Mbp

# Save result once — reload instantly in future sessions
saveROH(roh, "sheep_roh.roh")
roh <- loadROH("sheep_roh.roh")

# All downstream functions accept the ROH object — no file paths needed
summaryRuns(roh, Class = 2)
Froh_inbreeding(roh)
plot_Runs(roh)
plot_manhattanRuns(roh)
```

## PED/MAP input

```r
pedFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped",
                        package = "detectRUNS")
roh <- scanRUNS(pedFile, method = "consecutive",
                minSNP = 15, minLengthBps = 100000)
```

## Build an ROH object from existing results

```r
# From a CSV saved by a previous session or from readExternalRuns()
runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.sliding.csv",
                         package = "detectRUNS")
mapFile  <- system.file("extdata", "Kijas2016_Sheep_subset.map",
                         package = "detectRUNS")
runs <- readExternalRuns(runsFile, program = "detectRUNS")
roh  <- as_ROH(runs, mapFile = mapFile, method = "sliding", type = "ROHom")
```

## Documentation

See `vignette("detectRUNS.vignette", package = "detectRUNS")` for a full tutorial.
