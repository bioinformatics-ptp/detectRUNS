<!-- README.md is generated from README.Rmd. Please edit that file -->
detectRUNS
==========

detectRUNS is a R package for the detection of runs of homozygosity (ROH/ROHom) and of heterozygosity (ROHet, a.k.a. "heterozygosity-rich regions") in diploid genomes. Besides runs detection, it implements several functions to summarize and plot results.

Installation
------------

detectRUNS is installed as a standard R package. Some core functions are written in C++ to increase efficieny of calculations: this makes use of the R library Rcpp. detectRUNS uses other R packages for data manipulation and plots. These packages are set as *Imports*, and detectRUNS will try to install any missing packages upon installation.

Dependencies
------------

detectRUNS imports: plyr, iterators, itertools, ggplot2, reshape2, Rcpp, gridExtra, data.table detectRUNS suggests: testthat, knitr, rmarkdown, prettydoc

Documentation
-------------

Please see the package vignette for a complete tutorial. What follows is a minimal working example to give the gist of the tool.

Example
-------

This is a basic example which shows you how to detect runs of homozygosity (ROH):

``` r
#1) detectRUNS (sliding-windows method)
genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
# calculating runs with sliding window approach
\dontrun{
 # skipping runs calculation
 runs <- slidingRUNS.run(genotypeFile, mapFile, windowSize = 15, threshold = 0.1,
 minSNP = 15, ROHet = FALSE,  maxOppWindow = 1, maxMissWindow = 1, maxGap=10^6,
 minLengthBps = 100000,  minDensity = 1/10000)
}
# loading pre-calculated data
runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.sliding.csv", package="detectRUNS")
colClasses <- c(rep("character", 3), rep("numeric", 4)  )
runs <- read.csv2(runsFile, header = TRUE, stringsAsFactors = FALSE,  colClasses = colClasses)

#2) summarise results
summaryList <- summaryRuns(runs = runs, mapFile = mapFilePath, genotypeFile = genotypeFilePath, Class = 6, snpInRuns = TRUE)

#3) plot results
plot_Runs(runs = runs)
```
