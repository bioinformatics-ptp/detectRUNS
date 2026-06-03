
# detectRUNS 0.9.7

## Minor changes and bug fixes

* `tableRuns()` reimplemented in C++ (Rcpp) for improved performance; added `avg_pct` output column tracking the average percentage of samples in which each SNP is in a run
* `Froh_inbreedingClass()`: run size class bins are now exclusive intervals, fixing incorrect class label assignment (issue #41)
* `plot_manhattanRuns()`: new parameters for output file type, plot width/height, reference threshold line, and font size
* `plot_DistributionRuns()` and `summaryRuns()`: fixed class labels for run size bins (issue #41)
* `plot_InbreedingChr()` and `plot_DistributionRuns()`: fixed ggplot2 deprecation warning
* `tableRuns()`: fixed threshold filtering bug
* `snpInsideRunsCpp()`: performance improvement

# detectRUNS 0.9.3

## Major changes

* First submission to CRAN

## Bug fixes

* No bugs identified at the moment
