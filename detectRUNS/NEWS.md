
# detectRUNS 0.9.7

## New features

* New `classifyRuns()` method to bin ROH by size classes, replacing the previous ad-hoc approach.
* `tableRuns()` reimplemented in C++ (`tableRunsCpp`) for improved performance; threshold filtering now applied before per-breed evaluation.
* `snpInsideRunsCpp()` now returns a `Number` column; population counts computed outside the C++ loop for efficiency.
* `plot_manhattanRuns()` gains parameters to control output file type, plot width/height, reference threshold, and font size.
* `summaryRuns()` output now includes corrected mean chromosome ROH dataset descriptions.
* Added GitHub Actions CI/CD workflows (R CMD check, test coverage); removed legacy Travis CI configuration.
* Added pull request and issue templates.

## Bug fixes

* Fixed R CMD check warnings: C++ format specifier, Rd math notation, and `.Rbuildignore` entries.
* Fixed class labels in `plot_DistributionRuns()` and `summaryRuns()` (closes #41).
* Fixed threshold check in `tableRuns()`.
* Fixed bug in `tableRuns()` filtering logic.
* Fixed CHR casting to string in `snpInsideRunsCpp()`.
* `tableRuns()` now skips gracefully when no data are available for a subset.
* `plot_InbreedingChr()` and `plot_DistributionRuns()` updated for compatibility with ggplot2 `guides(fill = "none")`.

# detectRUNS 0.9.6

## Bug fixes

* Removed compilation warnings in C++ code.
* Removed outdated vignette files.

# detectRUNS 0.9.5

## Minor changes

* Updated references in `DESCRIPTION`.

# detectRUNS 0.9.4

## Bug fixes

* Functions no longer write files to the current working directory (closes security concern).

## Minor changes

* Removed examples from non-exported functions to clean up documentation.
* Updated package description.
* Added CI badges and coverage integration.

# detectRUNS 0.9.3

## Major changes

* First submission to CRAN.

# detectRUNS 0.9.2

## Major changes

* Split main entry point into two dedicated functions: `slidingRUNS.run()` and `consecutiveRUNS.run()`.
* `snpInsideRunsCpp()` implemented in C++ for performance.
* New plotting function: `plot_InbreedingChr()`.
* Added pre-calculated example runs data for use in examples and vignette.
* Renamed variables for consistency across functions and plots.

## Bug fixes

* Fixed bug in `snpInsideRuns()`.
* Fixed examples and suppressed spurious NOTEs during tests.

# detectRUNS 0.6.6

## Major changes

* `consecutiveRuns()` reimplemented in C++ (`consecutiveRunsCpp()`).
* `snpInRun()` and related sliding window functions ported to C++.
* New `plot_DistributionRuns()` function for run length distribution by size class.

## Bug fixes

* Fixed bug in `snpInsideRunsCpp()`.
* Fixed `find_Chromosome_length()` and `reorderDF()`.

# detectRUNS 0.6.2

## Major changes

* `consecutiveRuns()` fully reimplemented in C++ (`consecutiveRunsCpp()`).
* Code refactoring across sliding window and consecutive runs functions.

## Bug fixes

* Fixed `createRunDf()` bugs.
* Fixed `slidingWindowCpp()`.
* Fixed opposite/missing genotype calculation in C++.

# detectRUNS 0.6.1

## Major changes

* New example data files (`.ped`/`.map`).
* Updated and fixed `createRUNdf()` with tests.
* Fixed `reorderDF()` using `dcast()` for long-to-wide conversion.

## Bug fixes

* Fixed `snpRunCpp()`.
* Fixed `homoZygotTest()` and `heteroZygotTest()`.

# detectRUNS 0.5.0

## Major changes

* Added checks for maximum opposite genotypes and maximum missing genotypes within runs (`maxOppRun`, `maxMissRun`).
* New functions: `summaryRuns()`, `Froh_inbreeding()`, `Froh_inbreedingClass()`, `tableRuns()`.
* `consecutiveRuns()` rewritten with improved logic.
* `reorderDF()` added to sort results by chromosome.

## Bug fixes

* Fixed `Stats.R` edge cases with empty subsets.
* Fixed `SnpInRun()` indexing bug in sliding window overlap.

# detectRUNS 0.4.5

## Major changes

* Implemented `consecutiveRuns()` as a second ROH detection method.
* Added Manhattan plot for ROH.
* Fixed bugs in `consecutiveRuns()` start/end position tracking and NA handling.
* Improved documentation throughout.

# detectRUNS 0.2.5

## New features

* New `readPOPCpp()` function to extract population/ID information from ped files.
* New `readFromPlink()` function.
* Improved handling of missing genotype values.

# detectRUNS 0.2.4

## Bug fixes

* Fixed small bug in `plotSNPinside()`.
* Removed `Rcpp::function` import causing errors and warnings.
* Replaced `data.table` ped scanning with line-by-line reading for robustness.

# detectRUNS 0.2.3

## New features

* Added `pedConvert()` function for ped format conversion.
* Input file validation: missing allele pairs in ped now throw an error.

## Bug fixes

* Fixed bugs in `RUNS.run()`.
* Fixed plot functions.

# detectRUNS 0.2.2

## Experimental

* Explored loading genotypes using `bigmemory` (not retained in later versions).

# detectRUNS 0.2.1

## Major changes

* Core functions ported to C++: `genoConvert()`, `homoZygotTest()`, `heteroZygotTest()`, `slidingWindow()`.
* Switched to using C++ functions in main `RUNS.run()` pipeline.
* `prints` replaced with `messages` throughout.
* Renamed `window` parameter to `windowSize` for clarity.
* Added tests and examples for core functions.

# detectRUNS 0.2.0

## Major changes

* First C++ function added (`genoConvertCpp()`).
* Added `snpInRunCpp()` with full row testing.
* New plotting functions for performance visualization.
* Added `detected.ROHet.csv` in `inst/extdata`.
* Package made compatible with `devtools::check()`.

# detectRUNS 0.1

## Major changes

* Initial development version.
