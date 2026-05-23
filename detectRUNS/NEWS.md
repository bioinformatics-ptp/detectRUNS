# detectRUNS 1.0.0

## Major new features

* **C++ BED engine**: native reader for PLINK binary (BED/BIM/FAM) files, parallelised
  via OpenMP. Replaces R-level PED parsing for large datasets. Results are validated
  to be byte-identical to PLINK `--homozyg` output.

* **`scanRUNS()`**: new unified entry point replacing the deprecated
  `slidingRUNS.run()` and `consecutiveRUNS.run()`. Auto-detects file format
  (BED or PED) from the file extension. Returns an `ROH` S3 object.

* **`ROH` S3 class**: `scanRUNS()` now returns a rich object that carries the
  run table, per-individual summary, chromosome lengths, sample info, and SNP map.
  All downstream functions (`summaryRuns`, `Froh_inbreeding`, all plot functions)
  accept an `ROH` object directly — no file paths needed after the initial scan.

* **`saveROH()` / `loadROH()`**: binary serialisation of scan results. Save once,
  reload instantly in future sessions without re-running the scan.

* **`as_ROH()`**: build an `ROH` object from pre-existing results (e.g. loaded via
  `readExternalRuns()` or from an earlier session).

## Dependency removal

* Removed `plyr`, `itertools`, `iterators`, and `reshape2` from `Imports`.
  All functionality replaced with base R and `data.table` equivalents.
  Remaining imports: `ggplot2`, `Rcpp`, `gridExtra`, `data.table`.

## API changes

* `snpInsideRuns()` signature changed: third argument is now `sample_info`
  (a `data.frame` with `group`/`id` columns) instead of `genotypeFile`.
  Callers updated accordingly.

* `plot_PatternRuns()`: removed unused `mapFile` parameter.

## Bug fixes

* Fixed silent wrong output in `snpInsideRuns()`: column was named `"GROUP"`
  but callers expected `"BREED"` (affected `plot_SnpsInRuns` and
  `plot_manhattanRuns`).

* Fixed `tableRuns()` crash on empty results (`seq(1, 0)` returns `c(1, 0)`
  not an empty vector).

* Fixed `guides(fill=FALSE)` deprecation warnings across plot functions
  (changed to `guides(fill="none")`).

## New functions

* **`runsAssociation()`**: tests the association between run presence/absence and
  a quantitative phenotype using linear regression. For each run region carried by
  at least `minFreq` fraction of individuals, fits `phenotype ~ run_presence` and
  returns effect sizes with Bonferroni and Benjamini-Hochberg FDR-adjusted p-values.
  Accepts an `ROH` object or a plain data.frame of runs.

## Backward compatibility

* `slidingRUNS.run()` and `consecutiveRUNS.run()` still exist as deprecated
  wrappers and continue to work unchanged.

* All statistics and plot functions still accept a plain `data.frame` alongside
  explicit `mapFile=` / `genotypeFile=` arguments.

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
