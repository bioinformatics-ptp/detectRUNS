# detectRUNS 1.0.0

## Major rewrite — new engine, new API, new output format

### C++ BED engine
Native reader for PLINK binary (BED/BIM/FAM) files, parallelised via OpenMP.
Replaces R-level PED parsing for large datasets (10–50x faster on typical livestock
arrays). Results are validated to be numerically identical to the previous R engine.

### Unified entry point: `scanRUNS()`
Replaces the deprecated `slidingRUNS.run()` and `consecutiveRUNS.run()`.
Auto-detects file format (BED or PED) from the file extension.
Accepts both methods via `method = "sliding"` or `method = "consecutive"`.
Returns an `ROH` S3 object.

### `ROH` S3 class
`scanRUNS()` now returns a rich object carrying the run table, per-individual
summary, chromosome map, sample info, SNP map, scan parameters, and session metadata.
All downstream functions (`summaryRuns`, `Froh_inbreeding`, all plot functions) accept
an `ROH` object directly — no file paths needed after the initial scan.
`print()`, `as.data.frame()`, and `summary()` methods provided.

### `saveROH()` / `loadROH()`
Binary serialisation of scan results. Save once, reload instantly in future sessions
without re-running the scan.

### `as_ROH()`
Construct an `ROH` object from pre-existing results (e.g. loaded via
`readExternalRuns()` or from an earlier session).

### `rohIslands()`
Permutation-based ROH island detection (method from Falchi et al. 2006,
*PLoS Genetics*). For each chromosome, builds a null SNP-in-ROH distribution by
randomly permuting sample identity `n_perm` times and re-running the ROH scan.
SNPs whose observed frequency exceeds the chromosome-specific permutation threshold are
declared ROH islands. The permutation loop is fully parallelised in C++ via OpenMP.
`print()`, `summary()`, and `plot()` S3 methods provided.

### `reportRUNS()`
Comprehensive report generator. Writes a Markdown, HTML, or PDF document from an
`ROH` object. Includes: executive summary, scan parameters, dataset overview,
per-chromosome coverage, run length class distribution, top SNPs and chromosomes in
ROH, common ROH regions, individual outlier flags, inbreeding coefficients (F_ROH),
and (when `islands` is supplied) a ROH island Manhattan plot.
HTML output is fully self-contained: all plots are base64-embedded, no external files
required.

## Dependency removal
Removed `plyr`, `itertools`, `iterators`, and `reshape2` from `Imports`.
All functionality replaced with base R and `data.table` equivalents.
Remaining `Imports`: `ggplot2`, `Rcpp`, `gridExtra`, `data.table`.

## API changes
* `snpInsideRuns()` third argument is now `sample_info` (a `data.frame` with
  `group`/`id` columns) instead of `genotypeFile`.
* `plot_PatternRuns()`: removed unused `mapFile` parameter.

## Bug fixes
* Fixed silent wrong output in `snpInsideRuns()`: column named `"GROUP"` but callers
  expected `"BREED"` (affected `plot_SnpsInRuns` and `plot_manhattanRuns`).
* Fixed `tableRuns()` crash on empty results.
* Fixed `guides(fill=FALSE)` deprecation warnings (changed to `guides(fill="none")`).

## Backward compatibility
`slidingRUNS.run()` and `consecutiveRUNS.run()` remain available as deprecated
wrappers and continue to work unchanged. They will be removed in the next release.
All statistics and plot functions still accept a plain `data.frame` alongside explicit
`mapFile=` / `genotypeFile=` arguments.

---

# detectRUNS 0.9.6

* Last release on CRAN before the 1.0.0 rewrite.

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

* First submission to CRAN.
