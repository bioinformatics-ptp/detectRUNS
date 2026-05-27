# detectRUNS

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**detectRUNS** detects **Runs of Homozygosity** (ROHom) and **Runs of
Heterozygosity** (ROHet) in diploid genomes from PLINK-formatted
genotype data. It implements two complementary detection algorithms —
sliding-window and consecutive-SNP — and provides functions to
summarise runs, compute inbreeding coefficients, identify ROH islands
via permutation testing, and produce publication-ready plots.

## Installation

Install the released version from CRAN:

```r
install.packages("detectRUNS")
```

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("bioinformatics-ptp/detectRUNS/detectRUNS")
```

## Quick start

```r
library(detectRUNS)

# Locate the bundled example dataset (Kijas 2016 sheep SNP array)
bed <- system.file("extdata", "Kijas2016_Sheep_subset.bed",
                   package = "detectRUNS")

# Scan for runs of homozygosity — sliding-window method
roh <- scanRUNS(bed, method = "sliding",
                minSNP = 15, minLengthBps = 100000)
print(roh)
#> ROH object  [method: sliding | type: ROHom]
#>   Samples  : 100  (2 groups)
#>   Runs     : 2388  across 2 chromosomes

# Summarise run lengths by class
summaryRuns(roh, Class = 2)

# Genomic inbreeding coefficients (F_ROH)
Froh_inbreeding(roh)

# Manhattan-style plot of run frequency across the genome
plot_manhattanRuns(roh)
```

## Key features

- **Two detection methods**: sliding-window (Purcell et al., 2007) and
  consecutive-SNP (Marras et al., 2015)
- **BED and PED/MAP input**: native PLINK formats; the BED engine uses
  C++ with optional OpenMP parallelism for large datasets
- **ROH islands**: permutation-based detection of genomic regions with
  unusually high ROH frequency (`runsIslands()`)
- **Inbreeding coefficients**: F_ROH by individual and chromosome
  (`Froh_inbreeding()`)
- **Plots**: run distribution, Manhattan frequency plot, per-chromosome
  inbreeding, ROH island manhattan plot
- **Persistent results**: save and reload ROH objects with `saveROH()`
  / `loadROH()` — no need to re-scan large files
- **Automated reports**: generate HTML or Markdown summaries with
  `reportRUNS()`

## Documentation

```r
vignette("detectRUNS.vignette", package = "detectRUNS")
```

Full API reference: <https://github.com/bioinformatics-ptp/detectRUNS>

## Reporting issues

<https://github.com/bioinformatics-ptp/detectRUNS/issues>

## References

Purcell S. et al. (2007). PLINK: a tool set for whole-genome
association and population-based linkage analyses. *American Journal of
Human Genetics*, 81(3), 559–575.
<https://doi.org/10.1086/519795>

Marras G. et al. (2015). Runs of homozygosity reveal genome-wide
autozygosity in Italian sheep breeds. *Animal Genetics*, 46(5),
484–491. <https://doi.org/10.1111/age.12259>

## License

GPL-3. See [LICENSE](https://www.gnu.org/licenses/gpl-3.0) for details.
