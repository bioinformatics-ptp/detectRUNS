
Detect Runs of Homozygosity and Runs of Heterozygosity in diploid genomes
=========================================================================

<!-- badges: start -->
[![R-CMD-check](https://github.com/bioinformatics-ptp/detectRUNS/workflows/R-CMD-check/badge.svg)](https://github.com/bioinformatics-ptp/detectRUNS/actions)
[![codecov.io](https://codecov.io/github/bioinformatics-ptp/detectRUNS/coverage.svg?branch=master)](https://codecov.io/github/bioinformatics-ptp/detectRUNS?branch=master)
[![CRAN version](http://www.r-pkg.org/badges/version/detectRUNS)](https://cran.r-project.org/package=detectRUNS)
<!-- badges: end -->

This repository contains the source code for the R package `detectRUNS` and related
`performance` tests. Here's the directory content:

```
├── detectRUNS
├── performance
├── README.md
└── TODO
```

`detectRUNS` implements two statistical approaches to runs' detection:
- a sequential approach (consecutive runs, as described in Marras et al. 2015, and implemented in the software package Zanardi: https://github.com/bioinformatics-ptp/Zanardi);
- an approach based on overlapping sliding windows (à la Plink: Purcell et al. 2007)

More info in the Google Doc (https://goo.gl/yXR7iA)
