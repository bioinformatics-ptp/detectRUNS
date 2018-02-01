
Detect Runs of Homozygosity and Runs of Heterozygosity in diploid genomes
=========================================================================

[![Travis-CI Build Status](https://travis-ci.org/bioinformatics-ptp/detectRUNS.svg?branch=master)](https://travis-ci.org/bioinformatics-ptp/detectRUNS)

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
