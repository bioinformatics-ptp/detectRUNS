
## test environments

Package tested using GitHub workflow [r-lib/actions](https://github.com/r-lib/actions)

* os: macOS-latest,   R: release
* os: windows-latest, R: release
* os: ubuntu-latest,  R: devel
* os: ubuntu-latest,  R: release
* os: ubuntu-latest,  R: oldrel-1
## R CMD check results
0 errors ✔ | 0 warnings ✔ | 1 note ✖

❯ checking installed package size ... NOTE
    installed size is  7.1Mb
    sub-directories of 1Mb or more:
      extdata   2.1Mb
      libs      4.0Mb

Installed package size is greater than 5Mb due to functions compiled with `RCpp`
and some example data used in vignettes

## Downstream dependencies
There are currently no downstream dependencies for this package
