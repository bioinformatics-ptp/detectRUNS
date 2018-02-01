
## test environments
* Local Debian 8, R 3.3.1
* Local Windows 7, R 3.3.2
* ubuntu 14.04 (on travis-ci), R 3.3.3
* ubuntu 14.04 (on travis-ci), R version 3.4.2
* ubuntu 14.04 (on travis-ci), R Under development (unstable) (2018-02-01 r74190)

## R CMD check results
0 errors | 0 warnings | 0 notes

* 1 note under devtools::release():

  Possibly mis-spelled words in DESCRIPTION:
    Heterozygosity (3:48)
    Homozygosity (3:23)
    heterozygosity (13:55)
    homozygosity (13:35)

  `Heterozygosity` and `Homozygosity` are not mis-spelled words.

## Downstream dependencies
There are currently no downstream dependencies for this package
