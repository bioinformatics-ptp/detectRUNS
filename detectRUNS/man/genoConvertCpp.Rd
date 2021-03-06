% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{genoConvertCpp}
\alias{genoConvertCpp}
\title{Convert 0/1/2 genotypes to 0/1}
\usage{
genoConvertCpp(genotype)
}
\arguments{
\item{genotype}{vector of 0/1/2 genotypes}
}
\value{
converted vector of genotypes (0/1)
}
\description{
This is a utility function, that convert 0/1/2 genotypes (AA/AB/BB) into 0/1
(either homozygous/heterozygous)
}
