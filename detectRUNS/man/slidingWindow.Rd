% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funktionen.R
\name{slidingWindow}
\alias{slidingWindow}
\title{Function to slide a window over a vector (individual's genotypes)}
\usage{
slidingWindow(
  data,
  gaps,
  windowSize,
  step,
  maxGap,
  ROHet = TRUE,
  maxOppositeGenotype = 1,
  maxMiss = 1
)
}
\arguments{
\item{data}{vector of 0/1/2 genotypes}

\item{gaps}{vector of differences between consecutive positions (gaps) in bps}

\item{windowSize}{size of window (n. of SNP)}

\item{step}{by which (how many SNP) is the window slid}

\item{maxGap}{max distance between consecutive SNP in a window to be still considered a potential run}

\item{ROHet}{shall we detect ROHet or ROHom?}

\item{maxOppositeGenotype}{max n. of homozygous/heterozygous SNP}

\item{maxMiss}{max. n. of missing SNP}
}
\value{
vector of TRUE/FALSE (whether a window is homozygous or NOT)
}
\description{
This is a core function. The functions to detect RUNS are slid over the genome
}
