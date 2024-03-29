% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_DistributionRuns}
\alias{plot_DistributionRuns}
\title{Plot Distribution of runs}
\usage{
plot_DistributionRuns(
  runs,
  mapFile,
  groupSplit = TRUE,
  style = c("MeanClass", "MeanChr", "RunsPCT", "RunsPCT_Chr", "All"),
  savePlots = FALSE,
  outputName = NULL,
  plotTitle = NULL,
  Class = 2
)
}
\arguments{
\item{runs}{R object (dataframe) with results on detected runs}

\item{mapFile}{Plink map file (for SNP position)}

\item{groupSplit}{plots split by group (defaults to TRUE)}

\item{style}{type of plot: MeanClass, MeanChr, RunsPCT, RunsPCT_Chr, All (all plots)}

\item{savePlots}{should plots be saved out to files or plotted in the graphical terminal (default)?}

\item{outputName}{title prefix (the base name of graph, if savePlots is TRUE)#'}

\item{plotTitle}{title in plot (default NULL)}

\item{Class}{group of length (in Mbps) by class (default: 0-2, 2-4, 4-8, 8-16, >16)}
}
\value{
plot Distribution Runs
}
\description{
This function the distribution of runs per group. The average run length per size-class,
the average run length per chromosome (and group), the percent distribution of runs
per size-class and group, and the proportion of runs per chromosome are plotted.
With \code{style="All"} all three plots are produced.
}
\examples{
# getting map and ped paths
genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")

# calculating runs of Homozygosity
\dontrun{
# skipping runs calculation
runs <- slidingRUNS.run(genotypeFile, mapFile, windowSize = 15, threshold = 0.1,  minSNP = 15,
ROHet = FALSE,  maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 100000,  minDensity = 1/10000)
}
# loading pre-calculated data
runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.sliding.csv", package="detectRUNS")
runsDF <- readExternalRuns(inputFile = runsFile, program = 'detectRUNS')

plot_InbreedingChr(runs = runsDF, mapFile = mapFile, style='All')

}
