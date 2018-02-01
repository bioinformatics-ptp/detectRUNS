## ---- eval=TRUE, include=TRUE, echo=FALSE--------------------------------
#just to be sure
library("detectRUNS")

## ------------------------------------------------------------------------
genotypeFilePath <- system.file(
  "extdata", "Kijas2016_Sheep_subset.ped", package="detectRUNS")
mapFilePath <- system.file(
  "extdata", "Kijas2016_Sheep_subset.map", package="detectRUNS")

## ----results='hide',message=FALSE, cache=FALSE, warning=FALSE------------
slidingRuns <- slidingRUNS.run(
  genotypeFile = genotypeFilePath, 
  mapFile = mapFilePath, 
  windowSize = 15, 
  threshold = 0.05,
  minSNP = 20, 
  ROHet = FALSE, 
  maxOppWindow = 1, 
  maxMissWindow = 1,
  maxGap = 10^6, 
  minLengthBps = 250000, 
  minDensity = 1/10^3, # SNP/kbps
  maxOppRun = NULL,
  maxMissRun = NULL
) 

## ----results='hide',message=FALSE, cache=FALSE, warning=FALSE------------
consecutiveRuns <- consecutiveRUNS.run(
  genotypeFile =genotypeFilePath,
  mapFile = mapFilePath,
  minSNP = 20,
  ROHet = FALSE,
  maxGap = 10^6,
  minLengthBps = 250000,
  maxOppRun = 1,
  maxMissRun = 1
)

## ----results='hide',message=FALSE, cache=FALSE, warning=FALSE------------
slidingRuns_het <- slidingRUNS.run(
  genotypeFile = genotypeFilePath, 
  mapFile = mapFilePath, 
  windowSize = 10, 
  threshold = 0.05,
  minSNP = 10, 
  ROHet = TRUE, 
  maxOppWindow = 2, 
  maxMissWindow = 1,
  maxGap = 10^6, 
  minLengthBps = 10000, 
  minDensity = 1/10^6, # SNP/kbps
  maxOppRun = NULL,
  maxMissRun = NULL
) 

## ----results='hide',message=FALSE, cache=FALSE---------------------------
consecutiveRuns_het <- consecutiveRUNS.run(
  genotypeFile =genotypeFilePath,
  mapFile = mapFilePath,
  minSNP = 10,
  ROHet = TRUE,
  maxGap = 10^6,
  minLengthBps = 10000,
  maxOppRun = 2,
  maxMissRun = 1
)

## ---- results='hide',message=FALSE, cache=FALSE--------------------------
summaryList <- summaryRuns(
  runs = slidingRuns, mapFile = mapFilePath, genotypeFile = genotypeFilePath, 
  Class = 6, snpInRuns = TRUE)

## ------------------------------------------------------------------------
summaryList$summary_ROH_count

## ------------------------------------------------------------------------
summaryList$summary_ROH_mean_chr

## ------------------------------------------------------------------------
head(summaryList$SNPinRun)

## ---- fig.show='hold', fig.cap="Plot runs per individual"----------------
plot_Runs(runs = slidingRuns)

## ---- fig.show='hold', fig.cap="Plot runs per individual", message=FALSE, cache=FALSE, warning=FALSE, results='hide'----
plot_StackedRuns(runs = slidingRuns)

## ---- fig.show='hold', message=FALSE, cache=FALSE, warning=FALSE, results='hide',fig.width=6,fig.height=4----
plot_SnpsInRuns(
  runs = slidingRuns[slidingRuns$chrom==2,], genotypeFile = genotypeFilePath, 
  mapFile = mapFilePath)

## ---- fig.show='hold', message=FALSE, cache=FALSE, warning=FALSE, results='hide',fig.width=6,fig.height=4----
plot_SnpsInRuns(
  runs = slidingRuns[slidingRuns$chrom==24,], genotypeFile = genotypeFilePath, 
  mapFile = mapFilePath)

## ----message=FALSE, cache=FALSE, results='hide'--------------------------
topRuns <- tableRuns(
  runs =  slidingRuns, genotypeFile = genotypeFilePath, mapFile = mapFilePath, 
  threshold = 0.7)

## ----echo=FALSE----------------------------------------------------------
print(topRuns)

## ---- message=FALSE, cache=FALSE, results='hide', fig.width=6,fig.height=4----
plot_manhattanRuns(
  runs = slidingRuns[slidingRuns$group=="Jacobs",], 
  genotypeFile = genotypeFilePath, mapFile = mapFilePath)

## ---- echo=FALSE, message=FALSE, cache=FALSE-----------------------------
head(
  Froh_inbreeding(runs = slidingRuns,mapFile = mapFilePath,genome_wide = TRUE))

## ---- echo=FALSE, message=FALSE, cache=FALSE, fig.width=5, fig.height=4----
plot_InbreedingChr(
  runs = slidingRuns, mapFile = mapFilePath, style = "FrohBoxPlot")

## ---- message=FALSE------------------------------------------------------
savedRunFile <- system.file(
  "extdata", "Kijas2016_Sheep_subset.sliding.csv", package="detectRUNS")
runs <- readExternalRuns(inputFile = savedRunFile, program = "detectRUNS")
head(runs)

