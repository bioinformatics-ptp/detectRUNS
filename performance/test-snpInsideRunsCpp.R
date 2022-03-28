
# get files path
mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package="detectRUNS")
genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package="detectRUNS")
runsfile <- system.file("extdata", "Kijas2016_Sheep_subset.sliding.csv", package="detectRUNS")

# read mapfile with custom methods
mappa <- readMapFile(mapFile)

# this is the chromosome I want to test
chrom <- "24"

# subsetting mapChrom (get x random snps from mapfile)
set.seed(42)
mapChrom <- mappa[mappa$CHR==chrom, ]
mapChrom <- mapChrom[sort(sample(nrow(mapChrom), 10)), ]

# loading pre-calculated data
runs <- readExternalRuns(runsfile, program="detectRUNS")
names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")
runsChrom <- runs[runs$CHROMOSOME==chrom, ]

# get snps inside runs
for (i in 1:100000) {
  test <- snpInsideRunsCpp(runsChrom, mapChrom, genotypeFile)  
}
