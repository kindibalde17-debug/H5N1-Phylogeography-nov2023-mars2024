install.packages("remotes")
remotes::install_github("sdellicour/seraphim/windows", upgrade = "never")
install.packages("diagram")
library(seraphim)
library(diagram)
library(sf)
ls("package:seraphim")

localTreesDirectory = "Extracted_trees_H5N1"

allTrees = scan(file = "Binome_8b.trees",
                what = "", sep = "\n", quiet = TRUE)

burnIn = 0
randomSampling = FALSE
nberOfTreesToSample = 100
mostRecentSamplingDatum = 2024.1693989071039
coordinateAttributeName = "location"
treeExtractions(localTreesDirectory,
                allTrees,
                burnIn,
                randomSampling,
                nberOfTreesToSample,
                mostRecentSamplingDatum,
                coordinateAttributeName)

nberOfExtractionFiles = 100
timeSlices = 100
onlyTipBranches = FALSE
showingPlots = FALSE
outputName = "H5N1"
nberOfCores = 1
slidingWindow = 1
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices,
                   onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)