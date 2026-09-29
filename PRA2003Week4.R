## PRA2003 - Monitoring Bacterial movement and populations
## Week 4 deliverable
#Runs the Deliverable 3 analysis on each of the 10 datasets and combines
#them into one overall average and standard deviation per strain.


# Code used to calculate mean and uncertainty of each bacteria ID
# Loads strain ID, strain name, settings and the analyzeFile() function
# Used as function/source to not repeat entire code
source("PRA2003Week3Corrected.R")

## which files to combine
#Names of output are 10 datasets (output-Set1.txt to output-Set10.txt)
#match inputPrefix pattern from week 3

cat("Please wait for a few minutes to see results...\n")

setFiles <- paste0(inputPrefix, "Set", 1:10, ".txt")


#protection in case one of files is not present
for (f in setFiles) {
  if (!file.exists(f)) {
    stop("Cant find file: ", f, " - Check its location.")
  }
}


nSets <- length(setFiles)

#total per strain across 10 files, and total number of events across all
# start at 0 and add on after every file
combinedTotalCount <- rep(0, nStrains)
combinedEvents <- 0
combinedOther <- 0
combinedMalformed <- 0

#average per event of every strain in every dataset
#one row per dataset, one column per strain, for the standard deviation
setAverages <- matrix(0, nrow = nSets, ncol = nStrains)

for (i in 1:nSets) {

  cat("Results for", setFiles[i], "...\n")

  result <- analyzeFile(setFiles[i]) #runs week 3 code on one file

  #adding specific files' numbers onto the running totals
  combinedTotalCount <- combinedTotalCount + result$totalCount
  combinedEvents <- combinedEvents + result$nEvents
  combinedOther <- combinedOther + result$nOther
  combinedMalformed <- combinedMalformed + result$nMalformed

  #save files average per event for each strain
  setAverages[i, ] <- result$totalCount / result$nEvents
}

cat("\nFinished reading all files.\n\n")


# Mean: divide summed totals by the summed events
    ## not averaging the 10 separate averages (weights each dataset by its events)
# Standard deviation: spread of the 10 dataset averages for each strain

combined_mean <- combinedTotalCount / combinedEvents
combined_sd <- apply(setAverages, 2, sd)

combinedResults <- data.frame(
  ID = strainID,
  Strain = strainName,
  TotalCount = combinedTotalCount,
  AveragePerEvent = round(combined_mean, averageDigits),
  StandardDeviation = signif(combined_sd, uncertaintySignifDigits)
)

write.csv(combinedResults, paste0(outputPrefix, "Combined", outputExtension), row.names = FALSE)

cat("Combined results across all", nSets, "datasets \n")
cat("total events processed:", combinedEvents, "\n")
cat("IDs outside the", nStrains, "tracked strains (ignored):", combinedOther, "\n")
cat("malformed lines (ignored):", combinedMalformed, "\n\n")
print(combinedResults)
