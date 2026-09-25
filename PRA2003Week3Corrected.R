#PRA2003 - Monitoring Bacterial movement and populations  
## Week 3 deliverable Corrections 
## The file PRA2003Week3 (submitted one) had errors which caused it to get incorrect results 
## This code will be used for the deliverable of Week 4 and onward

## This code finds the average count of bacteria in different experiments
## The code reads an output file, counts how many bacteria (specific ID) appear
## It then calculates the mean count and uncertainty for each bacteria ID

# How to run this script
## Open file in VS code or R studio
## Reads a simulation output file made of repeated blocks:
#   block header : <eventNumber> <nBacteriaInEvent>
#   then <nBacteriaInEvent> data lines : <px> <py> <pz> <strainID>
##For each of the 12 known strain IDs, reports:
#   - how many times it was recorded in total
#   - the average number of occurrences per valid event
        # - "Valid event" = header whose declared bacteria count > 0
#   - the statistical uncertainty on that average
# Every strain's average is divided by the SAME number of valid
# events, so an event where a strain never shows up still counts
# as a zero for that strain rather than being skipped.

#Chosen to look at "output-Set1.txt", so when running, ensure that is downloaded and in a folder where it can run.
#Output-Set1.txt will be used as default 
#Bottom will allow user input to calculate values for any output-Setx.txt (x = 0-10)


defaultFile <- "output-Set1.txt" #fallback filename if no other given 
chunkSize <- 1000000 #read max this many linrd
fieldSep <- " " # separator between the fields of a line
headerFields <- 2 # fields in header (2): eventNumber nBacteria
dataFields <- 4 # fields in a data line (4): px py pz bacterialID
inputPrefix <- "output-" # input "output-<name>.txt" gives ...
outputPrefix <- "results-" # ... output "results-<name>.csv"
outputExtension <- ".csv"
averageDigits <- 5 # decimals kept for the average per event
uncertaintySignifDigits <- 3 # significant digits kept for the uncertainty

#ID for all 12 bacteria, anything else ignored
strainID <- c(211, -211, 321, -321, 2212, -2212, 3122, -3122, 3312, -3312, 3334, -3334)


#Names of each strain in ID respective order
strainName <- c(
"E. coli WT",
"E. coli mutant",
"Bacillus subtilis WT",
"Bacillus subtilis mutant",
"Pseudomonas aeruginosa WT",
"Pseudomonas aeruginosa antibiotic-resistant",
"Streptococcus pneumoniae",
"Capsule-deficient S. pneumoniae",
"Mycobacterium tuberculosis",
"Drug-resistant M. tuberculosis",
"Salmonella enterica",
"Salmonella mutant"
)

#number of strains tracked based on IDs listed above
nStrains <- length(strainID)


# functions validating settings and strain list before running
isWholeNumber <- function(x, minimum) {
is.numeric(x) && length(x) == 1 && !is.na(x) && x >= minimum && x == round(x)
}
#true if x is single whole number (minimum value)

isText <- function(x) {
is.character(x) && length(x) == 1 && !is.na(x) && nzchar(x)
}
#true only if x is single, non-empty piece of text 

if (!isWholeNumber(chunkSize, 1)) {
stop("Setting chunkSize must be a whole number of at least 1, not: ", chunkSize)
}

if (!isWholeNumber(headerFields, 1) || !isWholeNumber(dataFields, 1) ||
headerFields == dataFields) {
stop("Settings headerFields and dataFields must be whole numbers and different from ",
"each other .")
}

if (!isText(fieldSep) || nchar(fieldSep) != 1) {
stop("Setting fieldSep must be a single character.")
}

if (!isWholeNumber(averageDigits, 0) ||
!isWholeNumber(uncertaintySignifDigits, 1)) {
stop("Settings averageDigits and uncertaintySignifDigits must be whole numbers.")
}

if (!isText(defaultFile) || !isText(outputPrefix) || !isText(outputExtension)) {
stop("Setting defaultFile, outputPrefix and outputExtension must be text.")
}

if (length(strainName) != nStrains) {
stop("strainID has ", nStrains, " values but strainName has ", length(strainName), ".")
}
#ensure that strainID and strainName line up, if they dont match stop 


if (anyNA(strainID) || anyDuplicated(strainID) > 0) {
stop("strainID must not contain missing or duplicate values.")
}
#guard from missing or duplicate values


#lineRegex(n) builds a pattern that matches a line with
#exactly n fields separated by fieldSep
lineRegex <- function(nFields) {
notSep <- paste0("[^", fieldSep, "]")
paste0("^", notSep, "+(", fieldSep, notSep, "+){", nFields - 1, "}$")
}

headerRegex <- lineRegex(headerFields)

dataRegex <- lineRegex(dataFields)

lastFieldRegex <- paste0("^.*", fieldSep)


#Process chunk cleans and tallies on batch of lines at a time
processChunk <- function(chunk) {
chunk <- sub("\r$", "", chunk, perl = TRUE) 
# returns true if strings are not empty, and false if they are
chunk <- chunk[nzchar(chunk)] # skip empty lines

#isHeader used to classify the line in chunk as header line or data line
    # if true, then looks like header, if false it does not
#perl = TRUE, used to tell R to use "per-stlye regen synthax"
    #(used to match patterns)
isHeader <- grepl(headerRegex, chunk, perl = TRUE)
isData <- !isHeader & grepl(dataRegex, chunk, perl = TRUE)
#line counts as data if not claimed as header 


#pulls out last field of each line
  # for header = bacteria count 
  # for data lines = bacterial ID
  # rest of columns are ignored 
lastField <- sub(lastFieldRegex, "", chunk, perl = TRUE)


#Look only at the header lines now. A header is considered malformed if its
#count is missing, negative, or not a whole number.

nBacteriaInEvent <- suppressWarnings(as.numeric(lastField[isHeader]))
#converts text values into numbers


badHeader <- is.na(nBacteriaInEvent) | nBacteriaInEvent < 0 |
#flagged as a bad header if it's missing, negative, or not a whole number

nBacteriaInEvent != round(nBacteriaInEvent)
nBacteriaInEvent <- nBacteriaInEvent[!badHeader]
#removes the malformed headers


#Now look at the data lines and pull out the bacterial ID. A non-numeric
#value here means the line is malformed

#lastField[isData] keeps only the data-line values, as.numeric converts them
#to numbers, and suppressWarnings hides the coercion warnings that would
#otherwise show up for anything that fails to convert
ids <- suppressWarnings(as.numeric(lastField[isData]))

#badID marks the ones that cant be used and removes them
badId <- is.na(ids)
ids <- ids[!badId]
row <- match(ids, strainID) # position in strainID, NA if not one of the strains

#known marks which IDs actually matched one of our 12 tracked strains
known <- !is.na(row)


#list used to show data in a structured way
list(
nEvents = sum(nBacteriaInEvent > 0), #counts headers with numbers more than 0
counts = tabulate(row[known], nbins = nStrains), #counts how many tracked strains showed up
nOther = sum(!known), #counts IDs not recognized
nMalformed = sum(!isHeader & !isData) + sum(badHeader) + sum(badId), #totals the lines that are not valid
nDeclared = sum(nBacteriaInEvent),
nData = length(ids) #records valid data rows after cleaning
)
}




#read in what was passed on command line after script name
args <- commandArgs(trailingOnly = TRUE)


#if lenght is more than 1 then it stops
if (length(args) > 1) {
stop(
"Give at most one argument (the input file), but got ", length(args), ": ",
paste(args, collapse = " ")
)
}

#Use argument as file path if one was given, otherwise back to default
filepath <- if (length(args) == 1) args[1] else defaultFile

if (!isText(filepath)) {
stop("The input file name is empty.")
}

#ensures that the file used exists, if not then states file not found
if (!file.exists(filepath)) {
stop("File not found: ", filepath, " - put it in the same folder as this script.")
}

#protection incase file is in folder and not file, allows correction
if (dir.exists(filepath)) {
stop(filepath, " is a folder, not a file.")
}

#mode=4 is read permision, if its not zero, then no permition to read file
if (file.access(filepath, mode = 4) != 0) {
stop("No permission to read: ", filepath)
}

#makes sure file size is not empty (if empty then show error "The file is empty")
if (is.na(file.size(filepath)) || file.size(filepath) == 0) {
stop("The file is empty: ", filepath)
}



totalCount <- rep(0, nStrains) # creates count per strain - starting at 0
nEvents <- 0 # start counter for events with at least one bacteria
nOther <- 0 # start counter for bacteria with an ID outside the strains
nDeclared <- 0 # start counter for total articles announced by the header lines
nData <- 0 # counter for valid data lines found
nMalformed <- 0 # counter for lines that are not a valid header or data line
linesRead <- 0 # lines read so far (used in error messages)

con <- file(filepath, "r") #open input file in read-only  and store it as con

tryCatch({
repeat {
chunk <- readLines(con, n = chunkSize)

#stop once there is nothing left 
if (length(chunk) == 0) {
  break
}

result <- processChunk(chunk)

totalCount <- totalCount + result$counts
nEvents    <- nEvents + result$nEvents
nOther     <- nOther + result$nOther
nMalformed <- nMalformed + result$nMalformed
nDeclared  <- nDeclared + result$nDeclared
nData      <- nData + result$nData
linesRead  <- linesRead + length(chunk)

}

#defines what to do if error occurs while reading
}, error = function(e) {
stop(
"Reading ", filepath, " failed after about ", linesRead, " lines: ",
conditionMessage(e),
call. = FALSE
)
}, finally = close(con))

if (nEvents == 0) {
stop("No events found in file - check the file path/format.")
}


#if nothing found and sum =0 then no strains found
if (sum(totalCount) == 0) {
stop("None of the ", nStrains, " strain IDs occur in the file - check the file format.")
}

if (nMalformed > 0) {
warning(
nMalformed, " lines were not a valid header (", headerFields, " fields) or ",
"data line (", dataFields, " fields with a numeric ID) and were ignored"
)
}

#check if number of bacteria declared by headers matches number of valid data lines
if (nDeclared != nData) {
warning("Headers announce ", nDeclared, " bacteria but ", nData, " were found")
}


#computing the statilistcs
#average per event = total count / number of events
#uncertainty = sqrt(total count) / number of events

avgPerEvent <- totalCount / nEvents
uncertainty <- sqrt(totalCount) / nEvents

results <- data.frame(
ID = strainID,
Strain = strainName,
TotalCount = totalCount,
AveragePerEvent = round(avgPerEvent, averageDigits),
Uncertainty = signif(
uncertainty,
uncertaintySignifDigits
) # significant digits: the smallest values are tiny
)

cat("Total events processed:", nEvents, "\n")
cat("Bacteria with an ID outside the", nStrains, "strains (ignored):", nOther, "\n")
cat("Lines that were not valid (ignored):", nMalformed, "\n\n")
print(results)
