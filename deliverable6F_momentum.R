## PRA2003 - momentum of our bacteria

## what this script does:
# 1. goes through the 10 subsample files (output-Set1.txt ...)
# 2. keeps only the bacterias of interest (12 IDs)
# 3. calculates the momentum of each of those bacteria: p = sqrt(px^2 + py^2 + pz^2)
# 4. saves ID + p in momenta-Set1.csv ... momenta-Set10.csv
#the python script then reads these csv files to plot


# the IDs we want: positive = wild type (WT), negative = mutant of the same bacterium
#   211 E. coli, 321 B. subtilis, 2212 P. aeruginosa,
#   3122 S. pneumoniae, 3312 M. tuberculosis, 3334 Salmonella
ids <- c(211, -211, 321, -321, 2212, -2212, 3122, -3122, 3312, -3312, 3334, -3334)

# read files at chunks of 2 million lineas at a time
chunkSize <- 2e6

# one loop round = one subsample file
for (i in 1:10) {
  # paste0 glues text together without spaces (avoid possible error)
  inFile  <- paste0("output-Set", i, ".txt")
  outFile <- paste0("momenta-Set", i, ".csv")

  # stop with a clear message if a file is missing
  # make sure that directory is correct 
  if (!file.exists(inFile)) stop("can't find ", inFile)
  cat("reading", inFile, "\n")   

  # make a new empty csv with only the column names (ID, p)
  # every chunk gets added to bottom of this file
  # overwrites any previous data 
  write.csv(data.frame(ID = integer(0), p = numeric(0)), outFile, row.names = FALSE)

  # open the file as connection (R remembers where it stopped reading)
  con <- file(inFile, "r")
  readLines(con, n = 1)   # he first line is header

  repeat {
    # next chunk of lines (as plain text for now)
    lines <- readLines(con, n = chunkSize)
    if (length(lines) == 0) break   # nothing left -> end of the file -> stop the loop

    # split every line on the spaces into 4 columns: px, py, pz, ID
    # read as textfirst, because not every line is a bacterium:
    # fill = TRUE: line <4 parts --> empty columns 
    # flush = TRUE: line > 4 parts just ignores extra bits
    cols <- scan(text = lines, what = list("", "", "", ""), fill = TRUE, flush = TRUE, quiet = TRUE)

    # 4th column -> headerline -->  NA (not numerical)
    # (suppressWarnings hides NA)
    id <- suppressWarnings(as.numeric(cols[[4]]))

    # TRUE for lines with one of 12 IDs , or else FALSE
    keep <- id %in% ids

    # only keep momentum components
    px <- as.numeric(cols[[1]][keep])
    py <- as.numeric(cols[[2]][keep])
    pz <- as.numeric(cols[[3]][keep])

    # total momentum
    p <- sqrt(px^2 + py^2 + pz^2)


    # reading table properly 
    write.table(data.frame(ID = id[keep], p = round(p, 6)), outFile, sep = ",",
                row.names = FALSE, col.names = FALSE, append = TRUE)
  }
  close(con)   # close the file before moving to next 
}


# count how many of each ID ended up in every csv (compare to week 4)
counts <- t(sapply(1:10, function(i) {
  m <- read.csv(paste0("momenta-Set", i, ".csv"))
  table(factor(m$ID, levels = ids))
}))
# sapply put every set in a column, t() flips it so every row is a set (easier to read)
rownames(counts) <- paste0("Set", 1:10)
print(counts)

cat("\ntotal over all 10 sets:\n")
print(colSums(counts))   # colSums = add up every column -> total per ID
