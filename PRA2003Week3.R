# PRA2003 Week 3 Deliverable
## Code finding the average count of bacteria in different experiments
## The code reads an output file, counts how many bacteria (specific ID) appear
## It then calculates the mean count and uncertainty for each bacteria ID

# How to run this script
## Open file in VS code or R studio
## You will be asked to type the filename (Ex. ouputSet0.txt or ouputSet1.txt)
## Program finds the presence of the bacteria in each block
## It then calculates its mean and standard error
## If a bacteria is not in a block, the block is not counted into the average
## Each bacteria ID (ex. 211), has a mutant (ex. -211)


#Main function
analyze_bacteria <- function(filename = "output-Set1.txt",
                   target_ids = c(211, -211, 321, -321, 2212, -2212, 3122, -3122, 3312, -3312, 3334, -3334)) {

  id_names <- as.character(target_ids)
  #target ids converted into strings(words) instead of numbers

  #Running total vectors for each bacteria id
  #Start at zero and then increases as it goes through file
  # setNames makes id names as labels

  blocks_present <- setNames(integer(length(target_ids)), id_names)
    #blocks_present = N --> counter for how many blocks the bacteria is present in

  total <- setNames(double(length(target_ids)), id_names)
  # total = running sum of counts for each id 
  # Only in blocks where value actually appeared (if id not present, block skipped)

  total_squared <- setNames(double(length(target_ids)), id_names)
  # used for variance calculation


  #Protection making sure the file actually exists
  if (!file.exists(filename)) {
    stop(paste("Error: could not find the file:", filename))
  }


  # scan() reads entire file and turns it into one long vector as values
  values <- scan(filename, what = numeric(), quiet = TRUE)

  # protection in case file is completely empty (length values = 0)
  # avoids crashing in case of table full of zeros
  if (length(values) == 0) {
    return(data.frame(
      particle_id = as.integer(target_ids),
      blocks_present = 0,
      mean_count_per_block = 0,
      uncertainty = 0
    ))
  }

  # i is index tracking where we are in values
  i <- 1
  while (i <= length(values)) { #loop runs once per block, continues until end of values
    # Header takes up two values, if not, then stop (protection)
    if (i + 1 > length(values)) {
      break
    }

    n_particles <- as.integer(values[i + 1]) #second number in header
      #how many particle lines in block 
    i <- i + 2 #move pointer i past 

    #resets every block
    current_block_counts <- setNames(integer(length(target_ids)), id_names)
      
      # Loop through exactly n. bacteria lines worth of data
    for (j in seq_len(n_particles)) { 
      # particle line has 4 values, if not then stop (protection)
      if (i + 3 > length(values)) {
        break
      }

    #focus only on 4th value (particle_id)
      particle_id <- as.integer(values[i + 3])
        #only count particle if its 12 id numbers, other is ignored
      if (particle_id %in% target_ids) {
        current_block_counts[as.character(particle_id)] <- #numeric id into text
          current_block_counts[as.character(particle_id)] + 1 #increase block count by 1
      }
      i <- i + 4 #move pointer by 4 to next bacteria line
    }
    
    #end of blocks so include results in total
    # go for each bacterial id
    for (id in id_names) {
      if (current_block_counts[id] > 0) { #if id didn't appear then ignore in count
        total[id] <- total[id] + current_block_counts[id]
          #Add block count for id to running sum
        total_squared[id] <- total_squared[id] + current_block_counts[id]^2
          #For variance squared
        blocks_present[id] <- blocks_present[id] + 1 #if id present add to block counter
      }
    }
  } #end of loop 

  mean_val <- ifelse(blocks_present > 0, total / blocks_present, 0)
  # mean = total / N, for each id
  # ifelse() protects against dividing by zero

  variance <- ifelse(
    # variance = (sum of squares - N * mean^2) / (N - 1)
    # requires at least two blocks
    # bacteria count fluctuation from one block to next
    blocks_present > 1,
    (total_squared - blocks_present * mean_val^2) / (blocks_present - 1),
    0
  )
  variance[variance < 0] <- 0 #force negative values to zero (protection)

    # stadard error (uncertainty) = std_dev / sqrt(N)
    # confidence of mean being close to average 
  error <- ifelse(blocks_present > 0, sqrt(variance) / sqrt(blocks_present), 0)

  # combine everything into one table
  results <- data.frame(
    particle_id = as.integer(target_ids),
    blocks_present = as.integer(blocks_present),
    mean_count_per_block = round(mean_val, 8),
    uncertainty = round(error, 8)
  )

  return(results)
}

# Section allowing for user input
run_analysis <- function() {
  filename <- readline(prompt = "Enter file name (e.g. output-Set0.txt or output-Set1.txt): ")
  if (filename == "") { #in case of not writing anything, inputs "output-Set1.txt"
    filename <- "output-Set1.txt"
  }
  
  cat("Please wait a few seconds or minutes...\n")
  cat("Running analysis on:", filename, "\n")
  print(analyze_bacteria(filename), digits = 10) #runs analysis
}

run_analysis()

#End of code with average count of bacteria id's per block
# and their uncertainties