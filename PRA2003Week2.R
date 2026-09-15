#Output0
#Formula for momentum
calculate_momentum <- function(px,py,pz) { #Function for value input
    p <- sqrt(px^2 + py^2 + pz^2)
    return(p) #Return calculated value with formula
}

#Load data 
data <- read.table("output-Set0.txt", #Reads data file into table
        skip = 1, #Skip first line as it is a header
        header = FALSE, #File has no column names yet
        col.names = c("px", "py", "pz", "particle_id")) #Name columns
    return(data)

#Loop through every row to calculate their momentums
momentum <- c() #Empty vector to enter values
for (i in 1:nrow(data)) { #Run line by line
    #Take specific value ($) from specific row [i]
    specif_px <- data$px[i]
    specif_py <- data$py[i]
    specif_pz <- data$pz[i]
    #Calculate momentum of specific [i] values
    specif_p <- calculate_momentum(specif_px, specif_py, specif_pz)
    #Append calculated momentum onto end of table
    momentum <- c(momentum, specif_p)
}
#Attach momentum resutls as new column in table aligning properly
data$momentum <- momentum

#Print final result
print("Output Set 0")
print(data)
print("The table above shows the magnitude of the momentum (movement) of each bacteria")
print("The bacteria with particle_id 211 in row 12 has the largest magnitude")