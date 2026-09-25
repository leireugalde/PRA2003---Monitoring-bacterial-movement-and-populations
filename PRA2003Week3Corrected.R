#PRA2003 - Monitoring Bacterial movement and populations  
## Week 3 deliverable Corrections 
## The file PRA2003Week3 (submitted one) had errors which caused it to get incorrect results 
## This code will be used for the deliverable of Week 4 and onward

## This code finds the average count of bacteria in different experiments
## The code reads an output file, counts how many bacteria (specific ID) appear
## It then calculates the mean count and uncertainty for each bacteria ID

# How to run this script
## Open file in VS code or R studio
## You will be asked to type the filename (Ex. ouputSet0.txt or ouputSet1.txt)
## Program finds the presence of the bacteria in each block
## It then calculates its mean and standard error
## If a bacteria is not in a block, the block is not counted into the average
## Each bacteria ID (ex. 211), has a mutant (ex. -211)
