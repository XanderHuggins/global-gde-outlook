### ---------------------\\ 
# Script objective:
# Set seed, working temp directory for spatial processes 
### ---------------------\\ 

# set temporary terra directory to external disk with storage availability
terraOptions(tempdir = "D://Geodatabase/Rtemp")
tmpFiles(current=TRUE, remove=TRUE) 

set.seed(1234)