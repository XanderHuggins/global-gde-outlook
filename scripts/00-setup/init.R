### ---------------------\\ 
# Script objective:
# Import custom functions, set seed, and set working temp directory for spatial processes 
### ---------------------\\ 

# import custom functions
invisible(sapply(paste0(here("scripts/0-functions"), "/", list.files(here("scripts/0-functions"))), source)) 

# set temporary terra directory to external disk with storage availability
terraOptions(tempdir = "D://Geodatabase/Rtemp")
tmpFiles(current=TRUE, remove=TRUE) 

set.seed(1234)