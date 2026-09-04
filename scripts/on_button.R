# Set working directory to the location of the project
library(here)

# Source setup scripts, including: wd args, plotting themes, custom functions, etc.
invisible(sapply(list.files(here("scripts/0-setup"), full.names = T), source))
