## Execture main.R and save output
# dev must have two environment vars in their .Renviron :
# FISHMAP_UPDATE_OUTPUTS=TRUE
# FISHMAP_OUTPUT_DIR=path/to/output/dir

# load lib to record execution time
library(tictoc)

# Set seed to make FishMap run reproducible
set.seed(29510)

# Setup output dir
output_dir <-  Sys.getenv("FISHMAP_OUTPUT_DIR")
dir.create(output_dir, recursive = TRUE)

# checking writing status
if (Sys.getenv("FISHMAP_UPDATE_OUTPUTS") == "TRUE") {
  message(sprintf("Saving output files in %s", output_dir))
} else {message("Running main.R without saving output")}

# run main.R
tic("run main.R")
source(here::here("dev/r_scripts_orig/main.R"))
toc()
