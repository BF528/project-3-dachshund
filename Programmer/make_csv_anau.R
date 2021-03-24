# Combines featureCounts_output files into one csv file

# Import libraries:
library(tidyverse)
library(janitor)

# Get list of applicable file names:
all_files <- list.files(pattern="^featureCounts_output.*txt$")
all_files
num_files <- length(all_files)

# For every file:
for (i in 1:num_files){
  my_file <- all_files[i]
  # Read in table:
  temp_ds <- read.table(my_file, header=TRUE)
  # Grab 1st (should be Geneid) and last (counts in Sample)
  temp_ds <- temp_ds[,c(1, ncol(temp_ds))]
  # Convert to tibble
  temp_tib <- as_tibble(temp_ds)
  # Start combined tibble if first file, otherwise merge with existing combined tibble:
  if (i==1){
    comb_tib <- temp_tib
  }else{
    comb_tib <- full_join(comb_tib, temp_tib, by="Geneid")
  }
}

# Shorten column names:
# WARNING: THIS STEP MAY CAUSE ERRORS IF NOT UPDATED CORRECTLY
# Example full column name:
# X.project.bf528.project_3.samples.SRR1178055__rn4_STAR__sortedByCoord.bam
# names(comb_tib) <- substring(names(comb_tib), first=27, last=45)
# TODO: adjust for our actual files: (can maybe instead remove first row?)
names(comb_tib) <- sub("X.project.bf528.project_3.samples.", "", names(comb_tib))
names(comb_tib) <- sub("__rn4_STAR__sortedByCoord.bam", "", names(comb_tib))

# Write csv:
write.csv(comb_tib, file="featureCounts_combined.csv", row.names=FALSE)
