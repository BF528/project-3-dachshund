# Combines featureCounts_output files into one csv file

# To run from command line:
# module load R/4.0.2
# Rscript make_csv_anau.R

#Install janitor:
install.packages("janitor", repos = "http://cran.us.r-project.org")

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

# Check column names:
names(comb_tib)

# Shorten column names:
# WARNING: THIS STEP MAY CAUSE ERRORS IF NOT UPDATED CORRECTLY
# Example full column name:
# X.projectnb2.bf528.users.dachshund.project_3.samples.STAR_output.SRR1178046Aligned.sortedByCoord.out.bam
names(comb_tib) <- sub("X.projectnb2.bf528.users.dachshund.project_3.samples.STAR_output.", "", names(comb_tib))
names(comb_tib) <- sub("Aligned.sortedByCoord.out.bam", "", names(comb_tib))

# Write csv:
# (Geneid column labeled "Geneid" as in input datasets)
write.csv(comb_tib, file="featureCounts_combined.csv", row.names=FALSE)
