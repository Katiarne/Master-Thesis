# Load intesect file created in bedtools.sh
intersect_data <- read.table("intersect_overlap_output")


# make sure no lncRNA
cds_only <- intersect_data[grep("ID=CDS", intersect_data$V9), ]

unique(intersect_data$V9) # to get the unique data points

# telle antall baser fra "reffrence" og "alternative"
library(tidyverse)

intersect_data <- intersect_data %>%
  mutate(reference_count = nchar(V13, type = "chars"),
         alternative_count = nchar(V14, type = "chars"))

# decide if deletion or insertion
intersect_data <- intersect_data %>%
  mutate(variant_type = ifelse(reference_count > alternative_count, 
                               "deletion", 
                               "insertion"))

# new datafraim with only deletions
deletion_data <- intersect_data %>%
  filter(variant_type == "deletion")

length(unique(deletion_data$V9))



