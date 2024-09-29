intersect_data <- read.table("intersect_overlap_output")

unique(intersect_data$V9) # to get the unique data points

# list of uniqe genes
attributions <- intersect_data$V9
unique_genes <- data.frame(table(attributions))
colnames(unique_genes)[colnames(unique_genes) == "Freq"] <- "SVs"
unique_genes

# count nr. of bases in Reference sequene and alternative sequence
library(tidyverse)

intersect_data <- intersect_data %>%
  mutate(reference_count = nchar(V13, type = "chars"),
         alternative_count = nchar(V14, type = "chars"))

# if reference have highest nr. of bases -> Deletions
# if alterantive have highest nr. of bases -> Insertion
intersect_data <- intersect_data %>%
  mutate(variant_type = ifelse(reference_count > alternative_count, 
                               "deletion", 
                               "insertion"))

# new tabel with only deletions
deletion_data <- intersect_data %>%
  filter(variant_type == "deletion")
