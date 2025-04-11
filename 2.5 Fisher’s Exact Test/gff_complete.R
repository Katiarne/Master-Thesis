library(tidyverse)

# filter bed file to only contain relevant colums
gff_complete <- read.delim("final_gff.bed", header = FALSE) %>%
  select(V1:V4)
