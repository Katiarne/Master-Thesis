library(tidyverse)

vfc_tissues <- read.delim("/mnt/SCRATCH/kristenl/katinka_master/Ssalv3.1_pangenome_pangenie_all_svs_aquafaang_samples.chr.vcf", 
                              header = FALSE, comment.char = "#", sep = "\t") %>%
  mutate(reference_count = nchar(V4, type = "chars"),
         alternative_count = nchar(V5, type = "chars")) %>%
  mutate(variant_type = ifelse(reference_count > alternative_count, 
                               "deletion", 
                               "insertion")) %>%
  mutate(SVid = paste(SVid = paste(V1, V2, sep = "-"))) %>%
  select(-V3:-V21, -alternative_count) 

bed_tissue <- function(tissue, result) {
  unified_peaks_tissue <- unified_peaks %>%
    filter(str_detect(X4, tissue))
  
  tissue_result <- result[[tissue]]
  
  eQTLs_bed_tissue <- tissue_result %>%
    mutate(chrom = str_extract(snps, "^[^-]+"),
           position = as.numeric(str_extract(snps, "[^-]+$"))) %>%
    mutate(start = position) %>%
    select(-statistic, -FDR, -beta, -position) %>%
    select(chrom, start, everything())
  
  print(head(eQTLs_bed_tissue))
  
  eQTLs_tissue <- eQTLs_bed_tissue %>%
    inner_join(vfc_tissues, by = join_by(snps == SVid) ) %>%
    mutate(end = ifelse(variant_type == "insertion", start, start+reference_count )) %>%
    select(-reference_count) %>%
    mutate(start = start-1) %>%
    select(-V1,-V2, -tissue) %>%
    select(chrom, start, end, everything())
  
  
  unified_peaks_file <- paste0("regions_", tolower(tissue), ".bed")
  eQTLs_file <- paste0("eQTLs_", tolower(tissue), ".bed")
  
  write_delim(unified_peaks_tissue, unified_peaks_file, delim = "\t", col_names = FALSE)
  write_delim(eQTLs_tissue, eQTLs_file, delim = "\t", col_names = FALSE)
  
  message(paste("Saved:", unified_peaks_file, "and", eQTLs_file))
}

result <- list(
  Brain = brain_results, DistalIntestine = distalintestine_results, Gill = gill_results,
  Gonad = gonad_results, HeadKidney = headkidney_results, Liver = liver_results,
  Muscle = muscle_results
)


tissues <- c("Brain", "DistalIntestine", "Gill", 
             "Gonad", "HeadKidney", "Liver", 
             "Muscle")

for (tissue_name in tissues) {
  bed_tissue(tissue_name, result)
}

print

## do the same for eQTL results 0bp ##
bed_tissue <- function(tissue, result) {
  tissue_result <- result[[tissue]]
  
  eQTLs_bed_tissue <- tissue_result %>%
    mutate(chrom = str_extract(snps, "^[^-]+"),
           position = as.numeric(str_extract(snps, "[^-]+$"))) %>%
    mutate(start = position) %>%
    select(-statistic, -FDR, -beta, -position) %>%
    select(chrom, start, everything())
  
  print(head(eQTLs_bed_tissue))
  
  eQTLs_tissue <- eQTLs_bed_tissue %>%
    inner_join(vfc_tissues, by = join_by(snps == SVid) ) %>%
    mutate(end = ifelse(variant_type == "insertion", start, start+reference_count )) %>%
    select(-reference_count) %>%
    mutate(start = start-1) %>%
    select(-V1,-V2, -tissue) %>%
    select(chrom, start, end, everything())
  
  eQTLs_file <- paste0("eQTLs_0bp_", tolower(tissue), ".bed")
  write_delim(eQTLs_tissue, eQTLs_file, delim = "\t", col_names = FALSE)
  
  message(paste("Saved:", eQTLs_0bp_file))
}

result <- list(
  Brain = brain_results_0bp, DistalIntestine = distalintestine_results_0bp, Gill = gill_results_0bp,
  Gonad = gonad_results_0bp, HeadKidney = headkidney_results_0bp, Liver = liver_results_0bp,
  Muscle = muscle_results_0bp
)


tissues <- c("Brain", "DistalIntestine", "Gill", 
             "Gonad", "HeadKidney", "Liver", 
             "Muscle")

for (tissue_name in tissues) {
  bed_tissue(tissue_name, result)
}

print

