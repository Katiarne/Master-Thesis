# bed file med gene navn
library(tidyverse)

gff_data <- read.delim("/mnt/SCRATCH/kristenl/katinka_master/Salmo_salar.Ssal_v3.1.107.chr.gff3", 
                       header = FALSE, comment.char = "#", sep = "\t") %>%
  select(-V2, -V3, -V6:-V8) %>%
  mutate(gene_name = str_replace(str_extract(V9, "description=([^;\\[]+)"), "description=", "")) %>%
  filter(!is.na(gene_name)) %>%
  mutate(gene_id = str_extract(V9, "ENSSSAG[0-9]+")) %>%
  select(V1, V4, V5, gene_id, gene_name) %>%
  mutate(V4 = V4-1)

eQTL_gene <- results_tpm_sum %>%
  mutate(chrom = str_extract(snps, "^[^-]+"),
         position = as.numeric(str_extract(snps, "[^-]+$"))) %>%
  mutate(start = position) %>%
  select(-statistic, -FDR, -beta, -position) %>%
  select(chrom, start, everything())

gene_name_eQTL <- eQTL_gene %>%
  left_join(gff_data, by = c("gene" = "gene_id"))
