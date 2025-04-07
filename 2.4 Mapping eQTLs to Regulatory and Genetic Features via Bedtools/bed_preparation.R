library(tidyverse)

# prepare the eQTL data
eQTLs_bed <- results_tpm_sum %>%
  mutate(chrom = str_extract(snps, "^[^-]+"),
         position = as.numeric(str_extract(snps, "[^-]+$"))) %>%
  mutate(start = position) %>%
  select(-statistic, -FDR, -beta, -position) %>%
  select(chrom, start, everything())

# Load vcf file
vcf_file <- read.delim("/mnt/SCRATCH/kristenl/katinka_master/Ssalv3.1_pangenome_pangenie_all_svs_aquafaang_samples.chr.vcf", 
                       header = FALSE, comment.char = "#", sep = "\t") %>%
  mutate(reference_count = nchar(V4, type = "chars"),
         alternative_count = nchar(V5, type = "chars")) %>%
  mutate(variant_type = ifelse(reference_count > alternative_count, 
                               "deletion", 
                               "insertion")) %>%
  mutate(SVid = paste(SVid = paste(V1, V2, sep = "-"))) %>%
  select(-V3:-V21, -alternative_count) 

vcf_only <- read.delim("/mnt/SCRATCH/kristenl/katinka_master/Ssalv3.1_pangenome_pangenie_all_svs_aquafaang_samples.chr.vcf", 
                       header = FALSE, comment.char = "#", sep = "\t") 

vcf_bed <- vcf_file %>%
  rename(chrom = V1) %>%
  rename(start = V2) %>%
  mutate(end = ifelse(variant_type == "insertion", start, start+reference_count )) %>%
  select(-reference_count) %>%
  mutate(start = start-1) %>%
  select(chrom, start, end, everything())

eQTLs <- eQTLs_bed %>%
  inner_join(vcf_file, by = join_by(snps == SVid) ) %>%
  mutate(end = ifelse(variant_type == "insertion", start, start+reference_count )) %>%
  select(-reference_count) %>%
  mutate(start = start-1) %>%
  select(-V1, -V2) %>%
  select(chrom, start, end, everything())


#regulatory elements / unified peaks

unified_peaks <- read_tsv("https://salmobase.org/datafiles/datasets/Aqua-Faang/robust_ATAC_peaks/unified_annotated_peaks/AtlanticSalmon_unified_peaks.bed",
                          col_names = FALSE)


# enhancers and CDS 

gff <- read.delim("/mnt/SCRATCH/kristenl/katinka_master/Salmo_salar.Ssal_v3.1.107.chr.gff3", 
                  header = FALSE, comment.char = "#", sep = "\t") %>%
  select(-V2, -V6:-V9) %>%
  filter(V3 %in% c("CDS", "exon", "five_prime_UTR", "three_prime_UTR")) %>%
  select(V1, V4, V5, V3) %>%
  mutate(V4 = V4-1)

gff_gene_cordinates <- read.delim("/mnt/SCRATCH/kristenl/katinka_master/Salmo_salar.Ssal_v3.1.107.chr.gff3", 
                               header = FALSE, comment.char = "#", sep = "\t") %>%
  filter(V3 %in% c("gene")) %>%
  mutate(gene_id = str_extract(V9, "ENSSSAG[0-9]+")) %>%
  select(chr = V1, start = V4, end = V5, gene_id) %>%
  mutate(start = start-1)

# save bedfiles
write_delim(gff_gene_cordinates, "gff_gene_cordinates.bed", delim = "\t", col_names = FALSE)
write_delim(gff, "gff.bed", delim = "\t", col_names = FALSE)
write_delim(vcf_bed, "vcf.bed", delim = "\t", col_names = FALSE)
write_delim(unified_peaks, "regions.bed", delim = "\t", col_names = FALSE)
write_delim(eQTLs, "eQTLs.bed", delim = "\t", col_names = FALSE)

# Run bedtools
#
#
