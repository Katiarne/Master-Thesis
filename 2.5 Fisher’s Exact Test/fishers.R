library(tidyverse)

# make a datafram with length of SV 
vcf <- vcf_file <- read.delim("/mnt/SCRATCH/kristenl/katinka_master/Ssalv3.1_pangenome_pangenie_all_svs_aquafaang_samples.chr.vcf", 
                              header = FALSE, comment.char = "#", sep = "\t") %>%
  mutate(reference_count = nchar(V4, type = "chars"),
         alternative_count = nchar(V5, type = "chars")) %>%
  mutate(variant_type = ifelse(reference_count > alternative_count, 
                               "deletion", 
                               "insertion")) %>%
  mutate(SVid = paste(SVid = paste(V1, V2, sep = "-"))) %>%
  mutate(SV_length = ifelse(variant_type == "deletion", reference_count, alternative_count)) %>%
  select(-V1:-V21, -alternative_count, -reference_count) 


eQTLresult <- results_tpm_sum %>%
  full_join(vcf, by = join_by(snps == SVid)) %>%
  mutate(eQTL_association = ifelse(is.na(gene), "NO", "YES")) %>%
  mutate(SV_category = ifelse(SV_length < 100, "short", "long"))%>%
  select(eQTL_association, SV_category)

contingency_table <- table(eQTLresult$SV_category, eQTLresult$eQTL_association)
print(contingency_table)

fisher_result <- fisher.test(contingency_table)
print(fisher_result)

## Deletion or insertion
type_test <- bedtools %>%
  mutate(eQTL = ifelse(!is.na(V6.y), "YES", "NO")) %>%
  # mutate(exon = ifelse(!str_detect(V4.x, "insertion") 
  #                      , FALSE, TRUE)) %>%
  select(eQTL, V4)

tabel_type <- table(type_test$eQTL, type_test$V4)
print(tabel_type)

fisher_result_type <- fisher.test(tabel_type)
print(fisher_result_type)

### fisher's exact test using ATAC-file information ###
bedtools_intersect <- read.delim("intersect_output.bed", header = FALSE) %>%
  distinct(V4,V5, .keep_all = TRUE)
bedtools_intersect_all <- read.delim("intersect_output_all.bed", header = FALSE) %>%
  distinct(V5, .keep_all = TRUE)

bedtools <- bedtools_intersect_all %>%
  left_join(bedtools_intersect, by = c("V5" = "V4"))

betools_TSS <- bedtools %>%
  filter(str_detect(V12.x, "Active TSS"))
  
bedtools_fishers <- bedtools %>%
  mutate(eQTL = ifelse(!is.na(V6.y), "YES", "NO")) %>%
  mutate(biological_region = ifelse(!str_detect(V12.x, "\\.") , TRUE, FALSE)) %>%
  select(eQTL, biological_region)

contingency_table_bio <- table(bedtools_fishers$eQTL, bedtools_fishers$biological_region)
print(contingency_table_bio)

fisher_result_bio <- fisher.test(contingency_table_bio)
print(fisher_result_bio)

# FISHERS on Active TSS
TSS_test <- bedtools %>%
  mutate(eQTL = ifelse(!is.na(V6.y), "YES", "NO")) %>%
  mutate(TSS = ifelse(!str_detect(V12.x, "Active TSS") , FALSE, TRUE)) %>%
  select(eQTL, TSS)

tabel_TSS <- table(TSS_test$eQTL, TSS_test$TSS)
print(tabel_TSS)

fisher_result_TSS <- fisher.test(tabel_TSS)
print(fisher_result_TSS)

# FISHER on ENHANCER
enhancer_test <- bedtools %>%
  mutate(eQTL = ifelse(!is.na(V6.y), "YES", "NO")) %>%
  mutate(enhancer = ifelse(!str_detect(V12.x, "Active enhancer|Primed enhancer|Poised enhancer 1|Poised enhancer 2") 
                           , FALSE, TRUE)) %>%
  select(eQTL, enhancer)

tabel_enhancer <- table(enhancer_test$eQTL, enhancer_test$enhancer)
print(tabel_enhancer)

fisher_result_enhancer <- fisher.test(tabel_enhancer)
print(fisher_result_enhancer)

# FISHER on accesible chromatin/ATAC
AC_test <- bedtools %>%
  mutate(eQTL = ifelse(!is.na(V6.y), "YES", "NO")) %>%
  mutate(AC = ifelse(!str_detect(V12.x, "Accessible chromatin"), FALSE, TRUE)) %>%
  select(eQTL, AC)

tabel_AC <- table(AC_test$eQTL, AC_test$AC)
print(tabel_AC)

fisher_result_AC <- fisher.test(tabel_AC)
print(fisher_result_AC)

# Fisher on CDS, Exon and UTR
bedtools_intersect_eQTL_gff <- read.delim("intersect_gff_complete.bed", header = FALSE) %>%
  distinct(V4,V5, V11, .keep_all = TRUE)
bedtools_intersect_gff <- read.delim("intersect_output_gff.bed", header = FALSE) %>%
  distinct(V5, V9, .keep_all = TRUE)

bedtools <- bedtools_intersect_gff %>%
  left_join(bedtools_intersect_eQTL_gff, by = c("V5" = "V4", "V9" = "V11")) %>%
  group_by(V4, V5.y) %>%
  mutate(region_priority = case_when(
    V9 == "CDS" ~ 1,
    grepl("UTR", V9) ~ 2, 
    V9 == "exon" ~ 3  
  )) %>%
  ungroup() %>% 
  arrange(region_priority) %>%
  distinct(V5, V5.y, .keep_all = TRUE)

CDS_all <- bedtools %>%
  filter(!is.na(V6.y) & V6.y != ".") %>%
  filter(str_detect(V9, "CDS"))

# Fishers on CDS
CDS_test <- bedtools %>%
  mutate(eQTL = ifelse(!is.na(V6.y), "YES", "NO")) %>%
  mutate(CDS = ifelse(str_detect(V9, "CDS"), TRUE, FALSE)) %>%
  select(eQTL, CDS)


tabel_CDS <- table(CDS_test$eQTL, CDS_test$CDS)
print(tabel_CDS)

fisher_result_CDS <- fisher.test(tabel_CDS)
print(fisher_result_CDS)


# Fishers on exon
exon_test <- bedtools %>%
  mutate(eQTL = ifelse(!is.na(V6.y), "YES", "NO")) %>%
  mutate(exon = ifelse(!str_detect(V9, "exon") 
                      , FALSE, TRUE)) %>%
  select(eQTL, exon)

tabel_exon <- table(exon_test$eQTL, exon_test$exon)
print(tabel_exon)

fisher_result_exon <- fisher.test(tabel_exon)
print(fisher_result_exon)

# Fisher on UTR
UTR_test <- bedtools %>%
  mutate(eQTL = ifelse(!is.na(V6.y), "YES", "NO")) %>%
  mutate(UTR = ifelse(!str_detect(V9, "five_prime_UTR|three_prime_UTR") 
                       , FALSE, TRUE)) %>%
  select(eQTL, UTR)

tabel_UTR <- table(UTR_test$eQTL, UTR_test$UTR)
print(tabel_UTR)

fisher_result_UTR <- fisher.test(tabel_UTR)
print(fisher_result_UTR)
