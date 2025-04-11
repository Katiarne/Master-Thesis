library(tidyverse)
library("MatrixEQTL")

# format expression file
tissues <- list(
  brain = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/RNA/Brain/results/star_salmon/salmon.merged.gene_tpm.tsv",
  distalintestine = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/RNA/DistalIntestine/results/star_salmon/salmon.merged.gene_tpm.tsv",
  gill = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/RNA/Gill/results/star_salmon/salmon.merged.gene_tpm.tsv",
  gonad = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/RNA/Gonad/results/star_salmon/salmon.merged.gene_tpm.tsv",
  headkidney = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/RNA/HeadKidney/results/star_salmon/salmon.merged.gene_tpm.tsv",
  liver = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/RNA/Liver/results/star_salmon/salmon.merged.gene_tpm.tsv",
  muscle = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/RNA/Muscle/results/star_salmon/salmon.merged.gene_tpm.tsv"
)

read_and_process <- function(file) {
  data <- read.delim(file, header = TRUE, sep = "\t") %>%
    select(-2) %>%
    mutate(across(-gene_id, log1p))
  
  colnames(data)[-1] <- gsub("_(Brain|DistalIntestine|Gill|Gonad|HeadKidney|Liver|Muscle)_", "_", colnames(data)[-1])
  data
}

dfs <- lapply(tissues, read_and_process)

combined_GE <- dfs[[1]]
for (df in dfs[-1]) {
  combined_GE[-1] <- combined_GE[-1] + df[-1]
}

GE <- combined_GE %>%
  mutate(across(-gene_id, log1p)) %>%
  rowwise() %>%
  ungroup() %>%
  # filter(rowSums(select(., -gene_id) == 0) < 11) %>%
  mutate(total_expression = rowSums(across(-gene_id)))  %>%
  filter(total_expression >= quantile(total_expression, 0.3)) %>%
  select(-total_expression)


# format SV dataset
aqua_faang_table <- tibble(
  sample_name = c(
    "AtlanticSalmon_RNA_Brain_Immature_Female_R1",
    "AtlanticSalmon_RNA_Brain_Immature_Female_R2",
    "AtlanticSalmon_RNA_Brain_Immature_Female_R3",
    "AtlanticSalmon_RNA_Brain_Immature_Male_R1",
    "AtlanticSalmon_RNA_Brain_Immature_Male_R2",
    "AtlanticSalmon_RNA_Brain_Immature_Male_R3",
    "AtlanticSalmon_RNA_Brain_Mature_Female_R1",
    "AtlanticSalmon_RNA_Brain_Mature_Female_R2",
    "AtlanticSalmon_RNA_Brain_Mature_Female_R3",
    "AtlanticSalmon_RNA_Brain_Mature_Male_R1",
    "AtlanticSalmon_RNA_Brain_Mature_Male_R2",
    "AtlanticSalmon_RNA_Brain_Mature_Male_R3"
  ),
  sample_id = c(
    "C1310101", "C1310201", "C1310301",
    "C1320101", "C1320201", "C1320301",
    "C1330401", "C1330501", "C1330601",
    "C1340401", "C1340501", "C1340601"
  ),
  fish_id = c("J8", "J9", "J10", "J15", "J16", "J18", "A14", "A15", "A24", "A19", "A20", "A22")
)


vcf_headers <- read_lines("/mnt/SCRATCH/kristenl/katinka_master/Ssalv3.1_pangenome_pangenie_all_svs_aquafaang_samples.chr.vcf")

column_headers <- vcf_headers %>%
  str_subset("^#CHROM") %>%           # Find the line starting with '#CHROM'
  str_split_fixed("\t", n = Inf) %>%  # Split into column names by tab
  as.character()

vcf_file <- read.delim("/mnt/SCRATCH/kristenl/katinka_master/Ssalv3.1_pangenome_pangenie_all_svs_aquafaang_samples.chr.vcf", 
                       header = FALSE, comment.char = "#", sep = "\t") %>%
  set_names(column_headers)

fish_ids <- column_headers[10:length(column_headers)] 

new_headers <- aqua_faang_table %>%
  filter(fish_id %in% fish_ids) %>%          
  arrange(match(fish_id, fish_ids)) %>%      
  pull(sample_name) 

final_headers <- c(column_headers[1:9], new_headers)

vcf_file <- vcf_file %>%
  set_names(final_headers)

SV_data <- vcf_file %>%
  mutate(SVid = paste(SVid = paste(`#CHROM`, POS, sep = "-"))) %>%
  select(SVid, final_headers[10:21]) %>%
  setNames(c("SVid", final_headers[10:21])) %>%
  mutate(across(-SVid, as.character)) %>%
  mutate(across(final_headers[10:21], ~ case_when(
    . == "0/0" ~ "0",      # Homozygous reference
    . %in% c("0/1", "1/0") ~ "1",   # Heterozygous
    . == "1/1" ~ "2",      # Homozygous alternate
    TRUE ~ NA_character_   # If non-standard, replace with NA
  )))

colnames(SV_data) <- gsub("_Brain", "", colnames(SV_data))

sv_samples <- colnames(SV_data)[-1]
ge_samples <- colnames(GE)[-1]

common_samples <- intersect(sv_samples, ge_samples)

SV_data <- SV_data %>%
  select(SVid, all_of(common_samples))

GE <- GE %>%
  select(gene_id, all_of(common_samples))  

# format covarities

headers_GE <- colnames(GE)[-1]
covariates_tibble <- tibble(
  sample_id = headers_GE,
  gender = if_else(str_detect(headers_GE, "Male"), 1, 0),          # 1 for Male, 0 for Female
  maturity = if_else(str_detect(headers_GE, "Mature"), 1, 0)       # 1 for Mature, 0 for Immature
)

covariates.txt <- covariates_tibble %>%
  pivot_longer(cols = c(gender, maturity), names_to = "id", values_to = "value") %>%  # Rename covariate to id
  pivot_wider(names_from = sample_id, values_from = value)

#SV location file

svpos <- SV_data %>%
  mutate(chromosome = str_extract(SVid, "^[^-]+"),
         position = str_extract(SVid, "[^-]+$")) %>%
  select(SVid, chromosome, position)

#Gene location file
gff_data <- read.delim("/mnt/SCRATCH/kristenl/katinka_master/Salmo_salar.Ssal_v3.1.107.chr.gff3", 
                       header = FALSE, comment.char = "#", sep = "\t") %>% 
  mutate(gene_id = str_extract(V9, "ENSSSAG[0-9]+")) %>%  
  select(gene_id, chr = V1, start = V4, end = V5)             

genpos <- GE %>%
  select(gene_id) %>%
  left_join(gff_data, by = "gene_id") %>%
  filter(!is.na(start) & !is.na(end)) 

# Save data
write.table(SV_data, file = "SV_data.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(GE, file = "GE_combined.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(covariates.txt, file = "covariates.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(svpos, file = "svpos.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(genpos, file = "genpos.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# eQTL Analysis
useModel = modelLINEAR 
SNP_file_name = "SV_data.txt"
expression_file_name = "GE_combined.txt"
covariates_file_name = "covariates.txt"
snps_location_file_name = "svpos.txt"
gene_location_file_name = "genpos.txt"

output_file_name_cis = tempfile()

pvOutputThreshold_cis = 0.01

errorCovariance = numeric()

cisDist = 0 #1e4 #10kb

snps = SlicedData$new()
snps$fileDelimiter = "\t"      # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 2000      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name )

gene = SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
gene$LoadFile( expression_file_name )

cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1       # one column of row labels
cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
cvrt$LoadFile( covariates_file_name )

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  output_file_name = NULL,
  pvOutputThreshol = 1e4, #10kb
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

## Make the histogram of local and distant p-values
plot(me)

results_tpm_sum_0bp <- me$cis$eqtls %>%
  filter(!is.infinite(statistic) & !is.nan(statistic))

unique_genes <- results_tpm_sum %>%
  distinct(gene, .keep_all = TRUE) 

write.csv(results_tpm_sum, "eQTL_tpm_sum_0bp.csv", row.names = FALSE)
