library(tidyverse)
library(MatrixEQTL)

# List of tissue URLs
tissues <- list(
  Brain = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/RNA/Brain/results/star_salmon/salmon.merged.gene_tpm.tsv",
  DistalIntestine = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/RNA/DistalIntestine/results/star_salmon/salmon.merged.gene_tpm.tsv",
  Gill = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/RNA/Gill/results/star_salmon/salmon.merged.gene_tpm.tsv",
  Gonad = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/RNA/Gonad/results/star_salmon/salmon.merged.gene_tpm.tsv",
  HeadKidney = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/RNA/HeadKidney/results/star_salmon/salmon.merged.gene_tpm.tsv",
  Liver = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/RNA/Liver/results/star_salmon/salmon.merged.gene_tpm.tsv",
  Muscle = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/RNA/Muscle/results/star_salmon/salmon.merged.gene_tpm.tsv"
)

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

thresholds <- list(
    Brain = 0.45,
    DistalIntestine = 0.45,
    Gill = 0.45,
    Gonad = 0.45,
    HeadKidney = 0.50,
    Liver = 0.6,
    Muscle = 0.55
    )


# Function to process each tissue
process_tissue <- function(tissue_name, tissue_url) {
  threshold <- thresholds[[tissue_name]]
  GE.txt <- read.delim(tissue_url, header = TRUE, sep = "\t") %>%
    select(-gene_name) %>%
    mutate(across(-gene_id, log1p))  %>%
    rowwise() %>%
    ungroup() %>%
    mutate(total_expression = rowSums(select(., -gene_id), na.rm = TRUE)) %>%
    filter(total_expression >= quantile(total_expression, threshold))
    
  print(hist(GE.txt$total_expression)) # to decide threshold

  GE.txt <- GE.txt %>%
    select(-total_expression)

  vcf_headers <- read_lines("/mnt/SCRATCH/kristenl/katinka_master/Ssalv3.1_pangenome_pangenie_all_svs_aquafaang_samples.chr.vcf")
  
  column_headers <- vcf_headers %>%
    str_subset("^#CHROM") %>%           
    str_split_fixed("\t", n = Inf) %>%  
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
    select(SVid, (final_headers[10:21])) %>%
    setNames(c("SVid", final_headers[10:21])) %>%
    mutate(across(-SVid, as.character)) %>%
    mutate(across(final_headers[10:21], ~ case_when(
      . == "0/0" ~ "0",      # Homozygous reference
      . %in% c("0/1", "1/0") ~ "1",   # Heterozygous
      . == "1/1" ~ "2",      # Homozygous alternate
      TRUE ~ NA_character_   # If non-standard, replace with NA
    ))) %>%
    mutate(across(final_headers[10:21], as.numeric)) 
  
  colnames(SV_data) <- gsub("_Brain", "", colnames(SV_data))
  tissue_pattern <- paste0("_", tissue_name)
  colnames(GE.txt) <- gsub(tissue_pattern, "", colnames(GE.txt))
  
  sv_samples <- colnames(SV_data)[-1]
  ge_samples <- colnames(GE.txt)[-1]      
  
  common_samples <- intersect(sv_samples, ge_samples)
  
  SV_data <- SV_data %>%
    select(SVid, all_of(common_samples))
  
  GE.txt <- GE.txt %>%
    select(gene_id, all_of(common_samples))
  
  headers_GE <- colnames(GE.txt)[-1]
  covariates_tibble <- tibble(
    sample_id = headers_GE,
    gender = if_else(str_detect(headers_GE, "Male"), 1, 0),
    maturity = if_else(str_detect(headers_GE, "Mature"), 1, 0)
  )
  covariates.txt <- covariates_tibble %>%
    pivot_longer(cols = c(gender, maturity), names_to = "id", values_to = "value") %>%
    pivot_wider(names_from = sample_id, values_from = value)
  
  svpos <- SV_data %>%
    mutate(chromosome = str_extract(SVid, "^[^-]+"),
           position = str_extract(SVid, "[^-]+$")) %>%
    select(SVid, chromosome, position)
  
  gff_data <- read.delim("/mnt/SCRATCH/kristenl/katinka_master/Salmo_salar.Ssal_v3.1.107.chr.gff3", 
                         header = FALSE, comment.char = "#", sep = "\t") %>%
    mutate(gene_id = str_extract(V9, "ENSSSAG[0-9]+")) %>%
    select(gene_id, chr = V1, start = V4, end = V5)
  
  genpos <- GE.txt %>%
    select(gene_id) %>%
    left_join(gff_data, by = "gene_id") %>%
    filter(!is.na(start) & !is.na(end))

  # eQTL analysis start
  write.table(SV_data, file = paste0("SV_data_", tissue_name, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(GE.txt, file = paste0("GE_", tissue_name, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(covariates.txt, file = paste0("covariates_", tissue_name, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(svpos, file = paste0("svpos_", tissue_name, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(genpos, file = paste0("genpos_", tissue_name, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  SNP_file_name = paste0("SV_data_", tissue_name, ".txt")
  expression_file_name = paste0("GE_", tissue_name, ".txt")
  covariates_file_name = paste0("covariates_", tissue_name, ".txt")
  snps_location_file_name = paste0("svpos_", tissue_name, ".txt")
  gene_location_file_name = paste0("genpos_", tissue_name, ".txt")
  
  output_file_name_cis = tempfile()
  
  snps = SlicedData$new()
  snps$fileDelimiter = "\t"  # Tab-separated file
  snps$fileOmitCharacters = "NA"  # Treat "NA" as missing
  snps$fileSkipRows = 1  # Skip the header row
  snps$fileSkipColumns = 1  # Skip the first column (SVid)
  snps$fileSliceSize = 2000  # Number of rows per slice
  
  snps$LoadFile(SNP_file_name)
  
  gene = SlicedData$new()
  gene$fileDelimiter = "\t"
  gene$fileOmitCharacters = "NA"
  gene$fileSkipRows = 1
  gene$fileSkipColumns = 1
  gene$fileSliceSize = 2000
  gene$LoadFile(expression_file_name)
  
  cvrt = SlicedData$new()
  cvrt$fileDelimiter = "\t"
  cvrt$fileOmitCharacters = "NA"
  cvrt$fileSkipRows = 1
  cvrt$fileSkipColumns = 1
  cvrt$fileSliceSize = 2000
  cvrt$LoadFile(covariates_file_name)
  
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)
  
  me <- Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    useModel = modelLINEAR,
    errorCovariance = numeric(),
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = 0.01,
    output_file_name = NULL,
    pvOutputThreshol = 0,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = 1e4, #10kb
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
  )
  
  # Extract significant results
  cis_results <- me$cis$eqtls
  
  # Save results to a file
  write.table(cis_results, file = paste0("cis_results_", tissue_name, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(cis_results)
}

# Loop through each tissue
results <- lapply(names(tissues), function(tissue_name) {
  process_tissue(tissue_name, tissues[[tissue_name]])
})

# `results` now contains the eQTL results for each tissue

brain_results <- results[[1]]  %>%
  filter(!is.infinite(statistic) & !is.nan(statistic)) %>%
  mutate(tissue = "brain")

distalintestine_results <- results[[2]] %>%
  filter(!is.infinite(statistic) & !is.nan(statistic)) %>%
  mutate(tissue = "distalintestine")

gill_results <- results[[3]] %>%
  filter(!is.infinite(statistic) & !is.nan(statistic)) %>%
  mutate(tissue = "gill")

gonad_results <- results[[4]] %>%
  filter(!is.infinite(statistic) & !is.nan(statistic)) %>%
  mutate(tissue = "gonad")

headkidney_results <- results[[5]] %>%
  filter(!is.infinite(statistic) & !is.nan(statistic)) %>%
  mutate(tissue = "headkidney")

liver_results <- results[[6]] %>%
  filter(!is.infinite(statistic) & !is.nan(statistic)) %>%
  mutate(tissue = "liver")

muscle_results <- results[[7]] %>%
  filter(!is.infinite(statistic) & !is.nan(statistic)) %>%
  mutate(tissue = "muscle")

results_tissues <- bind_rows(
  brain_results,
  distalintestine_results,
  gill_results,
  gonad_results,
  headkidney_results,
  liver_results,
  muscle_results
)

# Save file
write.csv(results_tissues, "eQTL_tissues.csv", row.names = FALSE)


