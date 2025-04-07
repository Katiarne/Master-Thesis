# List of tissues and their corresponding number of eQTLs
eQTL_numb <- c("brain" = 8855, "liver" = 5930, "gonad" = 7996, "muscle" = 6716)
eQTL_numb_CDS_UTR <- c("brain" = 6457, "liver" = 4060, "gonad" = 5747, "muscle" = 4659)

combined_data <- data.frame()

# Loop through each tissue and perform the analysis
for (tissue in names(eQTL_numb)) {
  
  # Read the data for each tissue
  bedtools_intersect <- read.delim(paste0("bedtools/intersect_output_", tissue, ".bed"), header = FALSE) %>%
    distinct(V5, .keep_all = TRUE)
  
  bedtools_intersect_eQTL <- read.delim(paste0("bedtools/intersect_output_eQTL_", tissue, ".bed"), header = FALSE) %>%
    distinct(V4, V5, .keep_all = TRUE)
  
  # Combine the data
  bedtools <- bedtools_intersect %>%
    left_join(bedtools_intersect_eQTL, by = c("V5" = "V4"))
  
  # Fisher's test on Biological Region
  bedtools_fishers <- bedtools %>%
    mutate(eQTL = ifelse(!is.na(V6.y), "YES", "NO")) %>%
    mutate(biological_region = ifelse(!str_detect(V12.x, "\\.") , TRUE, FALSE)) %>%
    select(eQTL, biological_region)
  
  contingency_table_bio <- table(bedtools_fishers$eQTL, bedtools_fishers$biological_region)
  print(contingency_table_bio)
  fisher_result_bio_brain <- fisher.test(contingency_table_bio)
  
  # Fisher's test on TSS
  TSS_test <- bedtools %>%
    mutate(eQTL = ifelse(!is.na(V6.y), "YES", "NO")) %>%
    mutate(TSS = ifelse(!str_detect(V12.x, "Active TSS") , FALSE, TRUE)) %>%
    select(eQTL, TSS)
  
  tabel_TSS <- table(TSS_test$eQTL, TSS_test$TSS)
  fisher_result_TSS_brain <- fisher.test(tabel_TSS_brain)
  
  # Fisher's test on Enhancer
  bedtools_enhancer <- bedtools %>%
    mutate(eQTL = ifelse(!is.na(V6.y), "YES", "NO")) %>%
    mutate(CDS = ifelse(!str_detect(V12.x, "Active enhancer|Primed enhancer|Poised enhancer 1|Poised enhancer 2") , FALSE, TRUE)) %>%
    select(eQTL, CDS)
  
  contingency_table_enhancer <- table(bedtools_enhancer $eQTL, bedtools_enhancer$CDS)
  fisher_result_enhancer_B <- fisher.test(contingency_table_enhancer_B)
  
  # Fisher's test on CDS
  intersect_eQTL <- read.delim(paste0("bedtools/intersect_gff_", tissue, ".bed"), header = FALSE) %>%
    distinct(V4, V5, V11, .keep_all = TRUE)
  
  bedtools_intersect_gff <- read.delim(paste0("intersect_output_gff.bed"), header = FALSE) %>%
    distinct(V5, V9, .keep_all = TRUE)
  
  bedtools <- bedtools_intersect_gff %>%
    left_join(intersect_eQTL, by = c("V5" = "V4", "V9" = "V11")) %>%
    group_by(V4, V5.y) %>%
    mutate(region_priority = case_when(
      V9 == "CDS" ~ 1,
      grepl("UTR", V9) ~ 2, 
      V9 == "exon" ~ 3  
    )) %>%
    ungroup() %>% 
    arrange(region_priority) %>%
    distinct(V5, V5.y, .keep_all = TRUE)
  
  bedtools_CDS <- bedtools %>%
    mutate(eQTL = ifelse(!is.na(V6.y), "YES", "NO")) %>%
    mutate(CDS = ifelse(str_detect(V9, "CDS"), TRUE, FALSE)) %>%
    select(eQTL, CDS)
  
  contingency_table_CDS <- table(bedtools_CDS$eQTL, bedtools_CDS$CDS)
  print(contingency_table_CDS_B)
  fisher_result_CDS_B <- fisher.test(contingency_table_CDS_B)

  # Fisher's test on UTR
  bedtools_UTR <- bedtools %>%
    mutate(eQTL = ifelse(!is.na(V6.y), "YES", "NO")) %>%
    mutate(UTR = ifelse(!str_detect(V9, "five_prme_UTR|three_prime_UTR") , FALSE, TRUE)) %>%
    select(eQTL, UTR)
  
  contingency_table_UTR <- table(bedtools_UTR$eQTL, bedtools_UTR$UTR)
  fisher_result_UTR_B <- fisher.test(contingency_table_UTR_B)
  
  # Figure #
  label <- c("ATAC", "TSS", "Enhancer", "CDS", "UTR")
  number_eQTL <- c(contingency_table_bio["YES", "TRUE"], tabel_TSS["YES", "TRUE"],
                   contingency_table_enhancer["YES", "TRUE"], contingency_table_CDS["YES", "TRUE"],
                   contingency_table_UTR["YES", "TRUE"])
  conf_int <- c(fisher_result_bio$conf.int,
                fisher_result_TSS$conf.int, fisher_result_enhancer$conf.int, fisher_result_CDS$conf.int, 
                fisher_result_UTR$conf.int)
  odds_ratio <- c(fisher_result_bio$estimate,
                  fisher_result_TSS$estimate, fisher_result_enhancer$estimate,
                  fisher_result_CDS$estimate, 
                  fisher_result_UTR$estimate)
  percent_eQTL <- round(
    c(number_eQTL[1:3] / eQTL_numb[tissue],  
      number_eQTL[4:5] / eQTL_numb_CDS_UTR[tissue]) * 100, 2)
  fig_lab <- c("ATAC", "TSS", "Enhancers", "CDS", "UTR")
  
  conf_matrix <- matrix(conf_int, ncol = 2, byrow = TRUE)
  lower <- conf_matrix[, 1]
  upper <- conf_matrix[, 2]
  
  # Combine all the data into a single dataframe for all tissues
  df <- data.frame(label, odds_ratio, lower, upper, percent_eQTL, fig_lab)
  df$eQTL_count <- eQTL_numb[tissue]
  df$tissue <- tissue 
  
  combined_data <- rbind(combined_data, df)  # Add the data for this tissue to the combined dataframe
}

## Figure ##

library(ggplot2)

# Convert 'label' to factor
combined_data$label <- factor(combined_data$label, levels = c("UTR", "CDS", "Enhancer", "TSS", "ATAC"))

# Create the plot using facet_wrap
fisher_tissues <- ggplot(combined_data, aes(x = label, y = odds_ratio, ymin = lower, ymax = upper)) +
  geom_pointrange() + 
  geom_hline(yintercept = 1, lty = 2) + 
  coord_flip() +  
  xlab("") + ylab("Odds Ratio (95% Confidence Interval)") + 
  facet_wrap(~ tissue, scales = "free_y") + 
  ggtitle("Enrichment of eQTL-Associated SVs Across Genomic and Regulatory Features - Tissues") +  
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 8),
    axis.title.y = element_text(size = 12), 
    strip.text = element_text(size = 10)  # Control the size of facet labels
  ) +
  geom_text(aes(label = paste0(percent_eQTL, "% ", fig_lab), 
                y = upper + 0.2), hjust = 0, size = 3) +
  ylim(min(combined_data$lower) * 0.8, 5.5)

print(fisher_tissues)

