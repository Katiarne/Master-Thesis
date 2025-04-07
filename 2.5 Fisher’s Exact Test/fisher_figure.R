### Plot fishers exact test #####
library(ggplot2)

label <- c("Short vs Long SVs", "Deletion vs. Insertion", "ATAC", "TSS", "Enhancer",
           "CDS", "UTR")

conf_int <- c(fisher_result$conf.int,fisher_result_type$conf.int, fisher_result_bio$conf.int,
              fisher_result_TSS$conf.int, fisher_result_enhancer$conf.int, fisher_result_CDS$conf.int, 
              fisher_result_UTR$conf.int)
odds_ratio <- c(fisher_result$estimate, fisher_result_type$estimate, fisher_result_bio$estimate,
                fisher_result_TSS$estimate, fisher_result_enhancer$estimate,
                 fisher_result_CDS$estimate, fisher_result_UTR$estimate)

number_eQTL <- c(contingency_table["short", "YES"],tabel_type["YES", "deletion"], contingency_table_bio["YES", "TRUE"], 
                 tabel_TSS["YES", "TRUE"],tabel_enhancer["YES", "TRUE"], tabel_CDS["YES", "TRUE"],
                 tabel_UTR["YES", "TRUE"])

denominator <- c(rep(10434, 5), rep(6955, 2)) 
percent_eQTL <- round((number_eQTL / denominator) * 100, 2)

fig_lab <- c("Short SVs", "Deletions", "ATAC", "TSS", "Enhancers", "CDS", "UTR")

conf_matrix <- matrix(conf_int, ncol = 2, byrow = TRUE)
lower <- conf_matrix[, 1]
upper <- conf_matrix[, 2]


df <- data.frame(label, odds_ratio, lower, upper, percent_eQTL, fig_lab)
df$label <- factor(df$label, levels = rev(df$label))

fp <- ggplot(data=df, aes(x=label, y=odds_ratio, ymin=lower, ymax=upper)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) + 
  coord_flip() +  
  xlab("") + ylab("odds-ratio (95% confidence interval)") +
  ggtitle("Enrichment of eQTL-Associated SVs Across Genomic and Regulatory Features") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold")  # Make title bold and centered
  ) +
  geom_text(aes(label = paste0(percent_eQTL, "% ", fig_lab), 
                y = upper + 0.2), hjust = 0, size = 3) +  
  ylim(min(lower) * 0.8, max(upper) * 1.5)
print(fp)
