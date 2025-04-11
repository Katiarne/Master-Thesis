## Regression ###
sv_id <- "17-12377744" # enter SV-id you want to look at
gene_id <- "ENSSSAG00000102446" # enter gene-id you want to look at

gene_expression <- read.delim("GE_Liver.txt") %>% # Gene expression file from eQTL analysis, change depending on the tissue you want to look at
  as.data.frame()

genotypes <- read.delim("SV_data_Liver.txt") %>% # SV file from eQTL analysis, change depending on the tissue you want to look at
  as.data.frame()

sv_row_index <- which(genotypes[, 1] == sv_id)
gene_row_index <- which(gene_expression[, 1] == gene_id)


e1 = as.numeric(gene_expression[gene_row_index, -1]) 
s1 = as.numeric(genotypes[sv_row_index, -1])
lm1 = lm(e1 ~ s1)

gene_name <- as.character(gene_expression[gene_row_index, 1])
sv_name <- as.character(genotypes[sv_row_index, 1])
colors <- c("lightblue", "lightgreen", "lightpink")

plot(e1 ~ jitter(s1), col = colors[s1+1], pch = 16, xaxt="n",
     xlab="Genotype", ylab="Expression", 
     ylim = range(c(e1, lm1$fitted)))
axis(1,at=c(0:2),labels=c("0/0","0/1","1/1"))
lines(lm1$fitted ~ s1,type="b",pch=15,col="darkgrey")
title(main = paste(gene_name, "vs", sv_name), cex.main = 1.2)
