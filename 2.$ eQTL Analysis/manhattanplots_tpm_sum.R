library(readr)
library(ggplot2)
library(ggtext)

resutls_tpm_sum <- read_csv("results_tpm_sum")
cis_results_all <- results_tpm_sum

cis_results <- cis_results_all %>%
  mutate(chromosome = as.numeric(sub("-.*", "", snps)),
         position = as.numeric(sub(".*-", "", snps))) %>%
  rename(SVid = snps)

cis_results_cum <- cis_results %>%
  group_by(chromosome) %>%
  summarise(max_bp = max(position)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(chromosome, bp_add)

cis_results <- cis_results %>%
  inner_join(cis_results_cum, by = "chromosome") %>%
  mutate(bp_cum = position + bp_add)

axis_set <- cis_results %>%
  group_by(chromosome) %>%
  summarise(center = mean(bp_cum))

ylim <- cis_results %>%
  filter(pvalue == min(pvalue)) %>%
  mutate(ylim = abs(floor(log10(pvalue))) + 2) %>%
  pull(ylim)

sig <- 0.01 / nrow(cis_results)


manhplot <- ggplot(cis_results, aes(
  x = bp_cum, y = -log10(pvalue),
  color = as_factor(chromosome), size = -log10(pvalue)
)) +
  # geom_hline(
  #   yintercept = -log10(sig), color = "grey40",
  #   linetype = "dashed"
  # ) +
  geom_point(alpha = 0.75, size = 1.5) +
  scale_x_continuous(
    label = axis_set$chromosome,
    breaks = axis_set$center
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(2, ylim)) +
  scale_color_manual(values = rep(
    c("#276FBF", "#183059"),
    unique(length(axis_set$chromosome))
  )) +
  scale_size_continuous(range = c(0.5, 3)) +
  labs(
    x = NULL,
    y = "-log<sub>10</sub>(pvalue)"
  ) +
  ggtitle("Genes Associated with eQTLs in All Tissues") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) 

print(manhplot)
