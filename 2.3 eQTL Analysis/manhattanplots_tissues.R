library(readr)
library(ggplot2)
library(ggtext)

# Retain the `tissue` column in cis_results
cis_results <- results_tissues %>%
  mutate(
    chromosome = as.numeric(sub("-.*", "", snps)),
    position = as.numeric(sub(".*-", "", snps))
  ) %>%
  rename(SVid = snps)

# Create cumulative positions
cis_results_cum <- cis_results %>%
  group_by(chromosome) %>%
  summarise(max_bp = max(position)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(chromosome, bp_add)

# Merge cumulative positions and compute `bp_cum`
cis_results <- cis_results %>%
  inner_join(cis_results_cum, by = "chromosome") %>%
  mutate(bp_cum = position + bp_add)

# Set up axis labels
axis_set <- cis_results %>%
  group_by(chromosome) %>%
  summarise(center = mean(bp_cum))

# Determine ylim
ylim <- cis_results %>%
  filter(pvalue == min(pvalue)) %>%
  mutate(ylim = abs(floor(log10(pvalue))) + 2) %>%
  pull(ylim)

# Set significance threshold
sig <- 0.05 / nrow(cis_results)

tissue_labels <- c(
  "brain" = "Brain",
  "distalintestine" = "Distal Intestine",
  "gill" = "Gill",
  "gonad" = "Gonad",
  "headkidney" = "Head Kidney",
  "liver" = "Liver",
  "muscle" = "Muscle"
)

manhplot <- ggplot(cis_results, aes(
  x = bp_cum,
  y = -log10(pvalue),
  color = as_factor(chromosome), 
)) +
  geom_point(alpha = 0.75) +
  scale_x_continuous(
    labels = axis_set$chromosome[seq(1, length(axis_set$chromosome), by = 2)], 
    breaks = axis_set$center[seq(1, length(axis_set$center), by = 2)]
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(2, ylim), breaks = pretty(c(2, ylim), n = 5)) +
  scale_color_manual(values = rep(
    c("#276FBF", "#183059"), length.out = length(unique(axis_set$chromosome))
  )) +
  scale_size_continuous(range = c(0.5, 3)) +
  facet_wrap(~tissue, ncol = 4, scales = "free_x", labeller = labeller(tissue = tissue_labels)) + 
  labs(
    x = NULL,
    y = "-log<sub>10</sub>(pvalue)"
  ) +
  ggtitle("Genes Associated with eQTLs Across Tissues") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 4, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    strip.placement = "outside",
    strip.background = element_rect(fill = "lightgrey", color = "black"),  
    strip.text = element_text(size = 10)  
  )

print(manhplot)
