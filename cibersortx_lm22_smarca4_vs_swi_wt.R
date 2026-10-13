suppressPackageStartupMessages({
  library(tidyverse); library(patchwork)
})

# Load data
key  <- read_csv("CIBERSORTx_rerun/CIBERSORTx_mixture_sample_key.csv", show_col_types = FALSE)
lm22 <- read_tsv("CIBERSORTx_rerun/CIBERSORTx_Job49_Adjusted.txt", show_col_types = FALSE) %>%
  inner_join(key, by = "Mixture")
cells <- setdiff(names(lm22), c("Mixture", "P-value", "Correlation", "RMSE", "patient_id", "group"))

# Groups and colours
S4 <- "SMARCA4-mutant"; WT <- "SWI/SNF WT"
N   <- table(key$group)
COL <- c("SMARCA4-mutant" = "#E63946", "SWI/SNF WT" = "#457B9D")
LAB <- c("SMARCA4-mutant" = paste0("SMARCA4-mutant (n = ", N[[S4]], ")"),
         "SWI/SNF WT"     = paste0("SWI/SNF WT (n = ", N[[WT]], ")"))
two <- lm22 %>% filter(group %in% c(S4, WT))

# Statistics
stars <- function(q) ifelse(is.na(q), "", ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", "ns"))))
wilcox_p <- function(a, b) {
  if (length(unique(c(a, b))) == 1) return(NA_real_)
  wilcox.test(a, b, exact = FALSE, correct = TRUE)$p.value
}
prev <- colMeans(two[cells] > 0)
fig_cells <- cells[prev >= 0.20]
stats <- map_dfr(fig_cells, function(c) {
  a <- two[[c]][two$group == S4]; b <- two[[c]][two$group == WT]
  tibble(cell_type = c, median_smarca4 = median(a), median_wt = median(b),
         mean_smarca4 = mean(a), mean_wt = mean(b), p_value = wilcox_p(a, b))
}) %>%
  mutate(q_value = p.adjust(p_value, method = "BH"),
         direction = case_when(median_smarca4 > median_wt ~ "higher in SMARCA4-mutant",
                               median_smarca4 < median_wt ~ "higher in SWI/SNF WT",
                               mean_smarca4 > mean_wt ~ "higher in SMARCA4-mutant (means)",
                               TRUE ~ "higher in SWI/SNF WT (means)"),
         significant = q_value < 0.05)
write_csv(stats, "fig15_lm22_smarca4_vs_swi_wt_stats.csv")

# Cell order
wt_med   <- sapply(fig_cells, function(c) median(two[[c]][two$group == WT]))
wt_order <- fig_cells[order(-wt_med)]
sig      <- wt_order[wt_order %in% stats$cell_type[stats$significant]]

# Box plot function
box_panel <- function(cell_set, base = 10) {
  d <- two %>% select(group, all_of(cell_set)) %>%
    pivot_longer(-group, names_to = "cell_type", values_to = "fraction") %>%
    mutate(cell_type = factor(cell_type, levels = cell_set),
           group = factor(group, levels = c(S4, WT)))
  med <- d %>% group_by(cell_type, group) %>% summarise(m = median(fraction), .groups = "drop")
  top <- max(d$fraction)
  st  <- stats %>% filter(cell_type %in% cell_set) %>%
    mutate(cell_type = factor(cell_type, levels = cell_set), lab = stars(q_value)) %>%
    filter(lab != "ns", lab != "")
  ggplot(d, aes(cell_type, fraction, fill = group, colour = group)) +
    geom_boxplot(position = position_dodge(0.8), width = 0.72, alpha = 0.35, outlier.shape = NA,
                 linewidth = 0.4) +
    geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8, seed = 1),
               size = 0.5, alpha = 0.45, stroke = 0) +
    geom_errorbar(data = med, aes(x = cell_type, ymin = m, ymax = m, group = group),
                  inherit.aes = FALSE, position = position_dodge(0.8), width = 0.72,
                  colour = "black", linewidth = 0.6) +
    geom_text(data = st, aes(x = cell_type, y = top * 1.03, label = lab), inherit.aes = FALSE,
              vjust = 0, fontface = "bold", size = base * 0.4) +
    scale_fill_manual(values = COL, labels = LAB, name = NULL) +
    scale_colour_manual(values = COL, labels = LAB, name = NULL) +
    scale_y_continuous(limits = c(0, top * 1.12), expand = expansion(mult = c(0.02, 0))) +
    labs(x = NULL, y = "Fraction of immune cells") +
    theme_classic(base_size = base) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
          panel.grid.major.y = element_line(colour = "#e5e5e5", linewidth = 0.4),
          legend.position = "top", legend.justification = "right",
          plot.margin = margin(5.5, 5.5, 5.5, 45))
}

# Fig 15A: all cell types present in >= 20% of tumours
p_a <- box_panel(wt_order)
ggsave("fig15A_lm22_smarca4_vs_swi_wt_all_cell_types.png", p_a, width = 13, height = 5.8, dpi = 300, bg = "white")

# Fig 15B: significant cell types
p_b <- box_panel(sig)
ggsave("fig15B_lm22_smarca4_vs_swi_wt_significant.png", p_b,
       width = max(4.5, 1.05 * length(sig) + 2), height = 5.2, dpi = 300, bg = "white")

# Fig 15: combined
p_ab <- (box_panel(wt_order, 12) / box_panel(sig, 12)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 20, face = "bold"))
ggsave("fig15_lm22_smarca4_vs_swi_wt_combined.png", p_ab, width = 9, height = 11.5, dpi = 300, bg = "white")
