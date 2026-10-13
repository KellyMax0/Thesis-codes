suppressPackageStartupMessages({
  library(tidyverse)
})

# Load data
key  <- read_csv("CIBERSORTx_rerun/CIBERSORTx_mixture_sample_key.csv", show_col_types = FALSE)
lm22 <- read_tsv("CIBERSORTx_rerun/CIBERSORTx_Job49_Adjusted.txt", show_col_types = FALSE) %>%
  inner_join(key, by = "Mixture")
cells <- setdiff(names(lm22), c("Mixture", "P-value", "Correlation", "RMSE", "patient_id", "group"))

# Groups and colours
SM <- "Other SWI/SNF-mutant"; WT <- "SWI/SNF WT"
COL <- c("Other SWI/SNF-mutant" = "#E63946", "SWI/SNF WT" = "#457B9D")
two <- lm22 %>% filter(group %in% c(SM, WT))

# Statistics
stars <- function(q) ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", "*"))
wilcox_p <- function(a, b) {
  if (length(unique(c(a, b))) == 1) return(NA_real_)
  wilcox.test(a, b, exact = FALSE, correct = TRUE)$p.value
}
stats <- map_dfr(cells, function(c) {
  a <- two[[c]][two$group == SM]; b <- two[[c]][two$group == WT]
  tibble(cell_type = c, median_swi_mut = median(a), median_wt = median(b),
         mean_swi_mut = mean(a), mean_wt = mean(b), p_value = wilcox_p(a, b))
}) %>%
  mutate(q_value = p.adjust(p_value, method = "BH"),
         significant = q_value < 0.05,
         panel = paste0(cell_type, "\n",
                        ifelse(is.na(q_value), "all zero",
                               ifelse(q_value < 0.05, stars(q_value), sprintf("q = %.2f", q_value)))))
write_csv(stats %>% select(-panel), "figS14_lm22_swi_mut_vs_swi_wt_stats.csv")

# Fig S14: SWI MUT vs SWI/SNF WT, 22 LM22 cell types
d <- two %>% select(group, all_of(cells)) %>%
  pivot_longer(-group, names_to = "cell_type", values_to = "fraction") %>%
  left_join(stats %>% select(cell_type, panel), by = "cell_type") %>%
  mutate(panel = factor(panel, levels = stats$panel),
         group = factor(ifelse(group == SM, "SWI MUT", "WT"), levels = c("SWI MUT", "WT")))
med <- d %>% group_by(panel, group) %>% summarise(m = median(fraction), .groups = "drop")
p <- ggplot(d, aes(group, fraction, fill = group, colour = group)) +
  geom_boxplot(width = 0.55, alpha = 0.35, outlier.shape = NA, linewidth = 0.4) +
  geom_errorbar(data = med, aes(x = group, ymin = m, ymax = m), inherit.aes = FALSE,
                width = 0.55, colour = "black", linewidth = 0.5) +
  facet_wrap(~ panel, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = c("SWI MUT" = COL[[SM]], "WT" = COL[[WT]]), guide = "none") +
  scale_colour_manual(values = c("SWI MUT" = COL[[SM]], "WT" = COL[[WT]]), guide = "none") +
  labs(x = NULL, y = "Fraction of immune cells") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 8.5),
        axis.text = element_text(colour = "black", size = 7.5))
ggsave("figS14_lm22_swi_mut_vs_swi_wt.png", p, width = 12, height = 16, dpi = 300, bg = "white")
