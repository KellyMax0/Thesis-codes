suppressPackageStartupMessages({
  library(tidyverse); library(patchwork)
})

# Load data
key  <- read_csv("CIBERSORTx_rerun/CIBERSORTx_mixture_sample_key.csv", show_col_types = FALSE)
lm22 <- read_tsv("CIBERSORTx_rerun/CIBERSORTx_Job49_Adjusted.txt", show_col_types = FALSE) %>%
  inner_join(key, by = "Mixture")
cells <- setdiff(names(lm22), c("Mixture", "P-value", "Correlation", "RMSE", "patient_id", "group"))

# Shared colour scale
vmax <- quantile(as.matrix(lm22[cells]), 0.995, type = 7)
GROUPS <- c("SMARCA4-mutant" = "SMARCA4-mutant", "Other SWI/SNF-mutant" = "SWI MUT",
            "SWI/SNF WT" = "SWI/SNF WT")

# Heatmap function
heat_panel <- function(g) {
  d <- lm22 %>% filter(group == g) %>% mutate(row = row_number()) %>%
    select(row, all_of(cells)) %>%
    pivot_longer(-row, names_to = "cell_type", values_to = "fraction") %>%
    mutate(cell_type = factor(cell_type, levels = cells))
  ggplot(d, aes(cell_type, row, fill = fraction)) +
    geom_raster() +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "Reds"), limits = c(0, vmax),
                         oob = scales::squish, name = "Fraction of\nimmune cells") +
    scale_y_reverse(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    labs(title = paste0(GROUPS[[g]], " (n = ", max(d$row), ")"), x = NULL, y = "Tumours") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5, colour = "black"),
          axis.text.y = element_blank(), panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
}

# Fig S13: per-tumour LM22 heatmaps
p <- wrap_plots(lapply(names(GROUPS), heat_panel), nrow = 1, guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 18, face = "bold"))
ggsave("figS13_lm22_heatmaps.png", p, width = 16, height = 9, dpi = 300, bg = "white")
