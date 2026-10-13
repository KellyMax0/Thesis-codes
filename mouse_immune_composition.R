suppressPackageStartupMessages({
  library(tidyverse)
})

# Colours
CELL_COLS <- c("T cells" = "#2D153C", "CD4 T cells" = "#D0EE5B", "CD8 T cells" = "#4C61D0",
               "Tregs" = "#F0BD54", "NK cells" = "#59A7F2", "B cells" = "#E57335",
               "Macrophages" = "#6BE1B8", "Neutrophils" = "#BA3A1E", "Dendritic cells" = "#98FB75")
CELL_NAMES <- c("T Cell" = "T cells", "CD4 T Cell" = "CD4 T cells", "CD8 T Cell" = "CD8 T cells",
                "Treg" = "Tregs", "NK Cell" = "NK cells", "B Cell" = "B cells",
                "Macrophage" = "Macrophages", "Neutrophil" = "Neutrophils",
                "Dendritic cell" = "Dendritic cells")
BASE <- 14

# Fig 17B: immune composition of each mouse
counts <- read_csv("ScRNA/stacked_bar_data_scrna.csv", show_col_types = FALSE) %>%
  rename(cell_type = 1) %>%
  filter(cell_type %in% names(CELL_NAMES)) %>%
  mutate(cell_type = CELL_NAMES[cell_type])
long <- counts %>% pivot_longer(-cell_type, names_to = "mouse", values_to = "n") %>%
  group_by(mouse) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  mutate(genotype  = factor(ifelse(startsWith(mouse, "KPS"), "KPS", "KP"), levels = c("KP", "KPS")),
         cell_type = factor(cell_type, levels = names(CELL_COLS)))
p_b <- ggplot(long, aes(sub("^KPS?", "", mouse), pct, fill = fct_rev(cell_type))) +
  geom_col(width = 0.8, colour = "white", linewidth = 0.2) +
  facet_grid(~ genotype, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_manual(values = CELL_COLS, breaks = names(CELL_COLS), guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(c(0, 0.01))) +
  labs(x = "Mouse", y = "Proportion of immune cells") +
  theme_bw(base_size = BASE) +
  theme(strip.placement = "outside", strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = BASE + 1),
        axis.text.x = element_text(size = BASE - 2), panel.grid.major.x = element_blank(),
        panel.spacing = unit(0.8, "lines"))
ggsave("fig17B_mouse_immune_composition.png", p_b, width = 6, height = 5, dpi = 300, bg = "white")

# Fig 17C: mean proportion of total cells, KP vs KPS
raw <- read_csv("ScRNA/Cell_Types_ScRNA_2.csv", skip = 1, show_col_types = FALSE) %>% rename(cell_type = 1)
totals <- raw %>% filter(cell_type == "Total Cells") %>%
  pivot_longer(-cell_type, names_to = "sample", values_to = "total") %>%
  select(sample, total)
df <- raw %>%
  filter(cell_type != "Total Cells") %>%
  pivot_longer(-cell_type, names_to = "sample", values_to = "count") %>%
  left_join(totals, by = "sample") %>%
  mutate(prop  = count / total,
         group = if_else(str_starts(sample, "KPS"), "KPS", "KP"))
summ <- df %>%
  group_by(cell_type, group) %>%
  summarise(mean = mean(prop), sd = sd(prop), .groups = "drop")
lvl <- summ %>% group_by(cell_type) %>% summarise(t = sum(mean)) %>%
  arrange(desc(t)) %>% pull(cell_type)
summ$cell_type <- factor(summ$cell_type, lvl)
p_c <- ggplot(summ, aes(cell_type, mean, fill = group)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "white", linewidth = 0.25) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.8), width = 0.2, linewidth = 0.3) +
  scale_fill_manual(values = c(KP = "#4C72B0", KPS = "#DD8452"), name = "Group") +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, 0.05))) +
  labs(x = NULL, y = "Mean proportion of total cells", title = "Immune cells: KP vs KPS") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("fig17C_mouse_immune_kp_vs_kps.png", p_c, width = 14, height = 8, dpi = 300, bg = "white")
