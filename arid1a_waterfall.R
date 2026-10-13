suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(readr) })

# Load data
expr <- read_csv("smarca4_split_normalised/expression_subgroups_deseq2.csv", show_col_types = FALSE)

# Plot settings
BASE_FONT   <- if (Sys.info()[["sysname"]] == "Darwin") "Helvetica" else "Liberation Sans"
STATUS_COLS <- c("Lower than WT median" = "#6CB4E4", "Higher than WT median" = "#C8642A")
TYPE_COLS   <- c("Frame_Shift_Del" = "#3B75AF", "Frame_Shift_Ins" = "#6A3D9A",
                 "Missense_Mutation" = "#56A03C", "Nonsense_Mutation" = "#D33A2F",
                 "Splice_Site" = "#F08C35")
X_LAB <- "Difference from WT median (DESeq2-normalised counts)"

# Waterfall plot function
waterfall <- function(d, fill_var, cols, legend_title, gene, n_wt = 310) {
  d <- d %>% arrange(diff) %>% mutate(patient_id = factor(patient_id, levels = patient_id))
  ggplot(d, aes(x = diff, y = patient_id, fill = .data[[fill_var]])) +
    geom_col(width = 0.85, colour = "grey30", linewidth = 0.15) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    scale_fill_manual(values = cols, name = legend_title, drop = TRUE) +
    labs(title = paste0(gene, "-mutant tumours: expression relative to the SWI/SNF WT median"),
         subtitle = paste0("Baseline (0) = median of ", n_wt, " SWI/SNF WT primary tumours"),
         x = X_LAB, y = paste0(gene, "-mutant tumours (n = ", nrow(d), ")")) +
    theme_minimal(base_size = 14, base_family = BASE_FONT) +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(), legend.position = "bottom",
          plot.title = element_text(size = 15), legend.text = element_text(size = 12))
}

# ARID1A data
d <- expr %>% filter(gene == "ARID1A") %>%
  mutate(status = ifelse(diff > 0, "Higher than WT median", "Lower than WT median"))

# Fig S6A: ARID1A waterfall by expression group
ggsave("figS6A_waterfall_arid1a_status.png",
       waterfall(d, "status", STATUS_COLS, "Expression", "ARID1A"),
       width = 9, height = 7.5, dpi = 300, bg = "white")

# Fig S6B: ARID1A waterfall by mutation type
ggsave("figS6B_waterfall_arid1a_mutation_type.png",
       waterfall(d, "mutation_type", TYPE_COLS, "Mutation type", "ARID1A"),
       width = 9, height = 7.5, dpi = 300, bg = "white")
