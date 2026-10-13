suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(ggrepel)
})

# Groups, colours and cell compartments
S4 <- "SMARCA4-mutant"; WT <- "SWI/SNF WT"
COL <- c("SMARCA4-mutant" = "#E63946", "SWI/SNF WT" = "#457B9D")
IMMUNE <- c("NK", "CD4 Th", "Naive CD4 T", "Treg", "CD8 T", "B cells", "Plasma cells",
            "Monocytes", "mo-Mac", "Alveolar Mac", "Dendritic cells", "Mast cells")
COMPARTMENT <- c(setNames(rep("Lymphoid", 7), c("Naive CD4 T", "CD4 Th", "Treg", "CD8 T", "NK", "B cells", "Plasma cells")),
                 setNames(rep("Myeloid", 5), c("Monocytes", "mo-Mac", "Alveolar Mac", "Dendritic cells", "Mast cells")),
                 setNames(rep("Epithelial", 6), c("Malignant", "Transitional epi", "AT1", "AT2", "Ciliated", "Club")),
                 setNames(rep("Stromal", 3), c("Fibroblasts", "Endothelial", "Mural")))
CD8_STATES <- c("CD8_Naive", "CD8_Cytotoxic", "CD8_Transitional", "CD8_Exhausted", "CD8_ProgenitorExhausted")
CD8_LABEL  <- c("Naive", "Cytotoxic", "Transitional", "Exhausted", "Progenitor\nexhausted")
META <- c("Mixture", "P-value", "Correlation", "RMSE", "Absolute score (sig.score)", "patient_id", "group")

# Load data
key <- read_csv("CIBERSORTx_rerun/CIBERSORTx_mixture_sample_key.csv", show_col_types = FALSE)
N   <- table(key$group)
LAB <- c("SMARCA4-mutant" = paste0("SMARCA4-mutant (n = ", N[[S4]], ")"),
         "SWI/SNF WT"     = paste0("SWI/SNF WT (n = ", N[[WT]], ")"))
load_job <- function(f) read_tsv(f, show_col_types = FALSE) %>% inner_join(key, by = "Mixture") %>%
  filter(group %in% c(S4, WT)) %>% mutate(group = factor(group, levels = c(S4, WT)))
cd8 <- load_job("CIBERSORTx_rerun/CIBERSORTx_Job51_Adjusted.txt")

# Statistics
fmt_p <- function(p) ifelse(p < 0.001, sprintf("%.1e", p), ifelse(p < 0.1, sprintf("%.3f", p), sprintf("%.2f", p)))
wilcox_p <- function(a, b) {
  if (length(unique(c(a, b))) == 1) return(NA_real_)
  wilcox.test(a, b, exact = FALSE, correct = TRUE)$p.value
}
compare <- function(d, cols) {
  map_dfr(cols, function(c) {
    a <- d[[c]][d$group == S4]; b <- d[[c]][d$group == WT]
    tibble(cell_type = c, median_smarca4 = median(a), median_wt = median(b),
           mean_smarca4 = mean(a), mean_wt = mean(b), p_value = wilcox_p(a, b))
  }) %>% mutate(q_value = p.adjust(p_value, method = "BH"), significant = q_value < 0.05)
}

# Shared plot elements
theme_panel <- function(base = 10) theme_classic(base_size = base) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid.major.y = element_line(colour = "#e5e5e5", linewidth = 0.4),
        legend.position = "top", legend.justification = "right")
box_layers <- function(d, x) {
  med <- d %>% group_by(.data[[x]], group) %>% summarise(m = median(value), .groups = "drop")
  list(geom_boxplot(position = position_dodge(0.8), width = 0.72, alpha = 0.35, outlier.shape = NA,
                    linewidth = 0.4),
       geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8, seed = 1),
                  size = 0.6, alpha = 0.45, stroke = 0),
       geom_errorbar(data = med, aes(x = .data[[x]], ymin = m, ymax = m, group = group),
                     inherit.aes = FALSE, position = position_dodge(0.8), width = 0.72,
                     colour = "black", linewidth = 0.6),
       scale_fill_manual(values = COL, labels = LAB, name = NULL),
       scale_colour_manual(values = COL, labels = LAB, name = NULL))
}

# Pan-cell figure function
pancell_figure <- function(job_file, with_cd8) {
  pan   <- load_job(job_file)
  cells <- setdiff(names(pan), META)
  pan <- pan %>% mutate(immune_total = rowSums(across(all_of(IMMUNE))),
                        immune_share = immune_total / `Absolute score (sig.score)`)

  # A: log2 ratio vs adjusted p
  st <- compare(pan, cells) %>%
    mutate(log2FC_medians = log2((median_smarca4 + 2e-4) / (median_wt + 2e-4)),
           log2FC_means = log2(mean_smarca4 / mean_wt),
           compartment = COMPARTMENT[cell_type])
  use_means <- any(st$median_smarca4 == 0 | st$median_wt == 0)
  st <- st %>% mutate(fc = if (use_means) log2FC_means else log2FC_medians,
                      colour = ifelse(!significant, "ns",
                                      ifelse(fc > 0, "Higher in SMARCA4-mutant", "Higher in SWI/SNF WT")))
  p_a <- ggplot(st, aes(fc, -log10(q_value), colour = colour)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#888888", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "#888888", linewidth = 0.4) +
    geom_point(size = 3) +
    geom_text_repel(data = filter(st, significant), aes(label = cell_type), colour = "black",
                    size = 3.2, min.segment.length = 0, segment.colour = "#999999", seed = 1) +
    scale_colour_manual(values = c("Higher in SMARCA4-mutant" = COL[[S4]],
                                   "Higher in SWI/SNF WT" = COL[[WT]], "ns" = "#b8b8b8"),
                        breaks = c("Higher in SMARCA4-mutant", "Higher in SWI/SNF WT", "ns"),
                        labels = c("Higher in SMARCA4-mutant", "Higher in SWI/SNF WT", expression(q >= 0.05)),
                        name = NULL) +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    labs(x = paste0(ifelse(use_means, "log2 ratio of means", "log2 ratio of medians"),
                    "\n(SMARCA4-mutant / SWI/SNF WT)"),
         y = "-log10 adjusted p-value") +
    theme_panel() + theme(panel.grid.major.x = element_line(colour = "#eeeeee", linewidth = 0.4))

  # B: mean relative fraction
  rel <- pan %>% mutate(across(all_of(cells), ~ .x / `Absolute score (sig.score)`)) %>%
    group_by(group) %>% summarise(across(all_of(cells), mean), .groups = "drop") %>%
    pivot_longer(-group, names_to = "cell_type", values_to = "mean_rel") %>%
    mutate(compartment = factor(COMPARTMENT[cell_type], levels = c("Lymphoid", "Myeloid", "Epithelial", "Stromal")))
  order_b <- rel %>% filter(group == WT) %>% arrange(compartment, desc(mean_rel)) %>% pull(cell_type)
  rel <- rel %>% mutate(cell_type = factor(cell_type, levels = rev(order_b)),
                        group = factor(group, levels = c(WT, S4)))
  p_b <- ggplot(rel, aes(mean_rel, cell_type, fill = group)) +
    geom_col(position = position_dodge(0.8), width = 0.76) +
    facet_grid(compartment ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = COL, labels = LAB, name = NULL, breaks = c(S4, WT)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Mean relative fraction (score / total score)", y = NULL) +
    theme_classic(base_size = 10) +
    theme(axis.text = element_text(colour = "black"), legend.position = "top",
          legend.justification = "right", strip.background = element_blank(),
          strip.text.y = element_text(angle = 0, face = "italic", colour = "#555555", hjust = 0),
          panel.grid.major.x = element_line(colour = "#eeeeee", linewidth = 0.4),
          panel.spacing = unit(2, "pt"))

  # C: total and immune scores
  lab_total  <- paste0("Total score\n(all ", length(cells), " cell types)")
  lab_immune <- "Immune score\n(12 immune cell types)"
  st_c <- compare(pan, c("Absolute score (sig.score)", "immune_total", "immune_share"))
  d_c <- pan %>% select(group, `Absolute score (sig.score)`, immune_total) %>%
    pivot_longer(-group, names_to = "score", values_to = "value") %>%
    mutate(score = factor(ifelse(score == "immune_total", lab_immune, lab_total), levels = c(lab_total, lab_immune)))
  top_c <- max(d_c$value)
  txt_c <- tibble(score = factor(c(lab_total, lab_immune), levels = c(lab_total, lab_immune)),
                  lab = paste0("Wilcoxon p = ", fmt_p(st_c$p_value[1:2])))
  p_c <- ggplot(d_c, aes(score, value, fill = group, colour = group)) +
    box_layers(d_c, "score") +
    geom_text(data = txt_c, aes(score, top_c * 1.03, label = lab), inherit.aes = FALSE, size = 3.3, vjust = 0) +
    scale_y_continuous(limits = c(NA, top_c * 1.12)) +
    labs(x = NULL, y = "Absolute score") + theme_panel()

  out <- list(stats_a = st, stats_c = st_c, p_a = p_a, p_b = p_b, p_c = p_c)

  # D: CD8 T-cell states
  if (with_cd8) {
    st_d <- compare(cd8, CD8_STATES)
    d_d <- cd8 %>% select(group, all_of(CD8_STATES)) %>%
      pivot_longer(-group, names_to = "state", values_to = "value") %>%
      mutate(state = factor(state, levels = CD8_STATES, labels = CD8_LABEL))
    top_d <- max(d_d$value)
    txt_d <- st_d %>% mutate(state = factor(cell_type, levels = CD8_STATES, labels = CD8_LABEL),
                             lab = ifelse(significant, ifelse(q_value < 0.001, "***", ifelse(q_value < 0.01, "**", "*")),
                                          paste0("q = ", fmt_p(q_value))))
    out$p_d <- ggplot(d_d, aes(state, value, fill = group, colour = group)) +
      box_layers(d_d, "state") +
      geom_text(data = txt_d, aes(state, top_d * 1.03, label = lab), inherit.aes = FALSE,
                size = 3, colour = "#555555", vjust = 0) +
      scale_y_continuous(limits = c(NA, top_d * 1.12)) +
      labs(x = NULL, y = expression("Fraction of CD8"^"+"~"T cells")) + theme_panel()
    out$stats_d <- st_d
  }
  out
}

# Fig 16A-D: all-sample pan-cell matrix (Job 50) and CD8 matrix (Job 51)
f16 <- pancell_figure("CIBERSORTx_rerun/CIBERSORTx_Job50_Adjusted.txt", with_cd8 = TRUE)
p16 <- ((f16$p_a | f16$p_b) / f16$p_c / f16$p_d) +
  plot_layout(heights = c(1.15, 0.85, 0.9)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 18, face = "bold"))
ggsave("fig16_pancell_cd8.png", p16, width = 13, height = 15, dpi = 300, bg = "white")
write_csv(f16$stats_a, "fig16A_pancell_stats.csv")
write_csv(f16$stats_c, "fig16C_total_scores_stats.csv")
write_csv(f16$stats_d, "fig16D_cd8_states_stats.csv")

# Fig S15A-C: primary-tumour pan-cell matrix (Job 48)
fS15 <- pancell_figure("CIBERSORTx_rerun/CIBERSORTx_Job48_Adjusted.txt", with_cd8 = FALSE)
pS15 <- ((fS15$p_a | fS15$p_b) / fS15$p_c) +
  plot_layout(heights = c(1.15, 0.85)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 18, face = "bold"))
ggsave("figS15_pancell_primary_matrix.png", pS15, width = 13, height = 10.5, dpi = 300, bg = "white")
write_csv(fS15$stats_a, "figS15A_pancell_primary_stats.csv")
write_csv(fS15$stats_c, "figS15C_total_scores_primary_stats.csv")
