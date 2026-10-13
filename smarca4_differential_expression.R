suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(dplyr); library(readr)
})
BASE_FONT  <- if (Sys.info()[["sysname"]] == "Darwin") "Helvetica" else "Liberation Sans"
GROUP_COLS <- c("SMARCA4-mutant" = "#E63946", "SWI/SNF WT" = "#457B9D")

# Samples
groups <- read_csv("patient_groups.csv", show_col_types = FALSE) %>%
  filter(in_rna_cohort_505 == "yes", group %in% c("SMARCA4-mutant", "SWI/SNF WT")) %>%
  transmute(sample    = tumour_rnaseq_sample,
            condition = factor(ifelse(group == "SWI/SNF WT", "WT", "SMARCA4"),
                               levels = c("WT", "SMARCA4")))
print(table(groups$condition))

# Counts and DESeq2
cts <- read.csv("expr_tumour_dedup.csv", row.names = 1, check.names = FALSE)
cts <- as.matrix(cts[, groups$sample])
coldata <- data.frame(condition = groups$condition, row.names = groups$sample)
keep <- rowSums(cts >= 10) >= min(table(coldata$condition))
dds <- DESeqDataSetFromMatrix(countData = cts[keep, ], colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "SMARCA4", "WT"), alpha = 0.05)

# Results table
gene_map <- read.delim("gene_mapping.tsv") %>% distinct(gene_id, .keep_all = TRUE)
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  left_join(gene_map, by = "gene_id") %>%
  relocate(gene_name, .after = gene_id) %>%
  arrange(padj)
write_csv(res_df, "mutation_pathway.csv")

# Ranked list for GSEA
rnk <- res_df %>%
  filter(!is.na(stat), !is.na(gene_name), gene_name != "") %>%
  group_by(gene_name) %>% slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup() %>%
  arrange(desc(stat)) %>% select(gene_name, stat)
write_tsv(rnk, "smarca4_vs_swi_wt.rnk", col_names = FALSE)

# Fig 14A: PCA
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"), 1)
pca <- pca %>%
  mutate(group = ifelse(condition == "WT", "SWI/SNF WT", "SMARCA4-mutant")) %>%
  arrange(desc(group == "SWI/SNF WT"))
n_lab <- c("SMARCA4-mutant" = paste0("SMARCA4-mutant (n = ", sum(pca$group == "SMARCA4-mutant"), ")"),
           "SWI/SNF WT"     = paste0("SWI/SNF WT (n = ",     sum(pca$group == "SWI/SNF WT"), ")"))
p_pca <- ggplot(pca, aes(PC1, PC2, colour = group)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_colour_manual(values = GROUP_COLS, labels = n_lab, name = NULL) +
  labs(x = paste0("PC1 (", pv[1], "% variance)"), y = paste0("PC2 (", pv[2], "% variance)")) +
  theme_bw(base_size = 16, base_family = BASE_FONT) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = "black", size = 14),
        legend.position = "top", legend.text = element_text(size = 14))
ggsave("fig14A_pca_smarca4_vs_swi_wt.png", p_pca, width = 7.5, height = 6.5, dpi = 300, bg = "white")
ggsave("fig14A_pca_smarca4_vs_swi_wt.pdf", p_pca, width = 7.5, height = 6.5, bg = "white", device = cairo_pdf)

# Fig 14B: volcano plot
vol <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(status = case_when(padj < 0.05 & log2FoldChange >  1 ~ "Enriched",
                            padj < 0.05 & log2FoldChange < -1 ~ "Depleted",
                            TRUE ~ "Not significant"))
n_up <- sum(vol$status == "Enriched"); n_dn <- sum(vol$status == "Depleted")
p_vol <- ggplot(vol, aes(log2FoldChange, -log10(padj), colour = status)) +
  geom_point(size = 1.1, alpha = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5) +
  scale_colour_manual(values = c(Enriched = "#E41A1C", Depleted = "#1F3DD6", `Not significant` = "grey75"),
                      guide = "none") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1.6, size = 6.5, family = BASE_FONT,
           label = paste0("Depleted: ", format(n_dn, big.mark = ","))) +
  annotate("text", x = Inf, y = Inf, hjust = 1.15, vjust = 1.6, size = 6.5, family = BASE_FONT,
           label = paste0("Enriched: ", format(n_up, big.mark = ","))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = "log\u2082 fold change (SMARCA4-mutant vs SWI/SNF WT)",
       y = "\u2212log\u2081\u2080 adjusted p-value") +
  theme_bw(base_size = 18, base_family = BASE_FONT) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = "black", size = 15))
ggsave("fig14B_volcano_smarca4_vs_swi_wt.png", p_vol, width = 8, height = 6.5, dpi = 300, bg = "white")
