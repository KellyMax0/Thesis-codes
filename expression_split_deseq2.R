suppressPackageStartupMessages({
  library(DESeq2); library(data.table); library(dplyr); library(readr)
})

# Load data
cnt <- as.matrix(data.frame(fread("expr_tumour_dedup.csv"), row.names = 1, check.names = FALSE))
colnames(cnt) <- substr(colnames(cnt), 1, 12)
gm <- fread("gene_mapping.tsv")
pg <- read_csv("patient_groups.csv", show_col_types = FALSE)
mt <- read_csv("waterfall_mutation_types.csv", show_col_types = FALSE) %>%
  select(gene, patient_id, mutation_type)

# DESeq2 normalisation
dds <- DESeqDataSetFromMatrix(countData = round(cnt),
                              colData = data.frame(row.names = colnames(cnt), x = rep(1, ncol(cnt))),
                              design = ~1)
dds <- estimateSizeFactors(dds)
norm <- counts(dds, normalized = TRUE)
ens <- setNames(gm$gene_id, gm$gene_name)[c("SMARCA4", "ARID1A")]

# Split mutants at the SWI/SNF WT median
wt_rna <- intersect(pg$patient_id[pg$group == "SWI/SNF WT"], colnames(cnt))
expr <- bind_rows(lapply(names(ens), function(g) {
  mut <- intersect(pg$patient_id[pg$group == paste0(g, "-mutant")], colnames(cnt))
  med <- median(norm[ens[[g]], wt_rna])
  data.frame(patient_id = mut, gene = g, deseq2_count = norm[ens[[g]], mut],
             wt_median_deseq2 = med, diff = norm[ens[[g]], mut] - med)
})) %>%
  mutate(expression_group = ifelse(diff > 0, "Higher", "Lower")) %>%
  left_join(mt, by = c("gene", "patient_id"))
print(expr %>% count(gene, expression_group))

# Save
dir.create("smarca4_split_normalised", showWarnings = FALSE)
write_csv(expr, "smarca4_split_normalised/expression_subgroups_deseq2.csv")
