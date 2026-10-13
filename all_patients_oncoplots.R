suppressPackageStartupMessages({
  library(tidyverse); library(maftools)
})

# Patients (n = 505)
cohort_ids <- read.csv("patient_groups.csv") %>%
  filter(in_rna_cohort_505 == "yes") %>%
  pull(patient_id)

# MAF: one primary tumour per patient
maf_df <- data.table::fread("tcga-luad_maf.tsv", skip = "Hugo_Symbol") %>%
  mutate(Patient_ID  = substr(Tumor_Sample_Barcode, 1, 12),
         sample_type = substr(Tumor_Sample_Barcode, 14, 16)) %>%
  filter(Patient_ID %in% cohort_ids, sample_type %in% c("01A", "01B")) %>%
  group_by(Patient_ID) %>%
  filter(sample_type == "01A" | !any(sample_type == "01A")) %>%
  ungroup() %>%
  mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 16)) %>%
  distinct(Tumor_Sample_Barcode, Chromosome, Start_Position,
           Reference_Allele, Tumor_Seq_Allele2, .keep_all = TRUE) %>%
  select(-Patient_ID, -sample_type)
cohort_maf <- read.maf(maf = maf_df)

# Gene lists
common_oncogenes <- c(
  "ARID1A", "TP53", "KRAS", "EGFR", "STK11", "KEAP1",
  "RB1", "CDKN2A", "NF1", "MET", "BRAF", "PIK3CA", "ATM", "SMARCA4",
  "RBM10", "U2AF1", "SETD2", "SMAD4", "PTEN"
)
swi_snf_genes <- c(
  "SMARCA4", "ARID1A", "ARID2", "BCL11B", "ARID1B",
  "BCL11A", "SMARCA2", "PBRM1", "DPF3", "SMARCC2",
  "BRD9", "ACTL6B", "DPF2", "SMARCC1", "ACTL6A",
  "SMARCD2", "SS18", "SMARCD1", "DPF1", "SMARCD3",
  "BRD7", "SMARCE1", "SMARCB1", "PHF10", "BCL7C",
  "BCL7A", "BCL7B", "ACTB", "SS18L1", "BICRA",
  "BICRAL"
)

# Plot functions
draw_oncoplot <- function(genes, title) {
  oncoplot(
    maf                = cohort_maf,
    genes              = genes,
    sortByMutation     = TRUE,
    gene_mar           = 8,
    titleText          = title,
    fontSize           = 1.0,
    titleFontSize      = 1.4,
    legendFontSize     = 1.4,
    annotationFontSize = 1.3,
    legend_height      = 5
  )
}
save_oncoplot <- function(genes, title, file_stub, height = 8) {
  pdf(paste0(file_stub, ".pdf"), width = 12, height = height)
  draw_oncoplot(genes, title); dev.off()
  png(paste0(file_stub, ".png"), width = 12, height = height, units = "in", res = 300)
  draw_oncoplot(genes, title); dev.off()
}

# Fig 10A: common oncogenes
save_oncoplot(common_oncogenes, "All Patients - Common Oncogenes",
              "fig10A_oncoplot_all_patients_common_oncogenes")

# Fig 10B: SWI/SNF genes
save_oncoplot(swi_snf_genes, "All Patients - SWI/SNF Genes",
              "fig10B_oncoplot_all_patients_swisnf", height = 10)
