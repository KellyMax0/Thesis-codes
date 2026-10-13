suppressPackageStartupMessages({ library(maftools); library(dplyr); library(readr) })

# Load data
maf  <- read.maf(maf = "tcga-luad_maf.tsv")
expr <- read_csv("smarca4_split_normalised/expression_subgroups_deseq2.csv", show_col_types = FALSE)
grp  <- function(g, lv = c("Higher", "Lower")) expr$patient_id[expr$gene == g & expr$expression_group %in% lv]

# One tumour barcode per patient
get_maf_barcodes <- function(patient_ids, maf_obj) {
  getSampleSummary(maf_obj) %>% pull(Tumor_Sample_Barcode) %>% as.character() %>%
    { data.frame(Tumor_Sample_Barcode = .) } %>%
    mutate(Patient_ID = substr(Tumor_Sample_Barcode, 1, 12)) %>%
    filter(Patient_ID %in% patient_ids) %>% group_by(Patient_ID) %>%
    filter(!(grepl("11A", Tumor_Sample_Barcode) & any(grepl("01A", Tumor_Sample_Barcode))),
           !(grepl("01B", Tumor_Sample_Barcode) & any(grepl("01A", Tumor_Sample_Barcode)))) %>%
    ungroup() %>% distinct(Patient_ID, .keep_all = TRUE) %>% pull(Tumor_Sample_Barcode)
}

# Colours and gene lists
VC_COLS <- c(Nonsense_Mutation = "#D1382C", Missense_Mutation = "#559E3E", Frame_Shift_Del = "#3C76B0",
             Frame_Shift_Ins = "#644095", Splice_Site = "#EF8632", In_Frame_Del = "#FFFDA6",
             In_Frame_Ins = "#A6CEE3", Translation_Start_Site = "#FB9A99", Nonstop_Mutation = "#CAB2D6",
             Multi_Hit = "#000000")
swi_snf_genes <- c("SMARCA4", "ARID1A", "ARID2", "BCL11B", "ARID1B", "BCL11A", "SMARCA2", "PBRM1",
                   "DPF3", "SMARCC2", "BRD9", "ACTL6B", "DPF2", "SMARCC1", "ACTL6A", "SMARCD2", "SS18",
                   "SMARCD1", "DPF1", "SMARCD3", "BRD7", "SMARCE1", "SMARCB1", "PHF10", "BCL7C",
                   "BCL7A", "BCL7B", "ACTB", "SS18L1", "BICRA", "BICRAL")

# Oncoplot function
onco_png <- function(patients, genes, title, file, w = 9, h = 7) {
  m <- subsetMaf(maf, tsb = get_maf_barcodes(patients, maf))
  png(file, width = w, height = h, units = "in", res = 300)
  oncoplot(maf = m, genes = genes, sortByMutation = TRUE, gene_mar = 8, titleText = title,
           colors = VC_COLS, bgCol = "#EDF0F1")
  dev.off()
}
common_oncogenes <- c("ARID1A", "TP53", "KRAS", "EGFR", "STK11", "KEAP1", "RB1", "CDKN2A", "NF1", "MET",
                      "BRAF", "PIK3CA", "ATM", "SMARCA4", "RBM10", "U2AF1", "SETD2", "SMAD4", "PTEN")

# Fig S8A-C: common oncogenes
onco_png(grp("ARID1A"), common_oncogenes, "ARID1A Mutants - Common Oncogenes",
         "figS8A_oncoplot_arid1a_all_common_oncogenes.png")
onco_png(grp("ARID1A", "Higher"), common_oncogenes, "ARID1A Mutants Higher Than WT Median - Common oncogenes",
         "figS8B_oncoplot_arid1a_higher_common_oncogenes.png")
onco_png(grp("ARID1A", "Lower"), common_oncogenes, "ARID1A Mutants Lower Than WT Median - Common oncogenes",
         "figS8C_oncoplot_arid1a_lower_common_oncogenes.png")

# Fig S9A-C: SWI/SNF genes
onco_png(grp("ARID1A"), swi_snf_genes, "ARID1A Mutants - SWI/SNF Genes",
         "figS9A_oncoplot_arid1a_all_swisnf.png")
onco_png(grp("ARID1A", "Higher"), swi_snf_genes, "ARID1A Mutants Higher Than WT Median - SWI/SNF Genes",
         "figS9B_oncoplot_arid1a_higher_swisnf.png")
onco_png(grp("ARID1A", "Lower"), swi_snf_genes, "ARID1A Mutants Lower Than WT Median - SWI/SNF Genes",
         "figS9C_oncoplot_arid1a_lower_swisnf.png")
