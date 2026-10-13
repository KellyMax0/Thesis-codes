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

# Lollipop plot function
lollipop_png <- function(patients, file) {
  m <- subsetMaf(maf, tsb = get_maf_barcodes(patients, maf))
  png(file, width = 10, height = 7.5, units = "in", res = 300)
  lollipopPlot(maf = m, gene = "ARID1A", AACol = "HGVSp_Short", labelPos = NULL, repel = TRUE,
               showDomainLabel = TRUE, cBioPortal = FALSE)
  dev.off()
}

# Fig S7A-C: all, lower and higher expression
lollipop_png(grp("ARID1A"), "figS7A_lollipop_arid1a_all.png")
lollipop_png(grp("ARID1A", "Lower"), "figS7B_lollipop_arid1a_lower.png")
lollipop_png(grp("ARID1A", "Higher"), "figS7C_lollipop_arid1a_higher.png")
