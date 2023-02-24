library(BiocManager)
options(repos = BiocManager::repositories())

library(tidyverse)
library(patchwork)
library(DT)
library(gtsummary)
library(ComplexHeatmap)
library(MultiAssayExperiment)
sor <- list.files("R", recursive = T, full.names = T) 
map(sor[sor != "R/global.R"], source)


cbio <- cBioPortalData::cBioPortal()
clin <- vroom::vroom("data/clinical/clin_pancan.tsv") %>% 
  dplyr::rename(sample_id = "Sample ID")
clin_cleanname <- janitor::clean_names(clin)
clin_colnames <- names(clin)
tcga_codes = c("paad_tcga", "coadread_tcga","brca_tcga",
               "luad_tcga","ov_tcga","chol_tcga")
mae_list_onco <- readRDS("data/mae/mae_list_onco.rds")
mae_list_gw <- readRDS("data/mae/mae_list_gw.rds")

# tcga_code <- "paad_tcga"
# gene_list <- GetGene("onco")
# mols_id <- c("gistic", "rna_seq_v2_mrna_median_Zscores")

# Gene annotation
band <- read_csv("data/genes/biomart_table_37_band.csv") %>% select(gene, band)
targetable <- readxl::read_excel("data/genes/geneset_gsea_Feb17_2023.xlsx", sheet = "Treatment") %>% 
  select(gene, targetable = geneset)
annodata <- left_join(band, targetable)

