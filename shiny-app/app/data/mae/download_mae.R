tcga_codes = c("paad_tcga", "coadread_tcga","brca_tcga",
               "luad_tcga","ov_tcga","chol_tcga")
gene_list = GetGene("onco")
mols_id = c("gistic", "rna_seq_v2_mrna_median_Zscores", "rppa_Zscores")

mae_list_onco <- list()
for (tcga_code in tcga_codes) {
  mae <- cBioPortalData::cBioPortalData(
    api = cbio,
    by = "hugoGeneSymbol",
    studyId = tcga_code,
    genes = gene_list,
    molecularProfileIds = paste0(tcga_code, "_", mols_id)
  )
  mae_list_onco[[tcga_code]] = mae
}
saveRDS(mae_list_onco, "data/mae/mae_list_onco.rds")



# Genome-wide genes
tcga_codes = c("chol_tcga", "paad_tcga", "coadread_tcga","brca_tcga",
  "ov_tcga","luad_tcga")
gene_list = GetGene("gw")
mols_id = c("gistic", "rna_seq_v2_mrna_median_Zscores", "rppa_Zscores")


mae_list_onco <- list()
for (tcga_code in tcga_codes) {
  print(tcga_code)
  
  
  wk_dir <- "/Users/qixu/Library/CloudStorage/Box-Box/Kowalski Lab/Labrary/Expressed Amplified Signatures TCGA/data"
  path_cna <- sprintf("%s/%s/data_cna.txt", wk_dir, tcga_code)
  path_mrna <- sprintf("%s/%s/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt", wk_dir, tcga_code)
  path_mut <- sprintf("%s/%s/data_mutations.txt", wk_dir, tcga_code)
  
  
  # mutation
  print(sprintf("Readin mutation genes %s", tcga_code))
  long_mut <- vroom::vroom(path_mut)  %>% 
    select(gene = Hugo_Symbol, sample_id = Tumor_Sample_Barcode) %>% 
    mutate(value = 1)
  genes_mut <- unique(long_mut$gene); length(genes_mut)
  print(sprintf("Finish Readin Mutation genes %s", length(genes_mut)))
  
  
  
  # CNA
  print(sprintf("Readin CNA genes %s", tcga_code))
  cn <- vroom::vroom(path_cna) %>% select(-Entrez_Gene_Id) %>%
    dplyr::rename(gene = Hugo_Symbol)

  long_cn <- pivot_longer(cn, -gene, names_to = "sample_id", values_to = "value") %>% 
    filter(value >= 1)
  genes_amp <- unique(long_cn$gene); length(genes_amp)
  print(sprintf("Finish Readin CNA genes %s", length(genes_amp)))
  
  
  
  # MRNA
  print(sprintf("Readin MRNA genes %s", tcga_code))
  mrna <- vroom::vroom(path_mrna) %>% select(-Entrez_Gene_Id) %>%
    dplyr::rename(gene = Hugo_Symbol)
  
  long_mrna <- pivot_longer(cn, -gene, names_to = "sample_id", values_to = "value") %>% 
    filter(value >= 1)
  genes_mrna <- unique(long_mrna$gene); length(genes_mrna)
  print(sprintf("Finish Readin mRNA genes %s", length(genes_mrna)))

  
  
  # RPPA
  gene_list = gene_list
  mae_rppa <- cBioPortalData::cBioPortalData(
    api = cbio,
    by = "hugoGeneSymbol",
    studyId = tcga_code,
    genes = gene_list,
    molecularProfileIds = paste0(tcga_code, "_", "rppa_Zscores")
  )
  long_rppa <- MultiAssayExperiment::longFormat(mae_rppa) %>% as_tibble() %>% 
   select(sample_id = colname, gene = rowname, value)
  print(sprintf("finish %s", tcga_code))
  
  mae_list <- 
    list(
      "gistic" = long_cn,
      "rna_seq_v2_mrna_median_Zscores" = long_mrna,
      "rppa_Zscores" = long_rppa,
      "mutation" = long_mut
    )
  

  saveRDS(mae_list, sprintf("data/mae/mae_list_gw/%s.rds", tcga_code))
  print(sprintf("Finish saving %s", tcga_code))
}

mae_list_gw <- 
  map(tcga_codes, ~readRDS(sprintf("data/mae/mae_list_gw/%s.rds", .x)))
names(mae_list_gw) <- tcga_codes
saveRDS(mae_list_gw, "data/mae/mae_list_gw.rds")




