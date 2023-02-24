GetGeneSignatureGW <- function(
    tcga_code = "chol_tcga", mae_list_gw, 
    mols_id = c("gistic",  "rppa_Zscores", "rna_seq_v2_mrna_median_Zscores", "mutation"),
    cna_cut = 1, mrna_cut = 2, rppa_cut = 2, mut_cut = 1
) {
  

  ProfileIds = mols_id
  mae_list_4 <- mae_list_gw[[tcga_code]]
  mae_list <- mae_list_4[ProfileIds]
  
  
  for (assay in names(mae_list)) {
    if (str_detect(assay, "rppa")) {
      mae_list[[assay]] <- mae_list[[assay]] %>% filter(value >= rppa_cut)
    }
    if (str_detect(assay, "mrna")) {
      mae_list[[assay]] <- mae_list[[assay]] %>% filter(value >= mrna_cut)
    }
    if (str_detect(assay, "gistic")) {
      mae_list[[assay]] <- mae_list[[assay]] %>% filter(value >= cna_cut)
    }
    if (str_detect(assay, "mutation")) {
      mae_list[[assay]] <- mae_list[[assay]] %>% filter(value >= mut_cut)
    }
  }
  mae_merge <- 
    purrr::reduce(
      mae_list, 
      dplyr::inner_join,
      by=c('gene', 'sample_id')
    )
  
  # Rename 
  mae_final <- mae_merge %>% select(gene, sample_id)
  
  return(mae_final)
  
}


