#' @title GetGeneSignature
#'
#' @param tcga_code "paad_tcga"
#' @param sample_list GetSample(tcga_code)
#' @param gene_list GetGene("onco")
#' @param mols_id c( "gistic", "rna_seq_v2_mrna_median_Zscores", "rppa_Zscores")
#' @param cna_cut cutoff to filter cna amp (default: >=1)
#' @param mrna_cut cutoff to filter mrna high (z-score, default: >=2)
#' @param rppa_cut cutoff to filter rppa high (z-score, default: >=2)
#'
#' @return gene and sample with expected genes
#' @export
#'
#' @examples
#' 

GetGeneSignature <- function(
    tcga_code = "chol_tcga", gene_list,
    mols_id = c("gistic", "rna_seq_v2_mrna_median_Zscores"),
    cna_cut = 1, mrna_cut = 2, rppa_cut = 2
    ) {

  # mae <- cBioPortalData::cBioPortalData(
  #   api = cbio,
  #   by = "hugoGeneSymbol",
  #   studyId = tcga_code,
  #   genes = gene_list,
  #   molecularProfileIds = paste0(tcga_code, "_", mols_id)
  # )
  
  ProfileIds = paste0(tcga_code, "_", mols_id)
  mae_selected <- mae_list_onco[[tcga_code]]
  mae <- mae_selected[,,ProfileIds]
  mae_list <- as_tibble(longFormat(mae)) %>% split(., .$assay)

  
  
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
  }
  mae_merge <- 
    purrr::reduce(
      mae_list, 
      dplyr::inner_join,
      by=c('colname', 'rowname','primary')
    )
  
  # Rename 
  mae_final <- mae_merge %>% 
    select(
      sample_id = colname,
      gene = rowname
    )
  
  return(mae_final)
  
}



