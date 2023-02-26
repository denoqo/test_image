

if (F) {
  tcga_code = "chol_tcga"
  scope = "gw"
  samples <- GetSample(tcga_code)
  genes <- GetGene(scope)
  if (scope == "onco") {
    sig <- GetGeneSignature(tcga_code, genes)
  } else {
    sig <- GetGeneSignatureGW(tcga_code, mae_list_gw, c("gistic", "rppa_Zscores"))
  }
 
  
  

  sample_group <- GetSampleGroup(sample_all = samples, sig = sig)
  
  
  
  
  clin_selected <- c("TMB (nonsynonymous)", "Fraction Genome Altered", "Sex")
  clin_specific <- clin %>% filter(sample_id %in% GetSample(tcga_code)) %>% 
    right_join(sample_group)
  group_clin <- clin_specific %>% 
    select(sample_id, group, all_of(clin_selected))
  
  
  asso_plot <- PlotAssociation(group_clin)
  asso_gt <- GetGTsummary(group_clin)
  surv_plot <- PerformSurvival(sample_group, clinical = clin_cleanname)
  
  
  #counts

  BarplotCount(dplyr::count(sample_group, group))
  DrawOnco(sig,samples,T)
  
  
  ngenes <- unique(sig[["gene"]])
  nsampl <- unique(sig[["sample_id"]])
  totalSamples <- nrow(sample_group)
  
  # sig
  siganno <- sig %>% left_join(annodata)
  mat_band <- dplyr::count(siganno, sample_id, band) %>% 
    spread(band, n, fill = 0) %>% 
    column_to_rownames("sample_id")
  
}


