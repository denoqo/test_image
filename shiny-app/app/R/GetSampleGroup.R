# All samples
# Samples wtih Genes




#' @title GetSampleGroup
#' @description sample group for survival
#' @param sig 
#' @param sample_all 
#'
#' @return
#' @export
#'
#' @examples
#' sample_all <- GetSample("paad_tcga")
#' sig <- GetGeneSignature("paad_tcga", sample_all, GetGene("onco"))
#' GetSampleGroup(sample_all, sig)


GetSampleGroup <- function(
    sample_all = NULL, sig = NULL, tcga_code = NULL, gene_scope = NULL
    ) {

  if (!is.null(tcga_code) & !is.null(gene_scope)) {
    sample_all <- GetSample(tcga_code)
    sig <- GetGeneSignature(tcga_code, sample_all, GetGene(gene_scope))
  }
  sample_full <- tibble(sample_id = sample_all) %>% mutate(group = "Ctrl")
  sample_alt <- distinct(sig, sample_id) 
  sample_group <- sample_full %>% 
    mutate(group = ifelse(sample_id%in%sample_alt$sample_id, "Alt", group))
  
  return(sample_group)
}



