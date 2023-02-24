# Oncogenes from local
gene_path <- "data/genes/ma_list_onco_TSG.rds"
genes_onco <- readRDS(gene_path)[["oncogene"]]


# Genome-wide from cbiportal
gene_gw_coding <- readRDS("data/genes/gene_gw_coding.rds")



#' GetGene
#'
#' @param scope gw or onco
#'
#' @return gene list
#' @export
#'
#' @examples
GetGene <- function(scope) {
  query_gene <- switch(
    scope, 
    "gw" = gene_gw_coding, 
    "onco" = genes_onco
    ) 
  
  return(query_gene)
}
