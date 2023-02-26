# Genome-wide from cbiportal
gene_gw <- cBioPortalData::geneTable(cbio, pageSize = 50000)
gene_gw_coding <- gene_gw %>% filter(type == "protein-coding") %>% 
  pull(hugoGeneSymbol)
saveRDS(gene_gw_coding, "data/genes/gene_gw_coding.rds")