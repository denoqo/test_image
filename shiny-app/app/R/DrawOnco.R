


DrawOnco <- function(sig, samples, fullsample = T, top_n = 100, annorow = T, 
                     sortBand = F){
  dataeas <- sig %>% mutate(value = "Alt") 
  
  top_gene <- dplyr::count(dataeas, gene, sort = T) %>% 
    head(top_n) %>% select(gene)
  
  if (fullsample) {
    gene_sample_df <- 
      lapply(samples, function(sample){
        top_gene %>% mutate(sample_id = sample)
      }) %>% bind_rows()
    
    data <- gene_sample_df %>% 
      left_join(dataeas) 
  } else {
    data <- dataeas
  }

  mat <- data  %>% 
    spread(sample_id, value, fill = "") %>% 
    column_to_rownames("gene") 
  
  if (annorow) {
    anno <- tibble(gene = rownames(mat)) %>% 
      inner_join(annodata) %>% 
      column_to_rownames("gene")
    
    band_ele <- anno$band %>% unique()
    my_colors <- vector(length = length(band_ele))
    for (i in 1:13) {
      my_colors[i] <- rgb(runif(1), runif(1), runif(1))
    }
    names(my_colors) <- band_ele
    
    ha <- ComplexHeatmap::HeatmapAnnotation(
      df = anno, 
      which = c("row")
    )
    
    ha2 = rowAnnotation(
      foo = anno_text(anno$band),
      gp = gpar(fill = my_colors, col = "white", border = "black")
      )
    
    
  }
  
  
  ht <- ComplexHeatmap::oncoPrint(mat,
                                  right_annotation = ha2, left_annotation = ha) 
  
  if (sortBand) {
    ht <- ComplexHeatmap::oncoPrint(mat, row_split = anno$band,
                                    right_annotation = ha2, left_annotation = ha) 
  }
  return(ht)
}



