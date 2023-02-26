#' @title  PlotAssociation
#'
#' @param data a tbl. sample_id, group + clinicalf feature (category or value)
#' @param allinone output all in one plot or a list of plot
#'
#' @return a plot
#' @export
#'
#' @examples
#' data_association <- readRDS("data/data_association.rds")
#' PlotAssociation(data_association, T)
PlotAssociation <- function(data, allinone = T, wrap_ncol = 3) {
  
  data <- data %>% 
    mutate(group = str_replace_all(group, "_", " ") %>% str_wrap(5))
  
  data_cat <- data %>% 
    select(sample_id, group, where(is.character)) %>% 
    gather("metric", "value", -sample_id, -group) %>% 
    mutate(metric = str_to_upper(metric)) %>% 
    mutate(value = str_wrap(value, 10))
  data_val <- data %>% 
    select(sample_id, group, where(is.numeric)) %>% 
    gather("metric", "value", -sample_id, -group) %>% 
    mutate(metric = str_to_upper(metric)) 
  
  GetBarPlot <- function(df) {
    plot_title = df$metric %>% unique() %>% 
      substr(1, 30) %>% 
      str_wrap(15)
    data <- df %>% dplyr::count(group, value) 
    
    bars <- ggplot(data, aes(x=group, y=n, fill=value)) +
      geom_bar(stat = "identity", color='black') + 
      shadowtext::geom_shadowtext(
        aes(label=n), position = position_stack(vjust=.5))+
      labs(x  = "", y = "")
    
    bars + 
      theme_light() +
      scale_fill_brewer(palette = "Set1") +
      ggtitle(plot_title) +
      guides(fill=guide_legend(title="Group"))+
      theme(legend.position = "right")
      
  }
  GetBoxPlot <- function(df) {
    plot_title = df$metric %>% unique() %>% 
      substr(1, 30) %>% 
      str_wrap(15)
    
    limits <- quantile(df$value, c(0.01, 0.95), na.rm =T)
    
    
    boxes <- ggplot(df, aes(
      x = group,
      y = value,
      color = group,
      fill = group
    )) +
      geom_jitter(width = 0.15) +
      geom_boxplot(alpha = 0.2) +
      coord_cartesian(ylim = limits) +
      labs(x  = "", y = "") 
    
    boxes +
      theme_light() +
      scale_fill_brewer(palette = "Set1") +
      theme(legend.position = "none") +
      ggtitle(plot_title) 
    
  }
  
  
  plot_cat_list <- 
    data_cat %>% 
    split(.$metric) %>% 
    lapply(GetBarPlot) 
  plot_cat <- plot_cat_list %>% 
    patchwork::wrap_plots(nrow = 1) %>% 
    try()
  
  plot_val_list <- 
    data_val %>% 
    split(.$metric) %>% 
    lapply(GetBoxPlot) 
  plot_val <- plot_val_list %>% 
    patchwork::wrap_plots(nrow = 1) %>% 
    try()
 
  
  final_list <- c(plot_cat_list, plot_val_list)
  final <- final_list %>% 
    patchwork::wrap_plots(n = wrap_ncol) 
  
  
  if (allinone == T) {
    return(final)
  } else {
    final_list
  }
}


