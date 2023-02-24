# colnames
# group
# n

BarplotCount <- function(data, cols_label = NULL){
  p <- ggplot(data, aes(y = factor(group) %>% fct_rev(), x = n, fill = group, label = n)) +
    geom_bar(stat = "identity") +
    geom_text(position = position_stack(vjust = 0.5), color = "white", size = 7,
              fontface ="bold") +
    theme_light() +
    labs(x = "") +
    theme(legend.position = "none") + 
    theme(panel.grid.major = element_line(linetype = "blank"),
          panel.grid.minor = element_line(linetype = "blank")) +labs(y = NULL) + 
    theme(axis.ticks = element_line(linetype = "blank"))+
    theme( panel.border = element_blank())+
    ggtitle("Sample Numbers")
  
  if (!is.null(cols_label)) {
    p <- p + scale_fill_manual(values = cols_label)
  }
  
  return(p)
}
