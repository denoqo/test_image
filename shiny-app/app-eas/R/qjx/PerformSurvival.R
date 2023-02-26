library(survival)
library(patchwork)
library(survminer)


RunOS <- function(data, group, subset = NULL){
  
  total_lab <- data %>% pull(group) %>% unique() %>% sort()
  
  color <- scales::hue_pal()(length(total_lab))
  names(color) <- total_lab
  
  if (!is.null(subset)) {
    data = filter(data, group %in% subset)
  }
  
  fit <- survfit(
    Surv(os_months, as.numeric(substr(os_status, 1, 1))) ~ group,
    data = data
  )
  
  legend_labs <- data %>% pull(group) %>% unique() %>% sort()
  ggsurvplot(
    fit,
    data = data,
    legend.labs = legend_labs,
    palette = color,
    legend=c(0.8,0.9),
    legend.title = "",
    xlab = "Months", title="OS",
    pval = TRUE, risk.table = TRUE,
    ggtheme = theme_bw(),
    font.x = c(14, face = "bold"),
    font.y = c(14, face = "bold"),
    font.tickslab = c(14, face = "bold"),
    font.legend = c(10, face = "bold"),
    risk.table.y.text = FALSE
  )
}
RunDFS <- function(data, group, subset = NULL){
  
  total_lab <- data %>% pull(group) %>% unique() %>% sort()
  
  color <- scales::hue_pal()(length(total_lab))
  names(color) <- total_lab
  
  if (!is.null(subset)) {
    data = filter(data, group %in% subset)
  }
  
  
  fit <- survival::survfit(
    Surv(dfs_months, as.numeric(substr(dfs_status, 1, 1))) ~ group,
    data = data
  )
  
  legend_labs <- data %>% pull(group) %>% unique() %>% sort()
  ggsurvplot(
    fit,
    data = data,
    legend.labs = legend_labs,
    palette = color,
    legend=c(0.8,0.9),
    legend.title = "",
    xlab = "Months", title="DFS",
    pval = TRUE, risk.table = TRUE,
    ggtheme = theme_bw(),
    font.x = c(14, face = "bold"),
    font.y = c(14, face = "bold"),
    font.tickslab = c(14, face = "bold"),
    font.legend = c(10, face = "bold"),
    risk.table.y.text = FALSE
  )
}
RunPFS <- function(data, group, subset = NULL){
  
  total_lab <- data %>% pull(group) %>% unique() %>% sort()
  
  color <- scales::hue_pal()(length(total_lab))
  names(color) <- total_lab
  
  if (!is.null(subset)) {
    data = filter(data, group %in% subset)
  }
  
  
  fit <- survival::survfit(
    Surv(pfs_months, as.numeric(substr(pfs_status, 1, 1))) ~ group,
    data = data
  )
  
  legend_labs <- data %>% pull(group) %>% unique() %>% sort()
  ggsurvplot(
    fit,
    data = data,
    legend.labs = legend_labs,
    palette = color,
    legend=c(0.8,0.9),
    legend.title = "",
    xlab = "Months", title="PFS",
    pval = TRUE, risk.table = TRUE,
    ggtheme = theme_bw(),
    font.x = c(14, face = "bold"),
    font.y = c(14, face = "bold"),
    font.tickslab = c(14, face = "bold"),
    font.legend = c(10, face = "bold"),
    risk.table.y.text = FALSE
  )
}







# column need to be "group"
# sample colname needs to be "sample_id"
# can be used for pancan TCGA
PerformSurvival <- function(sample_group, subset_group = NULL, os_only = F, clinical = NULL){
  
  if (is.null(clinical)) {
    sur_clinical <-
      vroom::vroom("data/clinical/clin_pancan.tsv") 
  } else {
    sur_clinical <- clinical
  }
  
  sur_clinical <- sur_clinical %>% 
    janitor::clean_names() %>% 
    select(
      sample_id,
      os_status = overall_survival_status, os_months = overall_survival_months, 
      dfs_status = disease_free_status, dfs_months = disease_free_months,
      pfs_status = progression_free_status, pfs_months = progress_free_survival_months
    )
  
  clin_group <- sample_group %>% inner_join(sur_clinical)
  p_os <- RunOS(clin_group, "group", subset = subset_group)
  p_dfs <- try(RunDFS(clin_group, "group", subset =subset_group), silent = TRUE)
  p_pfs <- try(RunPFS(clin_group, "group", subset =subset_group), silent = TRUE)
  
  if(os_only) {
    p <- (p_os$plot) /
      (p_os$table) +
      plot_layout(heights = c(3,1))
    return(p)
  }
  
  if ('try-error' %in% class(p_dfs)) {
    p <- p_os$plot/p_os$table+
      plot_layout(heights = c(3,1))
  } else {
 
    p <- (p_os$plot + p_dfs$plot + p_pfs$plot ) /
      (p_os$table + p_dfs$table + p_pfs$table) +
      plot_layout(heights = c(3,1))
  }
  

  return(p)

  
}



RunEFS <- function(data, group = "group",  sur_status = "EFS", sur_value = "EFS_mnths",
                   subset = NULL){
  
  data <-  data %>%
    dplyr::rename(group = {{ group }}) %>% 
    dplyr::rename(sur_value = {{ sur_value }}) %>% 
    dplyr::rename(sur_status = {{ sur_status }}) 
  
  total_lab <- data %>% pull(group) %>% unique() %>% sort()
  
  color <- scales::hue_pal()(length(total_lab))
  names(color) <- total_lab
  
  if (!is.null(subset)) {
    data = filter(data, group %in% subset)
  }
  
  fit <- survfit(
    Surv(sur_value, as.numeric(substr(sur_status, 1, 1))) ~ group,
    data = data
  )
  
  legend_labs <- data %>% pull(group) %>% unique() %>% sort()
  ggsurvplot(
    fit,
    data = data,
    legend.labs = legend_labs,
    palette = color,
    legend=c(0.8,0.9),
    legend.title = "",
    xlab = sur_value, title=sur_status,
    pval = TRUE, risk.table = TRUE,
    ggtheme = theme_bw(),
    font.x = c(14, face = "bold"),
    font.y = c(14, face = "bold"),
    font.tickslab = c(14, face = "bold"),
    font.legend = c(10, face = "bold"),
    risk.table.y.text = FALSE
  )
}
