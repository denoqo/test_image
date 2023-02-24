GetGTsummary <- function(df) {
  
  if ("sample_id" %in% names(df)) {
    df = df %>% select(-sample_id)
  }
  
  table2 <- 
    gtsummary::tbl_summary(
      df,
      by = group, # split table by group
      missing = "no" # don't list missing data separately
    ) %>%
    gtsummary::add_n() %>% # add column with total number of non-missing observations
    gtsummary::add_p() %>% # test for a difference between groups
    gtsummary::modify_header(label = "**Variable**") %>% # update the column header
    gtsummary::bold_labels() %>% 
    as_gt()
  
  return(table2)
}