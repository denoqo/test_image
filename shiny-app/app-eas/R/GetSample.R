#' @title GetSample
#'
#' @param tcga_code 
#' c("paad_tcga", "ov_tcga",  "coadread_tcga", "brca_tcga", "luad_tcga", "chol_tcga")
#' 
#'
#' @return list of samples
#' @export
#'
#' @examples
#' GetSample(tcga_code)
#' 
#' 
GetSample <- function(tcga_code) {
  query_sample <- switch(
    tcga_code,
    "paad_tcga" = read_csv("data/samples/PAAD114_sample_id.csv")[["sample_id_15"]],
    "ov_tcga" = read_csv("data/samples/OV311_sample_id.csv")[["SAMPLE_ID"]],
    "coadread_tcga" = read_csv("data/samples/CRC377_sample_id.csv")[["sample_id"]],
    "brca_tcga" = read_csv("data/samples/BRCA164_sample_id.csv") %>% 
      filter(ER=="Negative" & PR=="Negative") %>% pull(pid) %>%  paste0("-01"),
    cBioPortalData::allSamples(cbio, tcga_code)[["sampleId"]]
  )
  return(query_sample)
}