#' Method for curating a variant given in variant_info
#'
#' @param variant_info named list with all the information required for curation
#' @param gnomad4 gnomad population dataset
#' @param kova_tommo kova and tommo east asian population dataset
#' @param bulk_curations bulk curations dataset for PM5 and PS1
#' @param config_values configurations
#'
#' @return list with curation information
#' \itemize{
#'   \item decision: Final curation
#'   \item rule_curations: Individual rule curations in a named list
#' }
#'
#' @export
curateVariant <- function(variant_info, gnomad4, kova_tommo, bulk_curations, config_values) {

  #Assign defaults
  variant_info$end <- NA
  variant_info$exon_info <- NA
  variant_info$spliceAI_delta_score <- NA
  variant_info$gnomAD4_ac <- NA
  variant_info$gnomAD4_an <- NA
  variant_info$gnomAD4_hc <- NA
  variant_info$gnomAD4_nfe_ac <- NA
  variant_info$gnomAD4_nfe_an <- NA
  variant_info$gnomAD4_eas_ac <- NA
  variant_info$gnomAD4_eas_an <- NA
  variant_info$kova_tommo_ac <- NA
  variant_info$kova_tommo_an <- NA
  variant_info$gnomAD4_GroupMaxFAF <- NA
  variant_info$alternate_hgvsc <- NA
  variant_info$alternate_missense <- NA
  variant_info$alternate_splice <- NA
  variant_info$var_nearby <- NA
  variant_info$spliceInfo <- NA
  variant_info$functional_domain <- NA
  variant_info$odds_pathogenicity <- NA
  variant_info$odds_ratio <- NA
  variant_info$lb_95ci <- NA
  variant_info$ub_95ci <- NA
  variant_info$curation <- list()
  variant_info$curation_log <- list()
  variant_info$combined_curation <- list()
  variant_info$apply_manual_curations <- F


  #process variant_info
  variant_info$consequence <- assignConsequence(variant_info$consequence_VEP)
  variant_info$hgvsp <- sub("[*]", "Ter", variant_info$hgvsp)

  #Exon information
  index <- which((exons$CodingStart <= variant_info$cds_start) & (exons$CodingEnd >= variant_info$cds_start))
  if (length(index) != 0) {
    exon <- exons[index, ]
    variant_info$exon_info <- exon
  }

  # SpliceAI delta score
  if (!all(is.na(variant_info$spliceAI_score)))
  {
    spliceAI_delta_score <- max(variant_info$spliceAI_score[c("DS_AG", "DS_AL", "DS_DG", "DS_DL")])
    variant_info$spliceAI_delta_score <- as.numeric(spliceAI_delta_score)
  }

  #process data and generate required info
  variant_info <- assignPopulationInfo(variant_info, gnomad4, kova_tommo, config_values)
  variant_info <- calculateOddsValues(variant_info, config_values)
  variant_info <- assignProximateCurations(variant_info, bulk_curations, config_values)

  #Curate
  variant_info <- runCurations(variant_info, config_values)
  decision <- pathogenicityDecision(variant_info$curation)

  return (list(decision=decision, rule_curations=variant_info$curation))
}
