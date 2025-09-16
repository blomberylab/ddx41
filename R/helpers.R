#' assignPopulationInfo function
#' Search for the variant in population datasets and extract required information
#'
#' @param variant_info named list with all the information required for curation
#' @param gnomad4 gnomad population dataset
#' @param kova_tommo kova and tommo east asian population dataset
#' @param config_values configurations
#'
#' @return variant_info named list with all population information added
assignPopulationInfo <- function(variant_info, gnomad4, kova_tommo, config_values)
{
  gnomad4_sub <- subset(gnomad4, gnomad4$NM_016222.4_HGVSc == variant_info$hgvsc)

  # gnomad 4
  if (nrow(gnomad4_sub) != 0) {
    alelle_count <- gnomad4_sub$Allele.Count
    alelle_number <- gnomad4_sub$Allele.Number
    variant_info$gnomAD4_hc <- gnomad4_sub$Homozygote.Count
    variant_info$gnomAD4_GroupMaxFAF <- (gnomad4_sub$GroupMax.FAF.frequency) * 100

    variant_info$gnomAD4_nfe_ac <- gnomad4_sub$Allele.Count.European..non.Finnish.
    variant_info$gnomAD4_nfe_an <- gnomad4_sub$Allele.Number.European..non.Finnish.
    variant_info$gnomAD4_eas_ac <- gnomad4_sub$Allele.Count.East.Asian
    variant_info$gnomAD4_eas_an <- gnomad4_sub$Allele.Number.East.Asian
  } else {
    # No matching rows. Search the closest
    alelle_count <- 0

    closest_index <- which.min(abs(gnomad4$hg38.Position - variant_info$start))
    alelle_number <- mean(gnomad4$Allele.Number[closest_index])

    variant_info$gnomAD4_nfe_ac <- 0
    variant_info$gnomAD4_nfe_an <- mean(gnomad4$Allele.Number.European..non.Finnish.[closest_index])
    variant_info$gnomAD4_eas_ac <- 0
    variant_info$gnomAD4_eas_an <- mean(gnomad4$Allele.Number.East.Asian[closest_index])

    cat(paste0("gnomad4 exact match not found. Using ", paste(gnomad4$NM_016222.4_HGVSc[closest_index], collapse = ", "), "\n"))
  }

  variant_info$gnomAD4_ac <- alelle_count
  variant_info$gnomAD4_an <- alelle_number

  ### KOVA & ToMMo ----
  kova_tommo_sub <- subset(kova_tommo, kova_tommo$Transcript.Consequence == variant_info$hgvsc)
  if (nrow(kova_tommo_sub) != 0) {
    variant_info$kova_tommo_an <- kova_tommo_sub$AN
    variant_info$kova_tommo_ac <- kova_tommo_sub$AC
  }
  else
  {
    variant_info$kova_tommo_an <- stats::median(kova_tommo$AN)
    variant_info$kova_tommo_ac <- 0
  }

  return (variant_info)
}

#' calculateOddsValues function
#' Calculate ODDs ratios for PS4 and PP4
#'
#' @param variant_info named list with all the information required for curation
#' @param config_values configurations
#'
#' @return variant_info named list with odds ratios
calculateOddsValues <- function(variant_info, config_values)
{
  #Calculate odds of pathogenicity
  values <- calculate_updated_Pr_multi(variant_info$total_pp4, variant_info$single_recurrent,
                                       variant_info$single_non_recurrent, variant_info$multi_somatic)

  variant_info$odds_pathogenicity <- values$odds_path

  ### Odds ratio calculations ----
  oddsInfo <- calulateOddsRatio(variant_info, config_values$PS4_ethnicity_matching)

  #Add ratios to the reactive variable
  if (length(oddsInfo) != 0)
  {
    variant_info$odds_ratio <- oddsInfo$ODDS
    variant_info$lb_95ci <- oddsInfo$LB
    variant_info$ub_95ci <- oddsInfo$UB
  }

  #print(variant_info)
  return (variant_info)
}

#' normalize_hgvsp function
#' normalize hgvsp for partial matches
#'
#' @param hgvsp  original hgvsp
#'
#' @return normalised hgvsp
normalize_hgvsp <- function(hgvsp) {
  # Handle NA or empty strings
  if (is.na(hgvsp) | hgvsp == "" | hgvsp == "p.Met1?") {
    return(hgvsp)
  }

  # Normalize frameshift variants
  normalized <- ifelse(
    grepl("fs\\?", hgvsp, ignore.case = TRUE),
    # Normalize frameshift variants with a question mark
    gsub(
      pattern = "p\\.([A-Za-z]+)([0-9]+)[A-Za-z]*fs\\?.*",
      replacement = "p.\\1\\2fs",
      x = hgvsp,
      ignore.case = TRUE
    ),
    ifelse(
      grepl("fs", hgvsp, ignore.case = TRUE),
      gsub(
        pattern = "p\\.([A-Za-z]+)([0-9]+)[A-Za-z]*fs.*",
        replacement = "p.\\1\\2fs",
        x = hgvsp,
        ignore.case = TRUE
      ),
      # Normalize missense variants
      ifelse(
        !grepl("del|ins|dup|ter|ext", hgvsp, ignore.case = TRUE),
        stringr::str_sub(hgvsp, end = -4),
        hgvsp # Return the original for other types
      )
    )
  )

  return(normalized)
}

#' assignProximateCurations function
#' Search for the variants in proximity of this variant in curated dataset
#'
#' @param variant_info named list with all the information required for curation
#' @param bulk_curations bulk curations dataset for PM5 and PS1
#' @param config_values configurations
#'
#' @return variant_info named list with nearby bariant info added
assignProximateCurations <- function(variant_info, bulk_curations, config_values)
{
  hgvsp_query <- stringr::str_extract(variant_info$hgvsp, "p[.].*")
  hgvsc_query <- variant_info$hgvsc
  spliceInfo <- spliceJunction(hgvsc_query)

  ## Curation output ----
  bulk_curations <- bulk_curations %>% mutate(
    hgvsp_extract = stringr::str_extract(.data$hgvsp, "p[.].*"),
    hgvsp_normalized = sapply(.data$hgvsp_extract, normalize_hgvsp),
    hgvsp_normalized = unname(.data$hgvsp_normalized)
  )

  bulk_curations <- bulk_curations %>%
    mutate(hgvsc_extract = case_when(
      grepl("splice", .data$consequence) & grepl("[+-]", hgvsc) ~ stringr::str_extract(hgvsc, "c\\.\\d+([+-]\\d+)?"),
      grepl("splice", .data$consequence) ~ stringr::str_extract(hgvsc, "c\\.\\d+"),
      !grepl("splice", .data$consequence) ~ stringr::str_extract(hgvsc, "c\\.\\d+"),
      TRUE ~ NA_character_
    ))

  bulk_curations <- bulk_curations %>%
    mutate(
      spliceAI_max = sapply(.data$spliceAI_score, extract_max_spliceAI),
      spliceAI_max = unname(.data$spliceAI_max),
      spliceAI_threshold = .data$spliceAI_max - 0.1
    )

  ### curated variant subsets for PS1 and PM5 datasets ----
  #### alternate missense variants ----
  hgvsp_query_normalized <- normalize_hgvsp(hgvsp_query)
  grantham_query <- variant_info$grantham_score

  variant_info$alternate_missense <- data.frame()
  if (grepl("missense", variant_info$consequence))
  {
    alternate_missense <- bulk_curations %>%
      # remove identical missense
      filter(.data$hgvsp_extract != hgvsp_query) %>%
      # match for same codon
      filter(.data$hgvsp_normalized == hgvsp_query_normalized) %>%
      mutate(
        higher_grantham = ifelse(
          grantham_query >= .data$grantham_score, TRUE, FALSE
        )
      ) %>%
      rename(
        HGVSc = .data$hgvsc,
        HGVSp = .data$hgvsp_extract,
        Grantham = .data$grantham_score
      ) %>%
      mutate(
        Curation = case_when(
          (config_values$PS4_ethnicity_matching == "none") & (!config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_none_revel_counts,
          (config_values$PS4_ethnicity_matching == "nfe") & (!config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_nfe_revel_counts,
          (config_values$PS4_ethnicity_matching == "eas") & (!config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_eas_revel_counts,
          (config_values$PS4_ethnicity_matching == "none") & (config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_none_alphamissense_counts,
          (config_values$PS4_ethnicity_matching == "nfe") & (config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_nfe_alphamissense_counts,
          (config_values$PS4_ethnicity_matching == "eas") & (config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_eas_alphamissense_counts,
          (config_values$PS4_ethnicity_matching == "none") & (!config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_none_revel_multinomial,
          (config_values$PS4_ethnicity_matching == "nfe") & (!config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_nfe_revel_multinomial,
          (config_values$PS4_ethnicity_matching == "eas") & (!config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_eas_revel_multinomial,
          (config_values$PS4_ethnicity_matching == "none") & (config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_none_alphamissense_multinomial,
          (config_values$PS4_ethnicity_matching == "nfe") & (config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_nfe_alphamissense_multinomial,
          (config_values$PS4_ethnicity_matching == "eas") & (config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_eas_alphamissense_multinomial
        )
      ) %>%
      select(.data$HGVSc, .data$HGVSp, .data$Grantham, .data$Curation, .data$higher_grantham)
    variant_info$alternate_missense <- alternate_missense
  }

  #### same hgvsp variants ----
  # Lookup bulk_curations static table
  variant_info$alternate_hgvsc <- data.frame()
  if (grepl("missense", variant_info$consequence))
  {
    alternate_hgvsc <- bulk_curations %>%
      # Exclude variants with the same HGVSc
      filter(.data$hgvsc != hgvsc_query) %>%
      # Look for same HGVSp
      filter(.data$hgvsp == hgvsp_query) %>%
      rename(
        HGVSc = .data$hgvsc
      ) %>%
      mutate(
        Curation = case_when(
          (config_values$PS4_ethnicity_matching == "none") & (!config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_none_revel_counts,
          (config_values$PS4_ethnicity_matching == "nfe") & (!config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_nfe_revel_counts,
          (config_values$PS4_ethnicity_matching == "eas") & (!config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_eas_revel_counts,
          (config_values$PS4_ethnicity_matching == "none") & (config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_none_alphamissense_counts,
          (config_values$PS4_ethnicity_matching == "nfe") & (config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_nfe_alphamissense_counts,
          (config_values$PS4_ethnicity_matching == "eas") & (config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_eas_alphamissense_counts,
          (config_values$PS4_ethnicity_matching == "none") & (!config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_none_revel_multinomial,
          (config_values$PS4_ethnicity_matching == "nfe") & (!config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_nfe_revel_multinomial,
          (config_values$PS4_ethnicity_matching == "eas") & (!config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_eas_revel_multinomial,
          (config_values$PS4_ethnicity_matching == "none") & (config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_none_alphamissense_multinomial,
          (config_values$PS4_ethnicity_matching == "nfe") & (config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_nfe_alphamissense_multinomial,
          (config_values$PS4_ethnicity_matching == "eas") & (config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_eas_alphamissense_multinomial
        )
      ) %>%
      select(.data$HGVSc, .data$Curation) %>%
      distinct()
    variant_info$alternate_hgvsc <- alternate_hgvsc
  }

  #### variants nearby and alternate_splice ----
  variant_info$alternate_splice <- data.frame()
  if (grepl("(splice|intron)", variant_info$consequence))
  {
    eval_nearby <- bulk_curations %>%
      filter(.data$hgvsc != hgvsc_query) %>%
      mutate(
        is_same_nucleotide = .data$hgvsc_extract == spliceInfo$nucleotide,
        is_nearby_canonical = .data$hgvsc_extract %in% nearby_canonical(hgvsc_query),
        is_nearby_motif = .data$hgvsc_extract %in% nearby_motif(hgvsc_query)
      ) %>%
      mutate(
        Curation = case_when(
          (config_values$PS4_ethnicity_matching == "none") & (!config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_none_revel_counts,
          (config_values$PS4_ethnicity_matching == "nfe") & (!config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_nfe_revel_counts,
          (config_values$PS4_ethnicity_matching == "eas") & (!config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_eas_revel_counts,
          (config_values$PS4_ethnicity_matching == "none") & (config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_none_alphamissense_counts,
          (config_values$PS4_ethnicity_matching == "nfe") & (config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_nfe_alphamissense_counts,
          (config_values$PS4_ethnicity_matching == "eas") & (config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_eas_alphamissense_counts,
          (config_values$PS4_ethnicity_matching == "none") & (!config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_none_revel_multinomial,
          (config_values$PS4_ethnicity_matching == "nfe") & (!config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_nfe_revel_multinomial,
          (config_values$PS4_ethnicity_matching == "eas") & (!config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_eas_revel_multinomial,
          (config_values$PS4_ethnicity_matching == "none") & (config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_none_alphamissense_multinomial,
          (config_values$PS4_ethnicity_matching == "nfe") & (config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_nfe_alphamissense_multinomial,
          (config_values$PS4_ethnicity_matching == "eas") & (config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_eas_alphamissense_multinomial
        )
      )
    variant_info$var_nearby <- eval_nearby
    alternate_splice <- eval_nearby %>%
      filter(
        .data$is_same_nucleotide == TRUE |
          .data$is_nearby_canonical == TRUE |
          .data$is_nearby_motif == TRUE
      ) %>%
      distinct() %>%
      rename(
        HGVSc = .data$hgvsc,
        SpliceAI = .data$spliceAI_max
      )  %>%
      mutate(
        Curation = case_when(
          (config_values$PS4_ethnicity_matching == "none") & (!config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_none_revel_counts,
          (config_values$PS4_ethnicity_matching == "nfe") & (!config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_nfe_revel_counts,
          (config_values$PS4_ethnicity_matching == "eas") & (!config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_eas_revel_counts,
          (config_values$PS4_ethnicity_matching == "none") & (config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_none_alphamissense_counts,
          (config_values$PS4_ethnicity_matching == "nfe") & (config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_nfe_alphamissense_counts,
          (config_values$PS4_ethnicity_matching == "eas") & (config_values$use_alphamissense) & (config_values$PP4_method == "counts") ~ decision_eas_alphamissense_counts,
          (config_values$PS4_ethnicity_matching == "none") & (!config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_none_revel_multinomial,
          (config_values$PS4_ethnicity_matching == "nfe") & (!config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_nfe_revel_multinomial,
          (config_values$PS4_ethnicity_matching == "eas") & (!config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_eas_revel_multinomial,
          (config_values$PS4_ethnicity_matching == "none") & (config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_none_alphamissense_multinomial,
          (config_values$PS4_ethnicity_matching == "nfe") & (config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_nfe_alphamissense_multinomial,
          (config_values$PS4_ethnicity_matching == "eas") & (config_values$use_alphamissense) & (config_values$PP4_method == "multinomial") ~ decision_eas_alphamissense_multinomial
        )
      ) %>%
      select(.data$HGVSc, .data$SpliceAI, .data$Curation)
    variant_info$alternate_splice <- alternate_splice
  }

  return (variant_info)
}

#' Main method for curating a ddx41 variant, this method expects the variant_info named list
#' to have information in a specific format
#'
#' @param variant_info named list with all the information required for curation
#' @param config_values configurations
#'
#' @return variant_info with curation results
runCurations <- function(variant_info, config_values) {
  curation <- list()

  curation <- evaluatePVS1(variant_info, curation, config_values)
  curation <- evaluatePS4(variant_info, curation, config_values)
  curation <- evaluatePM2(variant_info, curation, config_values)
  curation <- evaluatePM5(variant_info, curation, config_values)
  curation <- evaluatePP1(variant_info, curation, config_values)
  curation <- evaluatePP3(variant_info, curation, config_values)
  curation <- evaluatePP4(variant_info, curation, config_values)
  curation <- evaluatePS1(variant_info, curation)
  curation <- evaluateBA1(variant_info, curation, config_values)
  curation <- evaluateBS1(variant_info, curation, config_values)
  curation <- evaluateBP2(variant_info, curation, config_values)
  curation <- evaluateBP4(variant_info, curation, config_values)
  curation <- evaluateBP7(variant_info, curation, config_values)

  # Reordering step
  desired_order <- config_values$default_rules
  curation <- curation[desired_order]

  # Assign the final reordered output
  variant_info$curation <- curation

  return(variant_info) # For batch mode, not using reactive
}


#' evaluatePVS1 function
#' PVS1 curation
#' Null variant in a gene where loss of function is a known mechanism of disease
#' @param variant_info named list with all the information required for curation
#' @param curation named list with curation results
#' @param config_values configurations
#'
#' @return named list with curations
evaluatePVS1 <- function(variant_info, curation, config_values) {

  if (variant_info$hgvsc %in% config_values$PVS1_RNA_VStrong_Variants) {
    curation[["PVS1"]] <- "VeryStrong(RNA)"
  } else if (variant_info$hgvsc %in% config_values$PVS1_VStrong_Variants) {
    curation[["PVS1"]] <- "VeryStrong"
  } else if (variant_info$hgvsc %in% config_values$PVS1_Moderate_Variants) {
      curation[["PVS1"]] <- "Moderate"
  } else if (variant_info$hgvsc %in% config_values$PVS1_NotMet_Variants) {
    curation[["PVS1"]] <- "NotMet"
  } else if (variant_info$cds_start %in% c(1, 2, 3)) # Initiation codon
  {
    curation[["PVS1"]] <- config_values$PVS1_start_loss_decision
  } else if (grepl("(nonsense|frameshift|stop[_]gained)", variant_info$consequence)) # Nonsense/frameshift/stop_gained
  {
    # NMD(nonsense mediated decay) prediction
    # get terminal aa pos
    terminal_pos <- getTerminalPos(variant_info$hgvsp, variant_info$hgvsc, variant_info$exon, config_values)
    if (terminal_pos < config_values$nmd_pos) # NMD
    {
      curation[["PVS1"]] <- "VeryStrong"
    } else if (grepl("(Ter|X)[?]$", variant_info$hgvsp)) # Check whether this is a protein extension
    {
      curation[["PVS1"]] <- "VeryStrong"
    } else { #Not predicted to undergo NMD, check whether this is multi exon
      if (grepl("[-]", variant_info$exon)) #multi exon
      {
        curation[["PVS1"]] <- "Strong"
      }
      else
      {
        curation[["PVS1"]] <- "Moderate"
      }
    }
  } else {
    # Check whether the variant is in a splice junction +/- 1 or 2 only
    # Check for splice sites
    splice_locations <- as.numeric(names(PVS1_Splicing))
    if (spliceJunctionNotSupported(variant_info$hgvsc) && (grepl("(splice|intron)", variant_info$consequence))) {
      curation[["PVS1"]] <- "NotSupported"
    } else if (!is.na(spliceJunction(variant_info$hgvsc)$spliceInfo) &&
               spliceJunction(variant_info$hgvsc)$spliceInfo == "Canonical" &&
               variant_info$cds_start %in% splice_locations) {
      vals <- PVS1_Splicing[[as.character(variant_info$cds_start)]]
      if ((!is.null(vals$exception)) && (vals$exception == variant_info$hgvsc)) {
        curation[["PVS1"]] <- "NotMet"
      } else {
        curation[["PVS1"]] <- sub("PVS1_", "", vals$effect)
      }
    } else {
      curation[["PVS1"]] <- "NotApplicable"
    }
  }

  return(curation)
}

#' evaluatePS1 function
#' PS1 curation
#' Same amino acid change as a previously established P/LP variant regardless of the nucleotide change (must be reviewed by MM-VCEP)
#' For missense variants, RNA data or splicing predictor show no imapct on splicing
#' Not applicable for splice site variants, except for the
#'   canonical (+/-2) (DO NOT apply for +2G>C variants)
#'   U2 donor motif (last 3 bases of the exon and 6 nucleotides of the intron),
#'   or the U2 acceptor motif (20 nucleotides of the intron and 1st base of the exon).
#' Splicing predictions for the variant being evaluated and the known P/LP should match.
#' @param variant_info named list with all the information required for curation
#' @param curation named list with curation results
#'
#' @return named list with curations
#' @import dplyr
#' @import tidyr
evaluatePS1 <- function(variant_info, curation)
{
  hgvsc_query <- variant_info$hgvsc
  alternate_hgvsc <- variant_info$alternate_hgvsc
  alternate_splice <- variant_info$alternate_splice
  eval_PVS1 <- curation[["PVS1"]]
  eval_PP3 <- curation[["PP3"]]

  # For Missense Variants
  if (grepl("missense", variant_info$consequence))
  {
    ps1_status <- alternate_hgvsc %>%
      summarise(
        Pathogenic_count = sum(.data$Curation == "Pathogenic"),
        LikelyPathogenic_count = sum(.data$Curation == "Likely Pathogenic")
      ) %>%
      mutate(
        PS1_status = case_when(
          Pathogenic_count > 0 ~ "Strong",
          LikelyPathogenic_count > 0 ~ "Moderate",
          TRUE ~ "NotMet"
        )
      )
  }
  # For Splice Site Variants
  if (grepl(("splice|intron"), variant_info$consequence))
  {
    splice_result <- evaluate_splice(hgvsc_query, variant_info, eval_PVS1, eval_PP3)
  }

  if (grepl("missense", variant_info$consequence) && !grepl(("splice|intron"), variant_info$consequence))
  {
    curation[["PS1"]] <- ps1_status$PS1_status
  }
  else if (grepl(("splice|intron"), variant_info$consequence) && !grepl("missense", variant_info$consequence))
  {
    curation[["PS1"]] <- splice_result
  }
  else if (grepl(("splice|intron"), variant_info$consequence) && grepl("missense", variant_info$consequence))
  {
    #missense , Strong|Moderate|NotMet
    #splice|intron, NotApplicable|NotSupported|Strong|Moderate|Supporting|NotMet
    missense_outcome <- factor(ps1_status$PS1_status, levels = c("NotSupported", "NotApplicable", "NotMet", "Supporting", "Moderate", "Strong"), ordered = TRUE)
    splice_outcome <- factor(splice_result, levels = c("NotSupported", "NotApplicable", "NotMet", "Supporting", "Moderate", "Strong"), ordered = TRUE)

    if (missense_outcome >= splice_outcome)
    {
      curation[["PS1"]] <- ps1_status$PS1_status
    }
    else
    {
      curation[["PS1"]] <- splice_result
    }
  }
  else
  {
    curation[["PS1"]] <- "NotApplicable"
  }

  return(curation)
}

#' evaluatePS4 function
#' PS4 curation
#' Prevalence of the variant in affected individuals is significantly increased compared
#' with the prevalence in controls
#' @param variant_info named list with all the information required for curation
#' @param curation named list with curation results
#' @param config_values configurations
#'
#' @return named list with curations
evaluatePS4 <- function(variant_info, curation, config_values) {

  #If use known PS4 recurrent variants is active and the variant is in the lists
  if ((config_values$use_known_variants_ps4) && (variant_info$hgvsc %in% c(config_values$PS4_Moderate_Variants, config_values$PS4_Supporting_Variants)))
  {
    if (variant_info$hgvsc %in% config_values$PS4_Moderate_Variants) # Hard coded PS4_Moderate
    {
      curation[["PS4"]] <- "Moderate"
    } else if (variant_info$hgvsc %in% config_values$PS4_Supporting_Variants) # Hard coded PS4_Supporting
    {
      curation[["PS4"]] <- "Supporting"
    }
  }
  else #If use known PS4 recurrent variants is not selected or is selected and the variant is not in the lists
  {
    if (variant_info$literature_gcount == 0) {
      curation[["PS4"]] <- "NotMet"
    } else {
      ethnicity <- "None"
      if (config_values$PS4_ethnicity_matching == "eas")
        ethnicity <- "East Asian"
      else if (config_values$PS4_ethnicity_matching == "nfe")
        ethnicity <- "European non-Finnish"

      gnomad_ac <- variant_info$gnomAD4_ac
      gnomad_an <- variant_info$gnomAD4_an

      if (gnomad_ac > (gnomad_an/2)) # check whether this is a homozygous variant
      {
        curation[["PS4"]] <- "NotApplicable"
      }
      else
      {
        odds_ratio <- variant_info$odds_ratio
        lb_95ci <- variant_info$lb_95ci
        ub_95ci <- variant_info$ub_95ci

        if (lb_95ci >= config_values$PS4_Very_Strong) {
          curation[["PS4"]] <- "VeryStrong"
        } else if (lb_95ci >= config_values$PS4_Strong) {
          curation[["PS4"]] <- "Strong"
        } else if (lb_95ci >= config_values$PS4_Moderate) {
          curation[["PS4"]] <- "Moderate"
        } else if (lb_95ci >= config_values$PS4_Supporting) {
          curation[["PS4"]] <- "Supporting"
        } else {
          curation[["PS4"]] <- "NotMet"
        }
      }
    }
  }

  return(curation)
}


#' evaluatePM2 function
#' PM2 curation
#' GroupMax FAF applied for now
#' Possibly later change to not applicable VS absent/rare in population databases
#' @param variant_info named list with all the information required for curation
#' @param curation named list with curation results
#' @param config_values configurations
#'
#' @return named list with curations
evaluatePM2 <- function(variant_info, curation, config_values) {

  gnomad_groupmax <- variant_info$gnomAD4_GroupMaxFAF
  eval_PS4 <- curation[["PS4"]]

  if (!is.na(eval_PS4) && eval_PS4 != "NotMet") {
    curation[["PM2"]] <- "NotApplicable"
  } else if ((is.na(gnomad_groupmax)) && (variant_info$gnomAD4_ac == 0)) {
    curation[["PM2"]] <- "Supporting"
  } else if ((is.na(gnomad_groupmax)) && (variant_info$gnomAD4_ac != 0)) {
    curation[["PM2"]] <- "NotMet"
  } else if (gnomad_groupmax <= config_values$PM2_AF) {
    curation[["PM2"]] <- "Supporting"
  } else {
    curation[["PM2"]] <- "NotMet"
  }

  return(curation)
}


#' evaluatePM5 function
#' PM5 curation
#' Missense change at an amino acid residue where a different missense change which has
#'   been determined to be pathogenic before.
#' Nonsense/frameshift variants throughout the gene in transcript NM_016222.4
#' @param variant_info named list with all the information required for curation
#' @param curation named list with curation results
#' @param config_values configurations
#'
#' @return named list with curations
#'
#' @import dplyr
#' @import tidyr
evaluatePM5 <- function(variant_info, curation, config_values) {

  alternate_missense <- variant_info$alternate_missense
  grantham_query <- variant_info$grantham_score

  # PM5_Supporting for null variants
  if (grepl("(frameshift)", variant_info$consequence))
  {
    curation[["PM5"]] <- config_values$PM5_frameshift_decision
  }
  else if (grepl("(stop[_]gained)", variant_info$consequence))
  {
    curation[["PM5"]] <- config_values$PM5_stop_gained_decision
  }
  # Further logic branches if missense variant
  else if (grepl("missense", variant_info$consequence))
  {
    pm5_status <- alternate_missense %>%
      filter(.data$higher_grantham == TRUE) %>%
      summarise(
        Pathogenic_count = sum(.data$Curation == "Pathogenic"),
        LikelyPathogenic_count = sum(.data$Curation == "Likely Pathogenic")
      ) %>%
      mutate(
        PM5_status = case_when(
          Pathogenic_count > 1 ~ config_values$PM5_pathogenic_morethan1_decision,
          Pathogenic_count == 1 ~ config_values$PM5_pathogenic_equal1_decision,
          LikelyPathogenic_count > 0 ~ config_values$PM5_likelyPathogenic_decision,
          TRUE ~ "NotMet"
        )
      )

    curation[["PM5"]] <- if (nrow(alternate_missense) > 0) pm5_status$PM5_status else "NotMet"
  }
  else
  {
    curation[["PM5"]] <- "NotApplicable"
  }

  return(curation)
}


#' evaluatePP1 function
#' PP1 curation
#' Co-segregation with disease in multiple affected family members
#' @param variant_info named list with all the information required for curation
#' @param curation named list with curation results
#' @param config_values configurations
#'
#' @return named list with curations
evaluatePP1 <- function(variant_info, curation, config_values) {

  if (variant_info$hgvsc %in% config_values$PP1_Strong_Variants) {
    curation[["PP1"]] <- "Strong"
  } else if (variant_info$hgvsc %in% config_values$PP1_Supporting_variants) {
    curation[["PP1"]] <- "Supporting"
  } else {
    curation[["PP1"]] <- "NotMet"
  }

  return(curation)
}


#' evaluatePP3 function
#' PP3 curation
#' Multiple lines of computational evidence support a deleterious effect on the gene or gene product
#' @param variant_info named list with all the information required for curation
#' @param curation named list with curation results
#' @param config_values configurations
#'
#' @return named list with curations
evaluatePP3 <- function(variant_info, curation, config_values) {

  spliceAI_delta_score <- variant_info$spliceAI_delta_score

  if (grepl("missense", variant_info$consequence)) {
    if (!config_values$use_alphamissense) #use REVEL
    {
      if (((!is.na(variant_info$revel_score)) && (variant_info$revel_score >= config_values$PP3_REVEL)) &&
          ((!is.na(spliceAI_delta_score)) && (spliceAI_delta_score >= config_values$PP3_SpliceAI))) {
        curation[["PP3"]] <- "Supporting"
      } else if ((!is.na(variant_info$revel_score)) && (variant_info$revel_score >= config_values$PP3_REVEL)) {
        curation[["PP3"]] <- "Supporting"
      } else if ((!is.na(spliceAI_delta_score)) && (spliceAI_delta_score >= config_values$PP3_SpliceAI)) {
        curation[["PP3"]] <- "Supporting"
      } else {
        curation[["PP3"]] <- "NotMet"
      }
    }
    else #use AlphaMissense
    {
      if (((!is.na(variant_info$alphamissense_class)) && (variant_info$alphamissense_class == "likely_pathogenic")) &&
          ((!is.na(spliceAI_delta_score)) && (spliceAI_delta_score >= config_values$PP3_SpliceAI))) {
        curation[["PP3"]] <- "Supporting"
      } else if ((!is.na(variant_info$alphamissense_class)) && (variant_info$alphamissense_class == "likely_pathogenic")) {
        curation[["PP3"]] <- "Supporting"
      } else if ((!is.na(spliceAI_delta_score)) && (spliceAI_delta_score >= config_values$PP3_SpliceAI)) {
        curation[["PP3"]] <- "Supporting"
      } else {
        curation[["PP3"]] <- "NotMet"
      }
    }
  } else if ((grepl("(synonymous|intron)", variant_info$consequence)) && !canonicalSpliceSite(variant_info$hgvsc)) {
    if ((!is.na(spliceAI_delta_score)) && (spliceAI_delta_score >= config_values$PP3_SpliceAI)) {
      curation[["PP3"]] <- "Supporting"
    } else {
      curation[["PP3"]] <- "NotMet"
    }
  } else {
    curation[["PP3"]] <- "NotApplicable"
  }

  return(curation)
}

#' evaluatePP4 function
#' PP4 curation
#' Patient's phenotype or family history is highly specific for a disease with a single
#'  genetic etiology
#' @param variant_info named list with all the information required for curation
#' @param curation named list with curation results
#' @param config_values configurations
#'
#' @return named list with curations
evaluatePP4 <- function(variant_info, curation, config_values) {

  #check user choice
  if (config_values$PP4_method == "counts")
  {
    # Check whether we have somatic second hit evidence
    if (length(variant_info$literature_somatic) == 0) {
      curation[["PP4"]] <- "NotMet"
    } else {
      strong_index <- integer(0)
      for (var in config_values$PP4_Strong_Variants)
      {
        strong_index <- c(strong_index, which(grepl(var, variant_info$literature_somatic, fixed = T)))
      }

      var_strong <- variant_info$literature_somatic[strong_index]
      if (length(var_strong) != 0) {
        curation[["PP4"]] <- "Strong"
      } else {
        var_moderate <- variant_info$literature_somatic
        curation[["PP4"]] <- "Moderate"
      }
    }
  }
  else #multinomial prosterior probabilities
  {
    odds_path <- variant_info$odds_pathogenicity

    if ((!is.na(odds_path)) && (odds_path >= config_values$PP4_Very_Strong)) {
      curation[["PP4"]] <- "VeryStrong"
    } else if ((!is.na(odds_path)) && (odds_path >= config_values$PP4_Strong)) {
      curation[["PP4"]] <- "Strong"
    } else if ((!is.na(odds_path)) && (odds_path >= config_values$PP4_Moderate)) {
      curation[["PP4"]] <- "Moderate"
    } else if ((!is.na(odds_path)) && (odds_path >= config_values$PP4_Supporting)) {
      curation[["PP4"]] <- "Supporting"
    } else {
      curation[["PP4"]] <- "NotMet"
    }
  }

  return(curation)
}

#' evaluateBA1 function
#' BA1 curation
#' Allele frequency is > 5 percent in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium.
#' @param variant_info named list with all the information required for curation
#' @param curation named list with curation results
#' @param config_values configurations
#'
#' @return named list with curations
evaluateBA1 <- function(variant_info, curation, config_values) {

  gnomad_groupmax <- variant_info$gnomAD4_GroupMaxFAF

  if (is.na(gnomad_groupmax)) {
    curation[["BA1"]] <- "NotMet"
  } else if (gnomad_groupmax >= config_values$BA1_AF) {
    curation[["BA1"]] <- "Stand-alone_Benign"
  } else {
    curation[["BA1"]] <- "NotMet"
  }

  return(curation)
}


#' evaluateBS1 function
#' BS1 curation
#' Allele frequency is greater than expected for disorder.
#' @param variant_info named list with all the information required for curation
#' @param curation named list with curation results
#' @param config_values configurations
#'
#' @return named list with curations
evaluateBS1 <- function(variant_info, curation, config_values) {

  gnomad_groupmax <- variant_info$gnomAD4_GroupMaxFAF

  if (is.na(gnomad_groupmax)) {
    curation[["BS1"]] <- "NotMet"
  } else if (gnomad_groupmax > config_values$BA1_AF) {
    curation[["BS1"]] <- "NotMet"
  } else if (gnomad_groupmax >= config_values$BS1_AF) {
    curation[["BS1"]] <- "Strong"
  } else {
    curation[["BS1"]] <- "NotMet"
  }

  return(curation)
}

#' evaluateBP2 function
#' BP2 curation
#' Homozygosity in an unaffected individual
#' Not used when PP4 strong
#' @param variant_info named list with all the information required for curation
#' @param curation named list with curation results
#' @param config_values configurations
#'
#' @return named list with curations
evaluateBP2 <- function(variant_info, curation, config_values) {

  #Find whether the variant has been found with PP4 strong variants
  for (var in config_values$PP4_Strong_Variants)
  {
    strong_index <- c(strong_index, which(grepl(var, variant_info$literature_somatic$HGVSp_somatic, fixed = T)))
  }
  var_strong <- variant_info$literature_somatic[strong_index, ]

  if (nrow(var_strong) != 0)
  {
    curation[["BP2"]] <- "NotApplicable"
  }
  else
  {
    #Get number of Homozygous individuals
    gnomad_hc <- variant_info$gnomAD4_hc
    if (!is.na(gnomad_hc) && gnomad_hc>0)
    {
      curation[["BP2"]] <- "Supporting"
    }
    else
    {
      curation[["BP2"]] <- "NotMet"
    }
  }

  return(curation)
}

#' evaluateBP4 function
#' BP4 curation
#' Multiple lines of computational evidence suggest no impact on gene or gene product.
#' @param variant_info named list with all the information required for curation
#' @param curation named list with curation results
#' @param config_values configurations
#'
#' @return named list with curations
evaluateBP4 <- function(variant_info, curation, config_values) {

  spliceAI_delta_score <- variant_info$spliceAI_delta_score

  if (grepl("missense", variant_info$consequence)) {
    if (!config_values$use_alphamissense) #use REVEL
    {
      if ((!is.na(variant_info$revel_score)) && (variant_info$revel_score <= config_values$BP4_REVEL) &&
          (!is.na(spliceAI_delta_score)) && (spliceAI_delta_score <= config_values$BP4_SpliceAI)) {
        curation[["BP4"]] <- "Supporting"
      } else if ((!is.na(variant_info$revel_score)) && (variant_info$revel_score <= config_values$BP4_REVEL)) {
        curation[["BP4"]] <- "NotMet"
      } else if ((!is.na(spliceAI_delta_score)) && (spliceAI_delta_score <= config_values$BP4_SpliceAI)) {
        curation[["BP4"]] <- "NotMet"
      } else {
        curation[["BP4"]] <- "NotMet"
      }
    }
    else #use AlphaMissense
    {
      if (((!is.na(variant_info$alphamissense_class)) && (variant_info$alphamissense_class == "likely_benign")) &&
          (!is.na(spliceAI_delta_score)) && (spliceAI_delta_score <= config_values$BP4_SpliceAI)) {
        curation[["BP4"]] <- "Supporting"
      } else if ((!is.na(variant_info$alphamissense_class)) && (variant_info$alphamissense_class == "likely_benign")) {
        curation[["BP4"]] <- "NotMet"
      } else if ((!is.na(spliceAI_delta_score)) && (spliceAI_delta_score <= config_values$BP4_SpliceAI)) {
        curation[["BP4"]] <- "NotMet"
      } else {
        curation[["BP4"]] <- "NotMet"
      }
    }
  } else if ((grepl("(synonymous|intron)", variant_info$consequence)) && !canonicalSpliceSite(variant_info$hgvsc)) {
    if ((!is.na(spliceAI_delta_score)) && (spliceAI_delta_score <= config_values$BP4_SpliceAI)) {
      curation[["BP4"]] <- "Supporting"
    } else {
      curation[["BP4"]] <- "NotMet"
    }
  } else {
    curation[["BP4"]] <- "NotApplicable"
  }

  return(curation)
}

#' evaluateBP7 function
#' BP7 curation
#' A synonymous (silent) variant for which splicing prediction algorithms predict no impact to the splice consensus
#'      sequence nor the creation of a new splice site AND the nucleotide is not highly conserved.
#' @param variant_info named list with all the information required for curation
#' @param curation named list with curation results
#' @param config_values configurations
#'
#' @return named list with curations
evaluateBP7 <- function(variant_info, curation, config_values) {

  spliceAI_delta_score <- NA
  if (!all(is.na(variant_info$spliceAI_score))) {
    spliceAI_delta_score <- max(variant_info$spliceAI_score[c("DS_AG", "DS_AL", "DS_DG", "DS_DL")])
  }

  if ((grepl("synonymous", variant_info$consequence) && checkSynonymousSite(variant_info$cds_start, variant_info$exon_info)) ||
      (grepl("intron", variant_info$consequence) && checkIntronicSite(variant_info$hgvsc))) {
    if ((!is.na(spliceAI_delta_score)) && (spliceAI_delta_score <= config_values$BP7_SpliceAI)) {
      curation[["BP7"]] <- "Supporting"
    } else {
      curation[["BP7"]] <- "NotMet"
    }
  } else {
    curation[["BP7"]] <- "NotApplicable"
  }

  return(curation)
}

#' pathogenicityDecision function
#' Pathogenicity decision rules
#' @param lst_curation list of curations, one list item per rule
#'
#' @return decision
pathogenicityDecision <- function(lst_curation) {
  rules <- names(lst_curation)

  ### check evidence numbers
  num_pathogenic_veryStrong <- 0
  num_pathogenic_strong <- 0
  num_pathogenic_moderate <- 0
  num_pathogenic_supporting <- 0
  num_benign_stand_alone <- 0
  num_benign_strong <- 0
  num_benign_moderate <- 0
  num_benign_supporting <- 0

  for (i in 1:length(lst_curation))
  {
    if (grepl("VeryStrong", lst_curation[[i]]) && grepl("^P", rules[i])) {
      num_pathogenic_veryStrong <- num_pathogenic_veryStrong + 1
    } else if ((lst_curation[[i]] == "Strong") && grepl("^P", rules[i])) {
      num_pathogenic_strong <- num_pathogenic_strong + 1
    } else if ((lst_curation[[i]] == "Moderate") && grepl("^P", rules[i])) {
      num_pathogenic_moderate <- num_pathogenic_moderate + 1
    } else if ((lst_curation[[i]] == "Supporting") && grepl("^P", rules[i])) {
      num_pathogenic_supporting <- num_pathogenic_supporting + 1
    } else if ((lst_curation[[i]] == "Stand-alone_Benign") && (rules[i] == "BA1")) {
      num_benign_stand_alone <- num_benign_stand_alone + 1
    } else if ((lst_curation[[i]] == "Strong") && grepl("^B", rules[i])) {
      num_benign_strong <- num_benign_strong + 1
    } else if ((lst_curation[[i]] == "Moderate") && grepl("^B", rules[i])) {
      num_benign_moderate <- num_benign_moderate + 1
    } else if ((lst_curation[[i]] == "Supporting") && grepl("^B", rules[i])) {
      num_benign_supporting <- num_benign_supporting + 1
    }
  }

  ### check for conflicting evidence and adjust numbers
  pathogenic_score <- num_pathogenic_veryStrong*8 + num_pathogenic_strong*4 + num_pathogenic_moderate*2 + num_pathogenic_supporting*1
  benign_score <- num_benign_stand_alone*8 + num_benign_strong*4 + num_benign_moderate*2 + num_benign_supporting*1
  conflicting <- F
  score <- pathogenic_score - benign_score

  if ((pathogenic_score != 0) && (benign_score != 0))
    conflicting <- T

  ### final decision
  decision <- "Uncertain"
  if (!conflicting)
  {
    if ((num_pathogenic_veryStrong >= 2) ||
        ((num_pathogenic_veryStrong >= 1) && ((num_pathogenic_strong >= 1) ||
          ((num_pathogenic_moderate + num_pathogenic_supporting) >= 2))) ||
        (num_pathogenic_strong >= 2) ||
        ((num_pathogenic_strong >= 1) &&
         ((num_pathogenic_moderate >= 3) ||
          ((num_pathogenic_moderate >= 2) && (num_pathogenic_supporting >= 2)) ||
          ((num_pathogenic_moderate >= 1) && (num_pathogenic_supporting >= 4)))))
    {
      decision <- "Pathogenic"
    } else if (((num_pathogenic_veryStrong >= 1) && (num_pathogenic_moderate >= 1)) ||
               ((num_pathogenic_strong >= 1) && (num_pathogenic_moderate >= 1)) ||
               ((num_pathogenic_strong >= 1) && (num_pathogenic_supporting >= 2)) ||
               (num_pathogenic_moderate >= 3) ||
               ((num_pathogenic_moderate >= 2) && (num_pathogenic_supporting >= 2)) ||
               ((num_pathogenic_moderate >= 1) && (num_pathogenic_supporting >= 4))) {
      decision <- "Likely Pathogenic"
    }
    # Likely pathogenic exception
    else if ((lst_curation$PVS1 == "VeryStrong") && (lst_curation$PM2 == "Supporting"))
    {
      decision <- "Likely Pathogenic"
    }
    else if ((num_benign_stand_alone == 1) || (num_benign_strong >= 2)) # Benign
    {
      decision <- "Benign"
    }
    else if ((lst_curation$BS1 == "Strong") ||
             ((num_benign_strong >= 1) && (num_benign_supporting >= 1)) ||
             (num_benign_supporting >= 2)) # Likely Benign
    {
      decision <- "Likely Benign"
    }
  }
  else # conflicting evidence. Decide based on the points
  {
    if (score >= 10)
    {
      decision <- "Pathogenic"
    } else if (score >= 6)
    {
      decision <- "Likely Pathogenic"
    }
    else if (score >= 0)
    {
      decision <- "Uncertain"
    }
    else if (score >= -6) # Likely Benign
    {
      decision <- "Likely Benign"
    }
    else #(score >= 7) # Benign
    {
      decision <- "Benign"
    }
  }

  return(decision)
}


#' calulateOddsRatio function
#' Calculate Odds Ratio and CIs
#' @param variant_info named list with all the information required for curation
#' @param ethnicity ethnicity choice
#'
#' @return odds ratio, lower bound and upper bound of 95 percent confidence intervals
calulateOddsRatio <- function(variant_info, ethnicity) {

  if (ethnicity == "none")
  {
    gnomad_ac <- variant_info$gnomAD4_ac
    gnomad_an <- variant_info$gnomAD4_an
  }
  else if (ethnicity == "nfe")
  {
    gnomad_ac <- variant_info$gnomAD4_nfe_ac
    gnomad_an <- variant_info$gnomAD4_nfe_an
  }
  else if (ethnicity == "eas")
  {
    gnomad_ac <- sum(variant_info$gnomAD4_eas_ac, variant_info$kova_tommo_ac, na.rm=T)
    gnomad_an <- sum(variant_info$gnomAD4_eas_an, variant_info$kova_tommo_an, na.rm=T)
  }

  cases_with <- variant_info$literature_gcount
  cases_without <- variant_info$literature_cases - cases_with

  if (gnomad_ac < (gnomad_an/2))
  {
    #  cases_without <- sum(reference_info$Patients, na.rm = T) - cases_with
    controls_with <- gnomad_ac
    controls_without <- (gnomad_an / 2) - controls_with

    if ((cases_with == 0) || (controls_with == 0))
    {
      cases_with <- cases_with + 0.5
      cases_without <- cases_without + 0.5
      controls_with <- controls_with + 0.5
      controls_without <- controls_without + 0.5
    }
    #odds ratio
    odds_ratio <- (cases_with / cases_without) / (controls_with / controls_without)

    # standard error
    se_log_or <- sqrt((1/cases_with) + (1/cases_without) + (1/controls_with) + (1/controls_without))

    # Calculate the 95% confidence interval.
    lb_95ci <- exp(log(odds_ratio) - 1.96 * se_log_or)
    ub_95ci <- exp(log(odds_ratio) + 1.96 * se_log_or)
    log_message <- paste0("a:", cases_with, " b: ", cases_without, " c: ", controls_with, " d: ", controls_without)
    #cat(paste0("Odds Ratio parameters ", log_message, "\n"))

    return(list(ODDS = round(odds_ratio, digits = 2), LB = round(lb_95ci, digits = 2), UB = round(ub_95ci, digits = 2), log = log_message))
  }

  return(list())
}

#' getTerminalPos function
#' Get terminal position of a hgvsp and convert to coding position
#' If not found check hgvsc
#' @param hgvsp hgvsp
#' @param hgvsc hgvsc
#' @param exon exon
#' @param config_values configurations
#'
#' @return terminal position
getTerminalPos <- function(hgvsp, hgvsc, exon, config_values) {
  #If the variant involves multiple exons use hgvsc
  if (grepl("[-]", exon))
  {
    if (grepl("c[.][*]", hgvsc)) #3_prime_UTR_variant
      return (config_values$full_coding_length)
    else if (grepl("c[.][-]", hgvsc)) #3_prime_UTR_variant
      return (0)
    else
    {
      cds_match <- regexec("c\\.(\\d+)([+]|[_]|[-])?(\\d+)?(\\D+)?", hgvsc)
      cds_start_pos <- as.integer(regmatches(hgvsc, cds_match)[[1]][2])
      cds_stop_pos <- as.integer(regmatches(hgvsc, cds_match)[[1]][4])

      if ((is.na(cds_stop_pos)) || (cds_stop_pos < cds_start_pos)) {
        return(cds_start_pos)
      } else {
        return(cds_stop_pos)
      }
    }
  }
  else  #else use hgvsp first
  {
    hgvsp <- stringr::str_extract(hgvsp, "p[.].*")
    fs_match <- regexec("p\\.(\\D+)(\\d+)(\\D+)?fs(\\*|X|Ter)?(\\d+|\\?)?", hgvsp)
    stopgain_match <- regexec("p\\.(\\D+)(\\d+)(\\D+)?(\\*|X|Ter)", hgvsp)
    stopgain2_match <- regexec("p\\.(\\D+)(\\d+)[_](\\D+)(\\d+)(\\D+)(\\*|X|Ter)", hgvsp)
    ext_match <- regexec("p\\.(\\D+)(\\d+)(\\D+)?ext(\\D+)?(\\d+|\\?)", hgvsp)

    fs_start_pos <- as.integer(regmatches(hgvsp, fs_match)[[1]][3])
    fs_stop_pos <- as.integer(regmatches(hgvsp, fs_match)[[1]][6])

    stopgain_pos <- as.integer(regmatches(hgvsp, stopgain_match)[[1]][3])
    stopgain_delins_pos <- as.integer(regmatches(hgvsp, stopgain2_match)[[1]][5])

    ext_start_pos <- as.integer(regmatches(hgvsp, ext_match)[[1]][3])
    ext_stop_pos <- as.integer(regmatches(hgvsp, ext_match)[[1]][6])
    if ((!is.na(fs_start_pos)) && (!is.na(fs_stop_pos)))
      return ((fs_start_pos + fs_stop_pos)*3)
    else if ((!is.na(fs_start_pos)) && (is.na(fs_stop_pos)))
      return ((fs_start_pos)*3)
    else if (!is.na(stopgain_pos))
      return ((stopgain_pos)*3)
    else if (!is.na(stopgain_delins_pos))
      return ((stopgain_delins_pos)*3)
    else if ((hgvsp == "") || (hgvsp == "p.?"))
    {
      if (grepl("c[.][*]", hgvsc)) #3_prime_UTR_variant
        return (config_values$full_coding_length)
      else if (grepl("c[.][-]", hgvsc)) #3_prime_UTR_variant
        return (0)
      else
      {
        cds_match <- regexec("c\\.(\\d+)([+]?)(\\d+)?(\\D+)?", hgvsc)
        cds_start_pos <- as.integer(regmatches(hgvsc, cds_match)[[1]][2])
        cds_stop_pos <- as.integer(regmatches(hgvsc, cds_match)[[1]][4])

        if ((is.na(cds_stop_pos)) || (cds_stop_pos < cds_start_pos)) {
          return(cds_start_pos)
        } else {
          return(cds_stop_pos)
        }
      }
    } else if ((!is.na(ext_start_pos)) && (!is.na(ext_stop_pos))) {
      return((ext_start_pos + ext_stop_pos) * 3)
    } else if ((!is.na(ext_start_pos)) && (is.na(ext_stop_pos))) {
      return((ext_start_pos) * 3)
    } else # Something wrong?
    {
      return(NA)
    }
  }
}

#' spliceJunction function
#' Check whether the variant is in a splice junction and returns if "Canonical", "Motif" or "Exception"
#' @param  hgvsc hgvsc
#'
#' @return splice junction info
spliceJunction <- function(hgvsc) {
  spliceInfo <- data.frame()

  # Extract position and offset from HGVS notation
  position_match <- stringr::str_match(hgvsc, "c[.]([0-9]+)([+|-]?[0-9]*)[A-Z]>[A-Z]")
  position <- as.numeric(position_match[2]) # Base position (e.g., 27)
  offset <- position_match[3] # Offset (e.g., +2, -3)
  nucleotide <- substr(hgvsc, 1, nchar(hgvsc) - 3)

  PVS1_Splicing_df <- do.call(rbind, lapply(names(PVS1_Splicing), function(pos) {
    info <- PVS1_Splicing[[pos]]
    data.frame(
      position = as.integer(pos),
      site = info$site,
      role = info$role,
      effect = info$effect,
      exception = ifelse(!is.null(info$exception), info$exception, NA)
    )
  }))

  last_3_bases <- PVS1_Splicing_df %>%
    filter(.data$role == "donor") %>%
    select(position) %>%
    reframe(position = c(position, position - 1, position - 2)) # Exon/intron boundary (donor)

  if (grepl("\\+2T>C", position_match[1])) {
    spliceInfo <- "Exception"
  } else if (as.numeric(offset) %in% c(-2, -1, 1, 2)) {
    spliceInfo <- "Canonical"
  } else if (
    # U2 donor motif (last 3 bases of the exon)
    position %in% last_3_bases$position & offset == "" ||

    # U2 donor motif (6 nucleotides of the intron)
    as.numeric(offset) %in% c(3:6) ||

    # U2 acceptor motif (first base of the exon)
    position %in% PVS1_Splicing_df$position & offset == "" ||

    # U2 acceptor motif (20 nucleotides of the intron)
    as.numeric(offset) %in% c(-20:-3)
  ) {
    spliceInfo <- "U2 motif"
  } else if (is.na(as.numeric(offset)) || (as.numeric(offset) > 6 || as.numeric(offset) < -20)) {
    spliceInfo <- "Outside"
  } else {
    spliceInfo <- NA_character_
  }

  return(list(position = position, offset = offset, spliceInfo = spliceInfo, nucleotide = nucleotide))
}

#' spliceJunctionNotSupported function
#' Check whether the variant is in a splice junction +/- 1 or 2 only but del, ins, or delins
#' @param hgvsc hgvsc
#'
#' @return vartiant in a splice junction or not
spliceJunctionNotSupported <- function(hgvsc) {
  if (grepl("c\\..*(del|ins|delins|dup)", hgvsc)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' nearby_canonical function
#' Get the nearby canonical nucleotide positions
#' @param hgvsc hgvsc
#'
#' @return nearby canonical nucleotide positions
nearby_canonical <- function(hgvsc) {
  spliceInfo <- spliceJunction(hgvsc)
  if (is.na(spliceInfo$offset))
    return(NA)

  if (stringr::str_detect(spliceInfo$offset, "\\+")) {
    return(c(
      paste0("c.", spliceInfo$position, "+1"),
      paste0("c.", spliceInfo$position, "+2")
    ))
  } else if (stringr::str_detect(spliceInfo$offset, "\\-")) {
    return(c(
      paste0("c.", spliceInfo$position, "-1"),
      paste0("c.", spliceInfo$position, "-2")
    ))
  } else {
    return(NA)
  }
}

#' nearby_motif function
#' Get the nearby U2 motif nucleotide positions
#' @param hgvsc hgvsc
#'
#' @return nearby U2 motif nucleotide positions
nearby_motif <- function(hgvsc) {
  spliceInfo <- spliceJunction(hgvsc)
  if (is.na(spliceInfo$offset))
    return(NA)

  if (stringr::str_detect(spliceInfo$offset, "\\+")) {
    positions <- paste0("c.", spliceInfo$position, "+", 3:6)
    positions <- c(positions, paste0("c.", spliceInfo$position), paste0("c.", spliceInfo$position - 1), paste0("c.", spliceInfo$position - 2))
    return(positions)
  } else if (stringr::str_detect(spliceInfo$offset, "\\-")) {
    positions <- paste0("c.", spliceInfo$position, "-", 20:3)
    positions <- c(positions, paste0("c.", spliceInfo$position))
    return(positions)
  } else {
    return(NA)
  }
}

#' canonicalSpliceSite function
#' Check whether a canonical splice site variant
#' @param hgvsc hgvsc
#'
#' @return canonical splice site variant or not
canonicalSpliceSite <- function(hgvsc) {
  if (grepl("c[.][0-9]+[+][12]{1}[A-Za-z]+", hgvsc)) {
    return(T)
  } else if (grepl("c[.][0-9]+[-][12]{1}[A-Za-z]+", hgvsc)) {
    return(T)
  } else {
    return(F)
  }
}

#' checkSynonymousSite function
#' Check whether this synonymous variant is in last 3 bases of an exon or 1 st base of an exon
#' @param coding_start coding start
#' @param exon exon
#'
#' @return T if synonymous variant is in last 3 bases of an exon or 1 st base of an exon
checkSynonymousSite <- function(coding_start, exon) {
  if (!all(is.na(exon))) {
    if (exon$CodingStart == coding_start) { # first base of an exon
      return(F)
    } else if ((exon$CodingEnd - coding_start) < 3) { # last 3 bases of an exon
      return(F)
    } else {
      return(T)
    }
  }

  return(F)
}

#' checkIntronicSite function
#' Check whether this intronic variant is beyond +7 or -21
#' @param hgvsc hgvsc
#'
#' @return intronic variant is beyond +7 or -21
checkIntronicSite <- function(hgvsc) {
  intronic_pos <- as.numeric(stringr::str_extract(hgvsc, "([+]|[-])\\d*"))
  if ((!is.na(hgvsc)) && (!is.na(intronic_pos))) {
    if ((intronic_pos >= 7) && (intronic_pos > 0)) {
      return(T)
    }
    if ((intronic_pos < 0) && (intronic_pos <= -21)) {
      return(T)
    }
  }

  return(F)
}

#' extract_max_spliceAI function
#' Extract max spliceAI score
#' @param splice dataframe with splice AI scores
#'
#' @return max spliceAI score
extract_max_spliceAI <- function(splice) {
  # Split the string, remove non-numeric characters, and convert to numeric
  numbers <- as.numeric(unlist(strsplit(gsub("[^0-9\\.]", ",", splice), ",")))
  # Return NA if no numbers are present
  if (all(is.na(numbers))) {
    return(NA)
  }
  # Otherwise, return the maximum value
  max(numbers, na.rm = TRUE)
}

#' evaluate_splice function
#' evaluate splice site or intronic variant
#' @param hgvsc hgvsc
#' @param variant_info named list with all the information required for curation
#' @param eval_PVS1 PVS1 decision
#' @param eval_PP3 PP3 decision
#'
#' @return curation
#'
#' @import dplyr
#' @import tidyr
evaluate_splice <- function(hgvsc, variant_info, eval_PVS1, eval_PP3) {
  spliceInfo <- spliceJunction(hgvsc)
  eval_nearby <- variant_info$var_nearby
  spliceAI_delta_score <- variant_info$spliceAI_delta_score

  spliceNotSupported <- spliceJunctionNotSupported(hgvsc)
  if (!is.na(spliceNotSupported) && spliceNotSupported)
  {
    curation <- "NotSupported"
    return(curation)
  }

  eval_same_nucleotide <- eval_nearby %>%
    filter(.data$is_same_nucleotide) %>%
    filter(.data$spliceAI_threshold < spliceAI_delta_score)
  eval_nearby_canonical <- eval_nearby %>%
    filter(.data$is_nearby_canonical) %>%
    filter(.data$spliceAI_threshold < spliceAI_delta_score)
  eval_nearby_motif <- eval_nearby %>%
    filter(.data$is_nearby_motif & !(hgvsc %in% eval_same_nucleotide$hgvsc)) %>%
    filter(.data$spliceAI_threshold < spliceAI_delta_score)

  # Decision according to PS1 Table (Walker et al., 2023, PMID: 37352859)
  if (spliceInfo$spliceInfo == "Canonical")
  {
    canonical_result <- handle_canonical(eval_nearby_canonical, eval_nearby_motif, eval_PVS1)
    curation <- canonical_result$PS1
  }
  else if (spliceInfo$spliceInfo == "U2 motif")
  {
    motif_result <- handle_motif(eval_same_nucleotide, eval_nearby_canonical, eval_nearby_motif, eval_PP3)
    curation <- motif_result$PS1
  }
  else if (spliceInfo$spliceInfo == "Exception")
  {
    curation <- "NotApplicable"
  }
  else if (spliceInfo$spliceInfo == "Outside")
  {
    curation <- "NotApplicable"
  }

  return(curation)
}


#' handle_motif function
#' @param eval_same_nucleotide evaluate same nucleotide variants
#' @param eval_nearby_canonical evaluate nearby canonical variants
#' @param eval_nearby_motif evaluate nearby motif
#' @param eval_PP3 PP3 decision
#'
#' @return curation
handle_motif <- function(eval_same_nucleotide, eval_nearby_canonical, eval_nearby_motif, eval_PP3) {
  if (is.null(eval_PP3) || is.na(eval_PP3) || eval_PP3 != "Supporting")
  {
    return(list(PS1 = "NotMet"))
  }
  else
  {
    if (sum(eval_same_nucleotide$Curation == "Pathogenic") > 0)
    {
      return(list(PS1 = "Strong"))
    }
    else if (sum(eval_same_nucleotide$Curation == "Likely Pathogenic") > 0)
    {
      return(list(PS1 = "Moderate"))
    }
    else if (
      any(eval_nearby_canonical$Curation == "Pathogenic") ||
      any(eval_nearby_motif$Curation == "Pathogenic"))
    {
      return(list(PS1 = "Moderate"))
    }
    else if (
      any(eval_nearby_canonical$Curation == "Likely Pathogenic") ||
      any(eval_nearby_motif$Curation == "Likely Pathogenic"))
    {
      return(list(PS1 = "Supporting"))
    }
    else
    {
      return(list(PS1 = "NotMet"))
    }
  }
}

#' handle_canonical function
#' @param eval_nearby_canonical evaluate nearby canonical variants
#' @param eval_nearby_motif evaluate nearby motif
#' @param eval_PVS1 PVS1 decision
#'
#' @return list with curation and curation_log
handle_canonical <- function(eval_nearby_canonical, eval_nearby_motif, eval_PVS1) {
  if (!is.null(eval_PVS1) && !is.na(eval_PVS1) && eval_PVS1 == "VeryStrong")
  {
    if (sum(eval_nearby_canonical$Curation == "Pathogenic") > 0)
    {
      return(list(PS1 = "Supporting"))
    }
    if (any(eval_nearby_motif$Curation %in% c("Pathogenic", "Likely Pathogenic")))
    {
      return(list(PS1 = "Supporting"))
    }
  }

  if (eval_PVS1 %in% c("Strong", "Moderate", "Supporting"))
  {
    if (sum(eval_nearby_canonical$Curation == "Pathogenic") > 0)
    {
      return(list(PS1 = "Strong"))
    }
    else if (sum(eval_nearby_motif$Curation == "Pathogenic") > 0)
    {
      return(list(PS1 = "Moderate"))
    }
    else if (sum(eval_nearby_motif$Curation == "Likely Pathogenic") > 0)
    {
      return(list(PS1 = "Supporting"))
    }

    # If no pathogenic comparison variants are found
    return(list(PS1 = "NotMet"))
  }

  return(list(PS1 = "NotMet"))
}

#' assignConsequence function
#' Assign a main consequence from the consequences returned by VEP
#' @param consequence consequence from VEP
#'
#' @return main consequence
assignConsequence <- function(consequence)
{
  if (grepl("stop[_]lost", consequence))
    return ("stop_lost")
  if (grepl("frameshift[_]variant", consequence))
    return ("frameshift_variant")
  if (grepl("missense[_]variant", consequence))
    return ("missense_variant")
  if (grepl("synonymous[_]", consequence))
    return ("synonymous_variant")
  if (grepl("coding[_]", consequence))
    return ("splice_variant")
  if (grepl("start[_]lost", consequence))
    return ("start_lost")
  if (grepl("stop[_]gained", consequence))
    return ("stop_gained")
  if (grepl("splice[_]acceptor[_]variant", consequence))
    return ("splice_acceptor_variant")
  if (grepl("splice[_]donor[_]variant", consequence))
    return ("splice_donor_variant")
  if (grepl("intron[_]", consequence))
    return ("intron_variant")
  if (grepl("inframe[_]", consequence))
    return ("inframe_variant")
  if (grepl("5_prime_UTR[_]", consequence))
    return ("5_prime_UTR_variant")
  if (grepl("3_prime_UTR[_]", consequence))
    return ("3_prime_UTR_variant")

  return (consequence)
}

#' calculate_updated_Pr_multi
#' @param n is the denominator, no. of germline variant observed
#' @param k is the number of times a single recurrent somatic hotspot is seen
#' @param l is the number of times a single non-recurrent somatic variant is seen
#' @param m is the number of times a multi-somatic hit is seen
#' @param p is the expected probabilities of observing k, l, m (out of n) in pathogenic variants
#' @param q is the expected probabilities of observing k, l, m (out of n) in cases without a germline variant as a surrogate of it being benign
#' @param pr_A posterior probability of germline variant being pathogenic/deleterious
#'
#' @return calculated posterior probability
calculate_updated_Pr_multi <- function(n, k, l, m, p=c(0.58, 0.10, 0.007), q=c(0.002, 0.002, 0.0015), pr_A=0.1) {
  log_message <- paste0("n:", n, " k: ", k, " l: ", l, " m: ", m)
  #cat(paste0("Odds of Pathogenicity parameters ", log_message, "\n"))
  # Create all valid combinations of supplied indices
  grid <- expand.grid(
    n = n,
    k = k,
    l = l,
    m = m
  )

  # Keep only valid multinomial counts
  grid <- subset(grid, k + l + m <= n)

  # Complete probability vectors (category 4 = "none of the above")
  p_all <- c(p, 1 - sum(p))
  q_vec <- if (length(q) == 1) rep(q, 3) else q
  q_all <- c(q_vec, 1 - sum(q_vec))

  # Remaining category counts
  neu <- grid$n - grid$k - grid$l - grid$m

  # Log multinomial coefficient
  log_multinom_coef <- lgamma(grid$n + 1) -
    (lgamma(grid$k + 1) +
       lgamma(grid$l + 1) +
       lgamma(grid$m + 1) +
       lgamma(neu + 1))

  # Log probabilities under A and not-A
  log_pr_B_givenA <- log_multinom_coef +
    grid$k * log(p_all[1]) +
    grid$l * log(p_all[2]) +
    grid$m * log(p_all[3]) +
    neu     * log(p_all[4])

  log_pr_B_givenNotA <- log_multinom_coef +
    grid$k * log(q_all[1]) +
    grid$l * log(q_all[2]) +
    grid$m * log(q_all[3]) +
    neu     * log(q_all[4])

  # Prior log-odds
  log_prior_odds <- log(pr_A / (1 - pr_A))

  # Posterior log-odds (log version of Bayes theorem)
  log_odds <- log_prior_odds + (log_pr_B_givenA - log_pr_B_givenNotA)

  # Stable posterior probability using logistic transform
  posterior_multinom <- stats::plogis(log_odds)

  # Odds-path ratio (how much the odds changed)
  odds_path <- exp(log_odds - log_prior_odds)

  # Output
  grid$posterior_multinom <- posterior_multinom
  grid$pr_B_givenA_multinom <- exp(log_pr_B_givenA)
  grid$pr_B_givenNotA_multinom <- exp(log_pr_B_givenNotA)
  grid$odds_path <- odds_path

  return(grid)
}
