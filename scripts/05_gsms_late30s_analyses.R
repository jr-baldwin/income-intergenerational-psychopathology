# ============================================================
# 05: Supplementary analyses of late-30s psychiatric outcomes
# ============================================================
#
# This script conducts supplementary analyses examining:
#
# 1. Associations between early maternal depressive symptoms
#    and offspring psychiatric outcomes in their late 30s.
#
# 2. Moderation of these associations by cash-transfer receipt.
#
# Outcomes include depressive symptoms, alcohol-related
# problems, cannabis-related problems, and other drug-related
# problems.
#
# These analyses were conducted to examine outcomes measured
# at the late-30s assessment, which were not directly
# comparable with psychiatric measures at earlier waves.
#
# The script uses previously imputed GSMS datasets generated
# in script 01 and a separate intergenerational-wave dataset.
# ============================================================

# Load required packages ----
library(dplyr)
library(tidyr)
library(purrr)
library(haven)
library(mice)
library(flextable)

# Set up file paths ----
source("config.R")

data_dir <- file.path(gsms_dir, "restricted_data")
tables_dir <- file.path(gsms_dir, "results", "tables")

# Load imputed GSMS datasets ----
imp_list <- readRDS(
  file.path(data_dir, "imp_list_derivedVars_20260818.rds")
)

# Load late-30s dataset ----
# has a number of variables from the most recent intergenerational wave. One row per subject.
ig_wave <- read_sas(
  file.path(data_dir, "intergen.sas7bdat")
)

# Prepare late-30s outcomes for merging with imputed datasets ----
# Rename GSMS ID var for consistency with other dataset for merging
ig_wave <- ig_wave %>%
  dplyr::rename(GSMSID = gsmsid) %>%
  dplyr::rename(RACEAI = raceai)    

stopifnot(!anyDuplicated(ig_wave$GSMSID))

# Standardise age 38 outcomes
ig_wave <- ig_wave |>
  mutate(
    BCK19tot_z = as.numeric(scale(BCK19tot)),
    aud19tot_z = as.numeric(scale(aud19tot)),
    cud19tot_z = as.numeric(scale(cud19tot)),
    das19tot_z = as.numeric(scale(das19tot))
  )

## Derive person-level predictors and merge datasets ----
imp_list_38_unique <- lapply(imp_list, function(dat) {
  
  person_level <- dat |>
    arrange(GSMSID, WAVE) |>
    group_by(GSMSID) |>
    mutate(
      baseline_wave = min(WAVE[WAVE %in% 1:3], na.rm = TRUE)
    ) |>
    summarise(
      
      # Early maternal depressive symptoms
      exposure_early = if (
        all(is.na(mat_dep[WAVE %in% 1:4]))
      ) {
        NA_real_
      } else {
        mean(mat_dep[WAVE %in% 1:4], na.rm = TRUE)
      },
      
      # Cash-transfer indicator
      RACE_early = as.integer(any(RACEAI == 1, na.rm = TRUE)),
      
      # Covariates from first pre-transfer assessment
      singleP_baseline = DRISK4[WAVE == baseline_wave][1],
      income_baseline  = INCOME[WAVE == baseline_wave][1],
      lowSES_baseline  = DRISK1[WAVE == baseline_wave][1],
      bmi_baseline     = bmi[WAVE == baseline_wave][1],
      mh_baseline      = j4_core[WAVE == baseline_wave][1],
      
      # Fixed covariates
      sex     = SEX[1],
      low_bw  = DRISK38[1],
      AGEGRUP = AGEGRUP[1],
      WT2     = WT2[1],
      
      .groups = "drop"
    ) |>
    mutate(
      exposure_early_z = as.numeric(scale(exposure_early))
    )
  
  person_level |>
    left_join(
      ig_wave |>
        select(
          GSMSID,
          BCK19tot_z,
          aud19tot_z,
          cud19tot_z,
          das19tot_z
        ),
      by = "GSMSID"
    )
})

# Fit association and moderation models ----
fit_age38 <- function(outcome, moderation = FALSE) {
  
  exposure <- if (moderation) {
    "exposure_early_z * RACE_early"
  } else {
    "exposure_early_z"
  }
  
  formula <- as.formula(paste(
    outcome, "~", exposure,
    "+ sex + AGEGRUP + singleP_baseline",
    "+ income_baseline + lowSES_baseline",
    "+ low_bw + bmi_baseline + mh_baseline"
  ))
  
  lapply(imp_list_38_unique, function(dat) {
    lm(formula, data = dat, weights = WT2)
  })
}

# Associations with early maternal depressive symptoms
fit_dep_rq1a <- fit_age38("BCK19tot_z")
fit_aud_rq1a <- fit_age38("aud19tot_z")
fit_cud_rq1a <- fit_age38("cud19tot_z")
fit_das_rq1a <- fit_age38("das19tot_z")

# Moderation by cash transfer
fit_dep_rq3 <- fit_age38("BCK19tot_z", moderation = TRUE)
fit_aud_rq3 <- fit_age38("aud19tot_z", moderation = TRUE)
fit_cud_rq3 <- fit_age38("cud19tot_z", moderation = TRUE)
fit_das_rq3 <- fit_age38("das19tot_z", moderation = TRUE)

# Specify models for supplementary table ----

# Combine models and specify the coefficient to extract
model_info <- tribble(
  ~Research_question,        ~Outcome,                       ~Model,       ~Term,
  "Association",             "Depressive symptoms",          fit_dep_rq1a, "exposure_early_z",
  "Association",             "Alcohol-related problems",     fit_aud_rq1a, "exposure_early_z",
  "Association",             "Cannabis-related problems",    fit_cud_rq1a, "exposure_early_z",
  "Association",             "Other drug-related problems",  fit_das_rq1a, "exposure_early_z",
  "Cash-transfer moderation", "Depressive symptoms",         fit_dep_rq3,  "exposure_early_z:RACE_early",
  "Cash-transfer moderation", "Alcohol-related problems",    fit_aud_rq3,  "exposure_early_z:RACE_early",
  "Cash-transfer moderation", "Cannabis-related problems",   fit_cud_rq3,  "exposure_early_z:RACE_early",
  "Cash-transfer moderation", "Other drug-related problems", fit_das_rq3,  "exposure_early_z:RACE_early"
)

# Pool estimates across imputations ----

# Function to pool one model and extract the requested coefficient
extract_pooled_result <- function(model_list, term_name) {
  
  pooled_summary <- model_list |>
    mice::as.mira() |>
    mice::pool() |>
    summary(conf.int = TRUE)
  
  result <- pooled_summary |>
    filter(.data$term == term_name)
  
  stopifnot(nrow(result) == 1)
  
  result |>
    transmute(
      estimate = estimate,
      lower_ci = `2.5 %`,
      upper_ci = `97.5 %`,
      p_value = p.value
    )
}

# Extract and combine results
age38_results <- pmap_dfr(
  model_info,
  function(Research_question, Outcome, Model, Term) {
    
    result <- extract_pooled_result(
      model_list = Model,
      term_name = Term
    )
    
    bind_cols(
      tibble(
        Research_question = Research_question,
        Outcome = Outcome
      ),
      result
    )
  }
)

# Format supplementary results table ----
age38_table <- age38_results |>
  mutate(
    `β (95% CI)` = sprintf(
      "%.3f (%.3f, %.3f)",
      estimate,
      lower_ci,
      upper_ci
    ),
    `P value` = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  ) |>
  select(
    `Research question` = Research_question,
    Outcome,
    `β (95% CI)`,
    `P value`
  )

age38_table

# Create formatted Word table ----

ft_age38 <- flextable(age38_table) |>
  
  # Merge repeated research question labels
  merge_v(j = "Research question") |>
  
  # Header labels
  set_header_labels(
    `Research question` = "Research question",
    Outcome = "Late 30s outcome",
    `β (95% CI)` = "β (95% CI)",
    `P value` = "P value"
  ) |>
  
  # Formatting
  theme_booktabs() |>
  bold(part = "header") |>
  font(fontname = "Arial", part = "all") |>
  fontsize(size = 10, part = "all") |>
  
  align(
    j = c("Research question", "Outcome"),
    align = "left",
    part = "all"
  ) |>
  
  align(
    j = c("β (95% CI)", "P value"),
    align = "center",
    part = "all"
  ) |>
  
  valign(valign = "center", part = "all") |>
  valign(j = "Research question", valign = "top") |>
  
  width(j = "Research question", width = 2.1) |>
  width(j = "Outcome", width = 2.2) |>
  width(j = "β (95% CI)", width = 1.8) |>
  width(j = "P value", width = 0.8) |>
  
  autofit()

ft_age38

# Export supplementary table to Word ----
save_as_docx(
  `Table` = ft_age38,
  path = file.path(tables_dir, "age38_results.docx")
)
