# ============================================================
# 03: Additional analyses
# ============================================================
#
# This script conducts the following analyses:
#
# 1. Baseline descriptive statistics (Table 1), using
#    unimputed GSMS data and sampling weights.
#
# 2. Descriptive statistics for psychiatric outcomes in
#    adolescence and young adulthood (Supplementary Table 1),
#    using multiply imputed data.
#
# 3. Exploratory analyses (imputed data) examining:
#    - Whether cash-transfer receipt moderates concurrent
#      associations between maternal depressive symptoms and
#      offspring psychopathology before versus after the transfer.
#    - Whether cash-transfer receipt predicts subsequent
#      maternal depressive symptoms.
#
# The script uses the original GSMS dataset and previously
# imputed datasets generated in script 01.
# Multiple imputation is not rerun here.
# ============================================================

# Load required packages ----
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(haven)
library(lavaan)
library(survey)
library(lme4)
library(mice)
library(geepack)
library(flextable)
library(officer)

# Set up file paths ----

source("config.R")
data_dir <- file.path(gsms_dir, "restricted_data")
tables_dir <- file.path(gsms_dir, "results", "tables")

# Load original GSMS dataset ----

gsms_long <- read_sas(
  file.path(data_dir, "Jessie.sas7bdat")
)

# Load imputed datasets ----

imp_list <- readRDS(
  file.path(data_dir, "imp_list_derivedVars_20260818.rds")
)

# Table 1: Derive outcomes in unimputed data ----

## Derive internalising in unimputed dataset ##
gsms_long$int <- NA_real_
cells <- unique(gsms_long[, c("WAVE", "AGEGRUP")])

for (i in seq_len(nrow(cells))) {
  
  w <- cells$WAVE[i]
  a <- cells$AGEGRUP[i]
  
  idx <- gsms_long$WAVE == w & gsms_long$AGEGRUP == a
  
  Z <- scale(gsms_long[idx, c("j4nda", "jtotanx")]) # standardise by wave and cohort
  
  # Mean composite of depression and anxiety
  gsms_long$int[idx] <- rowMeans(Z, na.rm = TRUE)
}

gsms_long$int[is.nan(gsms_long$int)] <- NA_real_


## Derive externalising in unimputed dataset

## Create substance abuse
make_severity <- function(use, dep) {
  ifelse(use == 1 & dep == 1, 2,
         ifelse((use == 1 & dep == 0) | (use == 0 & dep == 1), 1,
                ifelse(use == 0 & dep == 0, 0, NA)))
}

# Code substance abuse in terms of severity
gsms_long$sev_nicotine <- make_severity(gsms_long$J4NUSE1, gsms_long$j4nsdp1)
gsms_long$sev_cocaine  <- make_severity(gsms_long$J4NUSE4, gsms_long$j4nsdp4)
gsms_long$sev_amphet   <- make_severity(gsms_long$J4NUSE5, gsms_long$j4nsdp5)
gsms_long$sev_inhal    <- make_severity(gsms_long$J4NUSE6, gsms_long$j4nsdp6)
gsms_long$sev_heroin   <- make_severity(gsms_long$J4NUSE7, gsms_long$j4nsdp7)
gsms_long$sev_halluc   <- make_severity(gsms_long$J4NUSE8, gsms_long$j4nsdp8)
gsms_long$sev_sed      <- make_severity(gsms_long$J4NUSE9, gsms_long$j4nsdp9)

gsms_long$sev_cannabis <- NA_real_
gsms_long$sev_cannabis[gsms_long$j5nthc == 0] <- 0
gsms_long$sev_cannabis[gsms_long$j5nthc %in% 1:2] <- 1
gsms_long$sev_cannabis[gsms_long$j5nthc >= 3] <- 2

gsms_long$sev_alcohol <- NA_real_
gsms_long$sev_alcohol[gsms_long$j5nalc == 0] <- 0
gsms_long$sev_alcohol[gsms_long$j5nalc %in% 1:2] <- 1
gsms_long$sev_alcohol[gsms_long$j5nalc >= 3] <- 2

## Specify substance model
library(lavaan)
substance_model <- '
Substance =~
  sev_alcohol +
  sev_cannabis +
  sev_nicotine +
  sev_cocaine +
  sev_amphet +
  sev_inhal + 
  sev_heroin +
  sev_halluc +
  sev_sed
'

# Fit the CFA for substance abuse

# Note: severity of alcohol and cannabis were derived differently
ordered_vars <- c("sev_alcohol", 
                  "sev_cannabis", "sev_nicotine", 
                  "sev_cocaine","sev_amphet",
                  "sev_inhal", "sev_heroin",
                  "sev_halluc","sev_sed"
)

gsms_long[ordered_vars] <-
  lapply(gsms_long[ordered_vars],
         function(x) ordered(x, levels = c(0,1,2)))

substance_fit <- cfa(
  substance_model,
  data = gsms_long,
  ordered = ordered_vars,
  estimator = "WLSMV"
)

# Extract substance use factor scores
fs <- lavPredict(substance_fit)[,1]
idx <- lavInspect(substance_fit, "case.idx")

gsms_long$substance_score <- NA_real_
gsms_long$substance_score[idx] <- fs

# Create externalising
ext_vars <- c(
  "j4ncd",
  "j4nopp",
  "c4nASP",
  "p4nadh",
  "substance_score"
)

gsms_long$ext_comp <- NA_real_

cells <- unique(gsms_long[, c("WAVE", "AGEGRUP")])

for (i in seq_len(nrow(cells))) {
  
  w <- cells$WAVE[i]
  a <- cells$AGEGRUP[i]
  
  idx <- gsms_long$WAVE == w &
    gsms_long$AGEGRUP == a
  
  X <- gsms_long[idx, ext_vars, drop = FALSE]
  
  keep <- sapply(X, function(x) {
    sdx <- sd(x, na.rm = TRUE)
    !is.na(sdx) && sdx > 0
  })
  
  if (sum(keep) < 2) next
  
  Z <- scale(X[, keep, drop = FALSE])
  
  n_observed <- rowSums(!is.na(Z))
  
  gsms_long$ext_comp[idx] <-
    ifelse(n_observed >= 2,
           rowMeans(Z, na.rm = TRUE),
           NA_real_)
}

gsms_long$ext_comp[is.nan(gsms_long$ext_comp)] <- NA_real_

# Table 1: Identify analytic sample ----
# Identify the 1,348 analytic participants
# At least one mental-health assessment from wave 5 onwards

analytic_ids <- gsms_long |>
  group_by(GSMSID) |>
  summarise(
    any_posttransfer_mh = any(
      WAVE %in% 5:18 & !is.na(j4nda)
    ),
    .groups = "drop"
  ) |>
  filter(any_posttransfer_mh) |>
  pull(GSMSID)

# Check
length(analytic_ids)

# Table 1: Prepare baseline dataset ----

baseline_unimp <- gsms_long |>
  filter(GSMSID %in% analytic_ids) |>
  arrange(GSMSID, WAVE) |>
  group_by(GSMSID) |>
  
  # First available pre-transfer assessment in waves 1-3
  mutate(
    baseline_wave = min(WAVE[WAVE %in% 1:3], na.rm = TRUE)
  ) |>
  
  summarise(
    
    # maternal depression: waves 1-4
    exposure_early = mean(DRISK22A[WAVE %in% 1:4], na.rm = TRUE),
    
    # fixed AI / cash-transfer indicator
    RACE_early = as.integer(any(RACEAI == 1, na.rm = TRUE)),
    
    # fixed number of AI parents
    AIpar = AIpar[1],
    
    # Covariates from first available pre-transfer assessment
    singleP_baseline = DRISK4[WAVE == baseline_wave][1],
    lowSES_baseline  = DRISK1[WAVE == baseline_wave][1],
    income_baseline  = INCOME[WAVE == baseline_wave][1],
    bmi_baseline     = bmi[WAVE == baseline_wave][1],
    int_baseline      = int[WAVE == baseline_wave][1],
    ext_baseline    = ext_comp[WAVE == baseline_wave][1],
    
    # Fixed covariates
    sex     = SEX[1],
    low_bw  = DRISK38[1],
    AGEGRUP = AGEGRUP[1],
    
    # baseline weight: keep as wave 1 unless you know GSMS averaged weights
    WT2_person = first(WT2[!is.na(WT2)]),
    
    .groups = "drop"
  ) |>
  mutate(
    exposure_early = ifelse(
      is.nan(exposure_early),
      NA_real_,
      exposure_early
    ),
    
    group = factor(
      RACE_early,
      levels = c(0, 1),
      labels = c("Non American-Indian / no transfer",
                 "American Indian / cash-transfer")
    ),
    
    AGEGRUP = factor(
      AGEGRUP,
      levels = c(1, 2, 3),
      labels = c("Youngest cohort", "Middle cohort", "Oldest cohort")
    ),
    
    AIpar = factor(
      AIpar,
      levels = c(2, 1, 0),
      labels = c("2", "1", "0")
    ),
    
    # 0 = girl, 1 = boy
    sex = factor(
      sex,
      levels = c(0, 1),
      labels = c("Female", "Male")
    ),
    
    singleP_baseline = factor(
      singleP_baseline,
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),
    
    lowSES_baseline = factor(
      lowSES_baseline,
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),
    
    low_bw = factor(
      low_bw,
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
  )

# Keep only rows with non-missing, positive baseline weights
baseline_unimp_wt <- baseline_unimp |>
  filter(!is.na(WT2_person), WT2_person > 0)

# Table 1: Calculate weighted descriptives ----

design <- svydesign(
  ids = ~1,
  weights = ~WT2_person,
  data = baseline_unimp_wt
)

# Helper functions
fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < .001) return("<0.001")
  sprintf("%.3f", p)
}

w_sd <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  x <- x[ok]
  w <- w[ok]
  if (length(x) <= 1) return(NA_real_)
  m <- weighted.mean(x, w)
  sqrt(sum(w * (x - m)^2) / sum(w))
}

p_cat <- function(var) {
  keep <- !is.na(design$variables[[var]]) & !is.na(design$variables[["group"]])
  dsub <- subset(design, keep)
  out <- tryCatch(
    svychisq(as.formula(paste0("~", var, " + group")), dsub),
    error = function(e) NULL
  )
  if (is.null(out)) return(NA_real_)
  out$p.value
}

p_cont <- function(var) {
  keep <- !is.na(design$variables[[var]]) & !is.na(design$variables[["group"]])
  dsub <- subset(design, keep)
  fit <- tryCatch(
    svyglm(as.formula(paste0(var, " ~ group")), design = dsub),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)
  coefs <- summary(fit)$coefficients
  if (nrow(coefs) < 2) return(NA_real_)
  coefs[2, 4]
}

# Multi-category variables
# unweighted n, weighted %, and p-values
multi_vars <- c("AGEGRUP", "AIpar", "sex")

multi <- map_dfr(multi_vars, function(v) {
  
  # Counts and weighted percentages by cash-transfer group
  dat <- baseline_unimp_wt |>
    filter(!is.na(.data[[v]]), !is.na(group)) |>
    group_by(group, level = .data[[v]]) |>
    summarise(
      n = n(),
      weighted_n = sum(WT2_person),
      .groups = "drop"
    ) |>
    group_by(group) |>
    mutate(pct = 100 * weighted_n / sum(weighted_n)) |>
    ungroup()
  
  # Counts and weighted percentages for the full sample
  total <- baseline_unimp_wt |>
    filter(!is.na(.data[[v]])) |>
    group_by(level = .data[[v]]) |>
    summarise(
      n = n(),
      weighted_n = sum(WT2_person),
      .groups = "drop"
    ) |>
    mutate(pct = 100 * weighted_n / sum(weighted_n))
  
  dat |>
    tidyr::complete(
      group,
      level,
      fill = list(n = 0, weighted_n = 0, pct = 0)
    ) |>
    mutate(
      variable = v,
      stat = sprintf("%d (%.1f)", n, pct)
    ) |>
    select(variable, level, group, stat) |>
    pivot_wider(names_from = group, values_from = stat) |>
    left_join(
      total |>
        mutate(Total = sprintf("%d (%.1f)", n, pct)) |>
        select(level, Total),
      by = "level"
    ) |>
    mutate(
      section = case_when(
        variable == "AGEGRUP" ~ "Age cohort",
        variable == "AIpar" ~ "Number of American Indian parents",
        variable == "sex" ~ "Sex"
      ),
      row_order = case_when(
        variable == "AGEGRUP" & level == "Youngest cohort" ~ 1,
        variable == "AGEGRUP" & level == "Middle cohort" ~ 2,
        variable == "AGEGRUP" & level == "Oldest cohort" ~ 3,
        variable == "AIpar" & level == "2" ~ 4,
        variable == "AIpar" & level == "1" ~ 5,
        variable == "AIpar" & level == "0" ~ 6,
        variable == "sex" & level == "Female" ~ 7,
        variable == "sex" & level == "Male" ~ 8,
        TRUE ~ 999
      )
    ) |>
    group_by(variable) |>
    mutate(
      p_value = if_else(row_number() == 1, fmt_p(p_cat(v)), "")
    ) |>
    ungroup()
})

# Binary variables
binary <- tribble(
  ~var, ~positive_level, ~label, ~section,
  "lowSES_baseline", "Yes", "Low SES", "Childhood risk factors (baseline)",
  "singleP_baseline", "Yes", "Family instability", "Childhood risk factors (baseline)",
  "low_bw", "Yes", "Low birth weight", "Childhood risk factors (baseline)"
)

bin_table <- pmap_dfr(binary, function(var, positive_level, label, section) {
  dat <- baseline_unimp_wt |>
    filter(!is.na(.data[[var]]), !is.na(group)) |>
    group_by(group) |>
    summarise(
      n = sum(.data[[var]] == positive_level, na.rm = TRUE),
      pct = 100 * weighted.mean(.data[[var]] == positive_level, WT2_person, na.rm = TRUE),
      .groups = "drop"
    )
  
  total <- baseline_unimp_wt |>
    filter(!is.na(.data[[var]])) |>
    summarise(
      n = sum(.data[[var]] == positive_level, na.rm = TRUE),
      pct = 100 * weighted.mean(.data[[var]] == positive_level, WT2_person, na.rm = TRUE)
    )
  
  dat |>
    mutate(stat = sprintf("%d (%.1f)", n, pct)) |>
    select(group, stat) |>
    pivot_wider(names_from = group, values_from = stat) |>
    mutate(
      label = label,
      section = section,
      Total = sprintf("%d (%.1f)", total$n, total$pct),
      p_value = fmt_p(p_cat(var))
    )
})

# Continuous variables

cont <- tribble(
  ~var,              ~label,                      ~section,
  "income_baseline", "Family income",             "Childhood risk factors (baseline)",
  "bmi_baseline",    "BMI",                       "Childhood risk factors (baseline)",
  "int_baseline",    "Internalising symptoms",    "Childhood risk factors (baseline)",
  "ext_baseline",    "Externalising symptoms",    "Childhood risk factors (baseline)",
  "exposure_early",  "Early depressive symptoms", "Maternal depressive symptoms (waves 1–4)"
)

cont_table <- pmap_dfr(cont, function(var, label, section) {
  dat <- baseline_unimp_wt |>
    group_by(group) |>
    summarise(
      m = weighted.mean(.data[[var]], WT2_person, na.rm = TRUE),
      sd = w_sd(.data[[var]], WT2_person),
      .groups = "drop"
    )
  
  total <- baseline_unimp_wt |>
    summarise(
      m = weighted.mean(.data[[var]], WT2_person, na.rm = TRUE),
      sd = w_sd(.data[[var]], WT2_person)
    )
  
  dat |>
    mutate(stat = sprintf("%.2f (%.2f)", m, sd)) |>
    select(group, stat) |>
    pivot_wider(names_from = group, values_from = stat) |>
    mutate(
      label = label,
      section = section,
      Total = sprintf("%.2f (%.2f)", total$m, total$sd),
      p_value = fmt_p(p_cont(var))
    )
})

# Table 1: Assemble and format table ----

table1 <- bind_rows(
  multi |>
    transmute(
      section = section,
      label = as.character(level),
      `Non American-Indian / no transfer`,
      `American Indian / cash-transfer`,
      Total,
      p_value
    ),
  bin_table |>
    select(section, label,
           `Non American-Indian / no transfer`,
           `American Indian / cash-transfer`,
           Total, p_value),
  cont_table |>
    select(section, label,
           `Non American-Indian / no transfer`,
           `American Indian / cash-transfer`,
           Total, p_value)
)

# Add section labels only on first row of each section
table1 <- table1 |>
  mutate(
    row_order = case_when(
      section == "Age cohort" & label == "Youngest cohort" ~ 1,
      section == "Age cohort" & label == "Middle cohort" ~ 2,
      section == "Age cohort" & label == "Oldest cohort" ~ 3,
      section == "Number of American Indian parents" & label == "2" ~ 4,
      section == "Number of American Indian parents" & label == "1" ~ 5,
      section == "Number of American Indian parents" & label == "0" ~ 6,
      section == "Sex" & label == "Female" ~ 7,
      section == "Sex" & label == "Male" ~ 8,
      section == "Childhood risk factors (baseline)" & label == "Low SES" ~ 9,
      section == "Childhood risk factors (baseline)" & label == "Family instability" ~ 10,
      section == "Childhood risk factors (baseline)" & label == "Family income" ~ 11,
      section == "Childhood risk factors (baseline)" & label == "Low birth weight" ~ 12,
      section == "Childhood risk factors (baseline)" & label == "BMI" ~ 13,
      section == "Childhood risk factors (baseline)" & label == "Internalising symptoms" ~ 14,
      section == "Childhood risk factors (baseline)" & label == "Externalising symptoms" ~ 15,
      section == "Maternal depression symptoms" ~ 16,
      TRUE ~ 999
    )
  ) |>
  arrange(row_order) |>
  group_by(section) |>
  mutate(section_print = if_else(row_number() == 1, section, "")) |>
  ungroup() |>
  select(
    section_print,
    label,
    `Non American-Indian / no transfer`,
    `American Indian / cash-transfer`,
    Total,
    p_value
  )

# Blank out repeated p-values within multi-row sections
table1 <- table1 |>
  mutate(
    p_value = case_when(
      label %in% c("Middle cohort", "Oldest cohort", "1", "0") ~ "",
      TRUE ~ p_value
    )
  )

# Flextable formatting
ft <- flextable(table1) |>
  set_header_labels(
    section_print = "",
    label = "",
    `Non American-Indian / no transfer` = "Non American-Indian / no transfer",
    `American Indian / cash-transfer` = "American Indian / cash-transfer",
    Total = "Total",
    p_value = "p"
  ) |>
  merge_v(j = "section_print") |>
  valign(j = "section_print", valign = "top") |>
  bold(j = "section_print", bold = TRUE) |>
  align(j = 1:2, align = "left", part = "all") |>
  align(j = 3:6, align = "center", part = "all") |>
  theme_booktabs() |>
  autofit()

ft <- add_footer_lines(
  ft,
  values = "Values are unweighted n (weighted %) for categorical variables and weighted mean (SD) for continuous variables. Participants with missing baseline weights were excluded from weighted analyses. P values are from weighted tests using WT2."
)

# Table 1: Export to Word ----
doc <- read_docx() |>
  body_add_flextable(ft)

print(doc, target = file.path(tables_dir, "Table1_weighted.docx"))

# Supplementary Table: Later psychiatric outcomes ----
# Calculate weighted descriptives within each imputed dataset,
# then average estimates across imputations

# Construct datasets needed for table
make_analysis_dat_period <- function(dat) {
  
  # Derive person-level early/baseline variables
  person_vars <- dat |>
    arrange(GSMSID, WAVE) |>
    group_by(GSMSID) |>
    
    # First available pre-transfer assessment
    mutate(
      baseline_wave = min(WAVE[WAVE %in% 1:3], na.rm = TRUE)
    ) |>
    
    summarise(
      
      # Exposure = mean maternal depression score across waves 1-4
      exposure_early = mean(
        mat_dep[WAVE %in% 1:4],
        na.rm = TRUE
      ),
      
      # Fixed race indicator
      RACE_early = as.integer(
        any(RACEAI == 1, na.rm = TRUE)
      ),
      
      # Baseline covariates:
      # first available assessment in waves 1-3
      singleP_baseline = DRISK4[WAVE == baseline_wave][1],
      income_baseline  = INCOME[WAVE == baseline_wave][1],
      lowSES_baseline  = DRISK1[WAVE == baseline_wave][1],
      bmi_baseline     = bmi[WAVE == baseline_wave][1],
      int_baseline     = int[WAVE == baseline_wave][1],
      ext_baseline     = ext_comp[WAVE == baseline_wave][1],
      
      # Fixed covariates
      sex     = SEX[1],
      low_bw  = DRISK38[1],
      AGEGRUP = AGEGRUP[1],
      
      .groups = "drop"
    ) |>
    mutate(
      # Replace NaN if all maternal depression values were missing
      exposure_early = ifelse(
        is.nan(exposure_early),
        NA_real_,
        exposure_early
      )
    )
  
  # Post-transfer outcomes
  # ----------------------------------------------------------
  analysis_dat <- dat |>
    filter(WAVE %in% 5:18) |>
    mutate(
      period = case_when(
        COHAGE <= 21 ~ "adolescence",
        COHAGE >= 22 & COHAGE <= 30 ~ "young_adulthood",
        TRUE ~ NA_character_
      ),
      period = factor(
        period,
        levels = c(
          "adolescence",
          "young_adulthood"
        )
      )
    ) |>
    select(
      GSMSID,
      WAVE,
      period,
      COHAGE,
      WT2,
      int,
      ext_comp,
      any_suicidality_broad
    ) |>
    left_join(
      person_vars,
      by = "GSMSID"
    )
  
  # Standardise maternal depression within analytic sample
  # using one observation per participant
  exposure_scaling <- analysis_dat |>
    distinct(
      GSMSID,
      exposure_early
    ) |>
    summarise(
      mean_exp = mean(
        exposure_early,
        na.rm = TRUE
      ),
      sd_exp = sd(
        exposure_early,
        na.rm = TRUE
      )
    )
  
  
  analysis_dat <- analysis_dat |>
    mutate(
      exposure_early_z =
        (exposure_early - exposure_scaling$mean_exp) /
        exposure_scaling$sd_exp
    )
  
  return(analysis_dat)
}


# Create developmental-period analytic datasets for all imputations
analysis_list <- lapply(imp_list, make_analysis_dat_period)

#  Weighted mean and SD function
w_mean_sd <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  x <- x[ok]
  w <- w[ok]
  
  m <- weighted.mean(x, w)
  v <- sum(w * (x - m)^2) / sum(w)
  s <- sqrt(v)
  
  c(mean = m, sd = s)
}

# Summarise one imputed dataset
summarise_dataset <- function(dat) {
  
  dat %>%
    group_by(period) %>%
    summarise(
      obs_n = n(),
      ind_n = n_distinct(GSMSID),
      
      int_mean = w_mean_sd(int, WT2)[["mean"]],
      int_sd   = w_mean_sd(int, WT2)[["sd"]],
      
      ext_mean = w_mean_sd(ext_comp, WT2)[["mean"]],
      ext_sd   = w_mean_sd(ext_comp, WT2)[["sd"]],
      
      suic_pct = 100 * weighted.mean(any_suicidality_broad == 1, WT2, na.rm = TRUE),
      
      .groups = "drop"
    )
}

# Apply across imputations + pool
desc_pooled <- lapply(analysis_list, summarise_dataset) %>%
  bind_rows(.id = "imp") %>%
  group_by(period) %>%
  summarise(
    obs_n = round(mean(obs_n)),
    ind_n = round(mean(ind_n)),
    
    int_mean = mean(int_mean, na.rm = TRUE),
    int_sd   = mean(int_sd, na.rm = TRUE),
    
    ext_mean = mean(ext_mean, na.rm = TRUE),
    ext_sd   = mean(ext_sd, na.rm = TRUE),
    
    suic_pct = mean(suic_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    period = recode(period,
                    "adolescence" = "Adolescence",
                    "young_adulthood" = "Young adulthood"
    )
  )

# Format table
desc_table <- desc_pooled %>%
  mutate(
    `Observations, n` = as.character(obs_n),
    `Individuals, n`  = as.character(ind_n),
    
    `Internalising composite score, mean (SD)` =
      sprintf("%.2f (%.2f)", int_mean, int_sd),
    
    `Externalising composite score, mean (SD)` =
      sprintf("%.2f (%.2f)", ext_mean, ext_sd),
    
    `Suicidality, %` =
      sprintf("%.2f", suic_pct)
  ) %>%
  select(
    period,
    `Observations, n`,
    `Individuals, n`,
    `Internalising composite score, mean (SD)`,
    `Externalising composite score, mean (SD)`,
    `Suicidality, %`
  ) %>%
  pivot_longer(
    cols = -period,
    names_to = "Characteristic",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = period,
    values_from = Value
  )

# Flextable (no title, clean formatting)
ft_desc <- flextable(desc_table) %>%
  autofit() %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  align(align = "left", part = "all") 


ft_desc

# Exploratory analyses: Concurrent associations ----

# Test if maternal depression was associated with child concurrent outcomes to 
# different extent in cash transfer group before and after cash transfer 

# Create the analysis dataset

make_analysis_dat_did <- function(dat) {
  
  # ----------------------------------------------------------
  # Person-level variables
  # ----------------------------------------------------------
  
  person_vars <- dat |>
    arrange(GSMSID, WAVE) |>
    group_by(GSMSID) |>
    
    # First available pre-transfer assessment
    mutate(
      baseline_wave = min(WAVE[WAVE %in% 1:3], na.rm = TRUE)
    ) |>
    
    summarise(
      
      # Cash-transfer / AI group
      RACE_early = as.integer(
        any(RACEAI == 1, na.rm = TRUE)
      ),
      
      # Baseline covariates
      singleP_baseline = DRISK4[WAVE == baseline_wave][1],
      income_baseline  = INCOME[WAVE == baseline_wave][1],
      lowSES_baseline  = DRISK1[WAVE == baseline_wave][1],
      bmi_baseline     = bmi[WAVE == baseline_wave][1],
      int_baseline = int[WAVE == baseline_wave][1],
      ext_baseline = ext_comp[WAVE == baseline_wave][1],
      
      # Fixed covariates
      sex     = SEX[1],
      low_bw  = DRISK38[1],
      AGEGRUP = AGEGRUP[1],
      WT2     = WT2[1],
      
      .groups = "drop"
    ) 
  
  # ----------------------------------------------------------
  # Repeated outcome data: pre and post cash transfer
  # ----------------------------------------------------------
  
  analysis_dat <- dat |>
    filter(WAVE %in% 1:8) |>
    mutate(
      
      # 0 = before cash transfer
      # 1 = after cash transfer
      post_transfer = ifelse(WAVE >= 5, 1, 0)
      
    ) |>
    select(
      GSMSID,
      WAVE,
      COHAGE,
      mat_dep,
      int,
      ext_comp,
      any_suicidality_broad,
      post_transfer
    ) |>
    left_join(person_vars, by = "GSMSID") |>
    mutate(
      RACE_early = factor(
        RACE_early,
        levels = c(0, 1)
      ),
      
      post_transfer = factor(
        post_transfer,
        levels = c(0, 1),
        labels = c("pre", "post")
      ),
      
      AGEGRUP = factor(
        AGEGRUP,
        levels = c(1, 2, 3)
      )
    )
  
  
  # ----------------------------------------------------------
  # Standardise maternal depression and psychiatric outcomes
  # across observations in waves 1–8, within each imputed dataset
  # ----------------------------------------------------------
  
  analysis_dat <- analysis_dat |>
    mutate(
      mat_dep_z = as.numeric(
        scale(mat_dep)
      ),
      int_z = as.numeric(
        scale(int)
      ),
      ext_comp_z = as.numeric(
        scale(ext_comp)
      )
    )
  
  return(analysis_dat)
}


analysis_list_did <- lapply(
  imp_list,
  make_analysis_dat_did
)

## Check age range
analysis_list_did[[1]] |>
  group_by(post_transfer) |>
  summarise(
    min_age = min(COHAGE, na.rm = TRUE),
    max_age = max(COHAGE, na.rm = TRUE)
  )


# Fit concurrent association models across imputed datasets
fit_concurrent <- function(outcome) {
  
  formula <- as.formula(
    paste0(
      outcome, " ~ ",
      "mat_dep_z * RACE_early * post_transfer + ",
      "COHAGE + sex + AGEGRUP + ",
      "singleP_baseline + income_baseline + lowSES_baseline + ",
      "low_bw + bmi_baseline + int_baseline + ext_baseline + ",
      "(1 | GSMSID)"
    )
  )
  
  models <- lapply(
    analysis_list_did,
    function(dat) {
      lmer(
        formula,
        data = dat,
        weights = WT2,
        REML = FALSE
      )
    }
  )
  
  pool(as.mira(models))
}

# Internalising symptoms
pool_concur_int <- fit_concurrent("int_z")
summary(pool_concur_int, conf.int = TRUE)

# Externalising symptoms
pool_concur_ext <- fit_concurrent("ext_comp_z")
summary(pool_concur_ext, conf.int = TRUE)

## Suicidality

## Run weighted GEE separately in each imputed dataset
concur_suic <- lapply(
  analysis_list_did,
  function(dat) {
    
    # geeglm requires complete data for all variables in the model
    dat_gee <- dat %>%
      filter(
        !is.na(any_suicidality_broad),
        !is.na(mat_dep_z),
        !is.na(RACE_early),
        !is.na(post_transfer),
        !is.na(AGEGRUP),
        !is.na(singleP_baseline),
        !is.na(income_baseline),
        !is.na(lowSES_baseline),
        !is.na(low_bw),
        !is.na(bmi_baseline),
        !is.na(int_baseline),
        !is.na(ext_baseline),
        !is.na(WT2),
        !is.na(GSMSID)
      ) %>%
      mutate(
        GSMSID_gee = as.integer(factor(GSMSID))
      ) %>%
      arrange(GSMSID_gee)
    
    geeglm(
      any_suicidality_broad ~
        mat_dep_z * RACE_early * post_transfer +
        AGEGRUP +
        singleP_baseline +
        income_baseline +
        lowSES_baseline +
        low_bw +
        bmi_baseline +
        int_baseline +
        ext_baseline,
      id = GSMSID_gee,
      data = dat_gee,
      weights = WT2,
      family = binomial(link = "logit"),
      corstr = "exchangeable"
    )
  }
)

pool_concur_suic <- pool(as.mira(concur_suic))
summary(pool_concur_suic, conf.int = TRUE)

### Extract the 3-way interaction from each model

term_concur <- "mat_dep_z:RACE_early1:post_transferpost"

get_pooled_term <- function(pool_object, target_term) {
  summary(pool_object, conf.int = TRUE) |>
    filter(.data$term == target_term) |>
    select(estimate, `2.5 %`, `97.5 %`, p.value)
}

int_concur_result  <- get_pooled_term(pool_concur_int, term_concur)
ext_concur_result  <- get_pooled_term(pool_concur_ext, term_concur)
suic_concur_result <- get_pooled_term(pool_concur_suic, term_concur)

# Check that exactly one interaction term was extracted per model
stopifnot(
  nrow(int_concur_result) == 1,
  nrow(ext_concur_result) == 1,
  nrow(suic_concur_result) == 1
)

## Format functions
format_est_ci <- function(est, low, high) {
  sprintf("%.3f (%.3f, %.3f)", est, low, high)
}

format_p <- function(p) {
  if (length(p) == 0 || is.na(p)) return(NA_character_)
  if (p < .001) return("<.001")
  sprintf("%.3f", p)
}

## Create the three GSMS concurrent-association rows
panel_gsms_concurrent <- bind_rows(
  tibble(
    Analysis = "Moderation of concurrent associations",
    Outcome = "Internalising symptoms",
    `Estimate (95% CI)` = format_est_ci(
      int_concur_result$estimate,
      int_concur_result$`2.5 %`,
      int_concur_result$`97.5 %`
    ),
    `P value` = format_p(int_concur_result$p.value)
  ),
  tibble(
    Analysis = "",
    Outcome = "Externalising symptoms",
    `Estimate (95% CI)` = format_est_ci(
      ext_concur_result$estimate,
      ext_concur_result$`2.5 %`,
      ext_concur_result$`97.5 %`
    ),
    `P value` = format_p(ext_concur_result$p.value)
  ),
  tibble(
    Analysis = "",
    Outcome = "Suicidality",
    `Estimate (95% CI)` = format_est_ci(
      suic_concur_result$estimate,
      suic_concur_result$`2.5 %`,
      suic_concur_result$`97.5 %`
    ),
    `P value` = format_p(suic_concur_result$p.value)
  )
)

panel_gsms_concurrent


# Exploratory analyses: Subsequent maternal depression ----
## Derive dataset with maternal depression post-cash-transfer
make_analysis_dat_matdep <- function(dat) {
  
  person_vars <- dat |>
    arrange(GSMSID, WAVE) |>
    group_by(GSMSID) |>
    # First available pre-transfer assessment
    mutate(
      baseline_wave = min(WAVE[WAVE %in% 1:3], na.rm = TRUE)
    ) |>
    
    summarise(
      
      # Cash-transfer / AI group
      RACE_early = as.integer(
        any(RACEAI == 1, na.rm = TRUE)
      ),
      
      # Baseline covariates
      singleP_baseline = DRISK4[WAVE == baseline_wave][1],
      income_baseline  = INCOME[WAVE == baseline_wave][1],
      lowSES_baseline  = DRISK1[WAVE == baseline_wave][1],
      bmi_baseline     = bmi[WAVE == baseline_wave][1],
      int_baseline = int[WAVE == baseline_wave][1],
      ext_baseline = ext_comp[WAVE == baseline_wave][1],
      
      # Fixed covariates
      sex     = SEX[1],
      low_bw  = DRISK38[1],
      AGEGRUP = AGEGRUP[1],
      
      .groups = "drop"
    )
  
  # Select maternal depression after the cash transfer
  analysis_dat <- dat |>
    filter(WAVE %in% 5:18) |>
    select(
      GSMSID,
      WAVE,
      mat_dep,
      WT2
    ) |>
    left_join(person_vars, by = "GSMSID")
  
  return(analysis_dat)
}

analysis_list_matdep <- lapply(imp_list, make_analysis_dat_matdep)
library(lme4)

## Test if cash-transfer predicts later maternal depression symptoms 
models_matdep <- lapply(
  analysis_list_matdep,
  function(dat) {
    lmer(
      scale(mat_dep) ~ RACE_early + WAVE +
        AGEGRUP +
        singleP_baseline +
        income_baseline +
        lowSES_baseline +
        low_bw +
        bmi_baseline +
        int_baseline +
        ext_baseline +
        (1 | GSMSID),
      data = dat,
      weights = WT2
    )
  }
)

# Extract cash-transfer group coefficient for subsequent maternal depression
pool_matdep <- pool(as.mira(models_matdep))

matdep_result <- summary(pool_matdep, conf.int = TRUE) %>%
  filter(term == "RACE_early")

matdep_result

# Exploratory analyses: Combine in panel and save results ----

# Check that the maternal depression coefficient was extracted
stopifnot(nrow(matdep_result) == 1)

# Format maternal depression result
panel_matdep <- tibble(
  Analysis = "Subsequent maternal depressive symptoms",
  Outcome = "Maternal depressive symptoms",
  `Estimate (95% CI)` = format_est_ci(
    matdep_result$estimate,
    matdep_result$`2.5 %`,
    matdep_result$`97.5 %`
  ),
  `P value` = format_p(matdep_result$p.value)
)

# Combine concurrent associations and maternal depression results
panel_gsms <- bind_rows(
  panel_gsms_concurrent,
  panel_matdep
)

panel_gsms

# Save results
saveRDS(
  panel_gsms,
  file.path(tables_dir, "panel_gsms_exploratory.rds")
)
