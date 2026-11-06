# ============================================================
# 04: Three-generation analyses
# ============================================================
#
# This script conducts exploratory analyses examining
# associations between psychopathology across three generations:
#
# 1. Association between G1 maternal depressive symptoms
#    and G3 offspring psychopathology.
#
# 2. Association between G1 maternal depressive symptoms
#    and G3 psychopathology, adjusting for G2 psychopathology.
#
# 3. Moderation of the G1-G3 association by cash-transfer
#    receipt, including a Bayes factor for the interaction.
#
# The script uses previously imputed GSMS datasets generated
# in script 01 and a separate G3 dataset.
# ============================================================

# Load required packages ----
library(dplyr)
library(tidyr)
library(haven)
library(lme4)
library(mice)
library(psych)

# Set up file paths ----
source("config.R")
data_dir <- file.path(gsms_dir, "restricted_data")

# ============================================================
# PART 1: Prepare G1-G3 dataset (run only when needed)
# ============================================================

# Load imputed GSMS datasets ----
imp_list <- readRDS(
  file.path(data_dir, "imp_list_derivedVars_20260818.rds")
)

# Load G3 dataset ----
g3 <- read_sas(
  file.path(data_dir, "intgen_gsmsyr19_g3long31aug2021.sas7bdat")
)

# Prepare G1/G2 person-level variables ----
g2_person_list <- lapply(imp_list, function(dat) {
  
  dat %>%
    arrange(GSMSID, WAVE) %>%
    group_by(GSMSID) %>%
    
    # Keep each person's first available pre-transfer assessment
    mutate(
      baseline_wave = min(WAVE[WAVE %in% 1:3], na.rm = TRUE)
    ) |>
    
    summarise(
      
      ## G1/G2 predictors used in the G3 analyses
      
      # Mean maternal depression across waves 1–4
      mat_dep_early = if (
        all(is.na(mat_dep[WAVE %in% 1:4]))
      ) {
        NA_real_
      } else {
        mean(mat_dep[WAVE %in% 1:4], na.rm = TRUE)
      },
      
      # Mean post-transfer internalising symptoms across waves 5–18
      int_post_mean = if (
        all(is.na(int[WAVE %in% 5:18]))
      ) {
        NA_real_
      } else {
        mean(int[WAVE %in% 5:18], na.rm = TRUE)
      },
      
      # Mean post-transfer externalising symptoms across waves 5–18
      ext_post_mean = if (
        all(is.na(ext_comp[WAVE %in% 5:18]))
      ) {
        NA_real_
      } else {
        mean(ext_comp[WAVE %in% 5:18], na.rm = TRUE)
      },
      
      # Any post-transfer broad suicidality across waves 5–18
      suicidality_post_ever = if (
        all(is.na(any_suicidality_broad[WAVE %in% 5:18]))
      ) {
        NA_integer_
      } else {
        as.integer(
          any(
            any_suicidality_broad[WAVE %in% 5:18] == 1,
            na.rm = TRUE
          )
        )
      },
      
      ## Covariates, defined consistently with the other analyses
      RACE_early = as.integer(any(RACEAI == 1, na.rm = TRUE)),
      singleP_baseline = DRISK4[WAVE == baseline_wave][1],
      income_baseline  = INCOME[WAVE == baseline_wave][1],
      lowSES_baseline  = DRISK1[WAVE == baseline_wave][1],
      bmi_baseline     = bmi[WAVE == baseline_wave][1],
      mh_baseline      = j4_core[WAVE == baseline_wave][1],
      
      sex     = SEX[1],
      low_bw  = DRISK38[1],
      AGEGRUP = AGEGRUP[1],
      WT2     = WT2[1],
      
      .groups = "drop"
    )
})

# Check mean age
describe(g3$c19dp25a)
# Check whether SDQ was provided
table(g3$KID19esd)

## Before merging, make a clean G3 analysis file with only:
# the key linking variables
# G3 outcome variables
# G3 covariates that could be relevant

## What to keep
# GSMSID — to link children to G2 parent
# a child identifier 
# G3 psychopathology outcome variables
# G3 child covariates such as age and sex
# any assessment variables, such as interview wave/date if relevant

# Check number of participants who are not biological children of G2 participants by relationship type
with(g3, table(child, c19dp25re, useNA = "always"))

g3_small <- g3 %>%
  dplyr::select(gsmsid, 
                child,  # B1=biological child 1, B2=biological child 2, etc... NB=non-biological child
                c19dp25re, # child relationship to parent: c19dp25re: childs relationship to tc, 0=bio, 1=adopt, 2=foster, 3=step, 4=partner's kid, 5=sibling, 6=niece/nephew, 7=other
                c19dp25a,  # child age
                c19dp25g, # child gender
                sdq19emo, # SDQ Emotional Problem Scale (Sum across items 3, 8, 13, 16, 24) .D=deceased .O=age>17 .S=Not elig .Y=Age<2
                sdq19con, # SDQ Conduct Problem Scale (Sum across items 5, 7r, 12, 18, 22) 
                sdq19hyp, # SDQ Hyperactivity Scale (Sum across items 2, 10 15, 21r, 25r)
                sdq19per, # SDQ Peer Problems Scale (Sum across items 6, 11r, 14r, 19, 23)
                sdq19pro, # SDQ Prosocial Scale (Sum across items 1, 4, 9, 17, 20)
                KID19esd, # Participant Eligible to Complete SDQ for this Child based on contact level .D=deceased .O=age>17 .R=Missing Contact Items .Y=Age<2
                KID19gsd # Participant recieved SDQ items for this child
  ) %>%
  filter(c19dp25re == 0) # Select only biological children of GSMS participants

## Rename gsmsid to GSMSID ahead of merging
names(g3_small)[names(g3_small) == "gsmsid"] <- "GSMSID"

## Check structure
nrow(g3_small) # 1809 rows 

table(g3_small$c19dp25re, useNA = "always")

# Check number of biological children per GSMS parent (ranges from 1-12)
g3_small %>%
  count(GSMSID) %>%
  count(n, name = "n_parents")

# Number of parents whose children completed SDQ
g3_small %>%
  filter(c19dp25re == 0) %>%
  summarise(
    n_parents_all = n_distinct(GSMSID),
    n_parents_sdq = n_distinct(GSMSID[KID19gsd == 1])
  )


# Check one row is returned per child
g3_small %>%
  count(GSMSID, child) %>%
  filter(n > 1)

table(g3_small$KID19esd, useNA = "always")
table(g3_small$KID19gsd, useNA = "always")

# Merge G1/G2 and G3 datasets ----

# This will create one row per G3 child and repeated G1/G2 predictors for siblings sharing GSMSID
analysis_list_g1g2g3 <- lapply(g2_person_list, function(dat) {
  g3_small %>%
    left_join(
      dat,
      by = "GSMSID",
      relationship = "many-to-one"
    )
})

# Checks after merging
# Number of rows (1809, the same)
nrow(g3_small)
nrow(analysis_list_g1g2g3[[1]])

# Check how many children matched a G2 record
analysis_list_g1g2g3[[1]] %>%
  summarise(
    matched = sum(!is.na(RACE_early)),
    total = n()
  )

# Make G3 SDQ total difficulties score
analysis_list_g1g2g3 <- lapply(analysis_list_g1g2g3, function(dat) {
  dat %>%
    mutate(
      sdq_n = rowSums(!is.na(cbind(sdq19emo, sdq19con, sdq19hyp, sdq19per))),
      sdq_total = ifelse(
        sdq_n == 4,
        rowSums(cbind(sdq19emo, sdq19con, sdq19hyp, sdq19per), na.rm = TRUE),
        NA)
    )
})

# Check G3 analytic sample ----

# Check number of children and parents
analysis_list_g1g2g3[[1]] |>
  summarise(
    n_children = n(),
    n_parents = n_distinct(GSMSID),
    n_valid_sdq = sum(!is.na(sdq_total))
  )

# Confirm only biological children represented
analysis_list_g1g2g3[[1]] %>%
  count(c19dp25re)

# Check final analytic sample across imputations
sample_counts <- lapply(analysis_list_g1g2g3, function(dat) {
  dat |>
    filter(
      !is.na(sdq_total),
      !is.na(mat_dep_early),
      !is.na(int_post_mean),
      !is.na(ext_post_mean),
      !is.na(suicidality_post_ever)
    ) |>
    summarise(
      n_children = n(),
      n_parents = n_distinct(GSMSID)
    )
})

# Save merged G1-G3 dataset ----
#saveRDS(analysis_list_g1g2g3, file.path(data_dir, "imp_g1g2g3_derived_20260819.rds"))

# ============================================================
# PART 2: Three-generation analyses
# ============================================================

# Load previously prepared G1-G3 dataset ----
analysis_list_g1g2g3 <- readRDS(
  file.path(data_dir, "imp_g1g2g3_derived_20260819.rds")
)

# Reload merged G1-G3 dataset ----
analysis_list_g1g2g3 <- readRDS(
  file.path(data_dir, "imp_g1g2g3_derived_20260819.rds")
)

## Main vars for test:
# G1 Maternal depression: early (prior to cash transfer)
# G2 psychopathology: post-cash transfer (int, ext, then suicidality)
# G3 psychopathology: SDQ total difficulties

# Standardise variables and define analytic sample ----
analysis_list_g1g2g3 <- lapply(analysis_list_g1g2g3, function(dat) {
  
  dat %>%
    filter(
      !is.na(sdq_total),
      !is.na(mat_dep_early),
      !is.na(int_post_mean),
      !is.na(ext_post_mean),
      !is.na(suicidality_post_ever)
    ) %>%
    mutate(
      sdq_total_z      = as.numeric(scale(sdq_total)),
      mat_dep_early_z  = as.numeric(scale(mat_dep_early)),
      int_post_mean_z  = as.numeric(scale(int_post_mean)),
      ext_post_mean_z  = as.numeric(scale(ext_post_mean))
    )
  
})

# Check mean age of children
analysis_list_g1g2g3[[1]] |>
  summarise(
    n_children = n(),
    mean_age = mean(c19dp25a, na.rm = TRUE),
    sd_age = sd(c19dp25a, na.rm = TRUE),
    min_age = min(c19dp25a, na.rm = TRUE),
    max_age = max(c19dp25a, na.rm = TRUE),
    missing_age = sum(is.na(c19dp25a))
  )

# G1-G3 association: Unadjusted for G2 psychopathology ----

# Maternal depression and G3 SDQ total difficulties
fit_list <- lapply(analysis_list_g1g2g3, function(dat) {
  lmer(sdq_total_z ~ mat_dep_early_z +
         (1 | GSMSID),
       data = dat,
       weights = WT2,
       REML = FALSE
  )
})

pool_g1_g3 <- pool(as.mira(fit_list))
summary(pool_g1_g3, conf.int = TRUE)
# check sample size 
sapply(fit_list, nobs)

# G1-G3 association: Adjusted for G2 psychopathology ----
fit_list <- lapply(analysis_list_g1g2g3, function(dat) {
  lmer(sdq_total_z ~ mat_dep_early_z + int_post_mean_z + ext_post_mean_z
       + suicidality_post_ever + 
         (1 | GSMSID),
       data = dat,
       weights = WT2,
       REML = FALSE
  )
})
pool_g1_g3_g2 <- pool(as.mira(fit_list))
summary(pool_g1_g3_g2, conf.int = TRUE)

# G1-G3 association: Cash-transfer moderation ----
fit_list <- lapply(analysis_list_g1g2g3, function(dat) {
  lmer(sdq_total_z ~ mat_dep_early_z*RACE_early +
         (1 | GSMSID),
       data = dat,
       weights = WT2,
       REML = FALSE
  )
})
pooled <- pool(as.mira(fit_list))
summary(pooled, conf.int = TRUE)

# Cash-transfer moderation: Bayes factor ----
bf_sdq_interaction <- lapply(analysis_list_g1g2g3, function(dat) {
  
  m0 <- lmer(
    sdq_total_z ~ mat_dep_early_z + RACE_early + (1 | GSMSID),
    data = dat,
    weights = WT2,
    REML = FALSE
  )
  
  m1 <- lmer(
    sdq_total_z ~ mat_dep_early_z * RACE_early + (1 | GSMSID),
    data = dat,
    weights = WT2,
    REML = FALSE
  )
  
  bic0 <- BIC(m0)
  bic1 <- BIC(m1)
  bf01 <- exp((bic1 - bic0) / 2)
  
  data.frame(
    bic_null = bic0,
    bic_full = bic1,
    bf01 = bf01
  )
})

bf_sdq_interaction_df <- do.call(rbind, bf_sdq_interaction)

bf_sdq_interaction_df
median(bf_sdq_interaction_df$bf01)
range(bf_sdq_interaction_df$bf01)
