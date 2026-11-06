# ============================================================
# 02: Main analyses
# ============================================================
#
# This script conducts the following analyses:
#
# 1. Associations between maternal depressive symptoms
#    and offspring psychopathology.
#
# 2. Differences in these associations between adolescence
#    and young adulthood.
#
# 3. Moderation of these associations by cash-transfer
#    receipt, including Bayes factors for interaction effects.
#
# 4. Heterogeneity analyses examining cash-transfer amount
#    and duration of exposure.
#
# The script uses previously imputed datasets generated
# in script 01. Multiple imputation is not rerun here.
# ============================================================

# Load required packages ----

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(lme4)
library(mice)
library(broom.mixed)
library(flextable)
library(officer)
library(stringr)

# Set up file paths ----

source("config.R")
data_dir <- file.path(gsms_dir, "restricted_data")
tables_dir <- file.path(gsms_dir, "results", "tables")

# Load imputed datasets ----

imp_list <- readRDS(
  file.path(data_dir, "imp_list_derivedVars_20260818.rds")
)

# Prepare analytic dataset for RQ1 ----

make_analysis_dat <- function(dat) {
  
  person_vars <- dat |>
    arrange(GSMSID, WAVE) |>
    group_by(GSMSID) |>
    
    # Keep each person's first available pre-transfer assessment
    mutate(
      baseline_wave = min(WAVE[WAVE %in% 1:3], na.rm = TRUE)
    ) |>
    
    summarise(
      
      # Mean maternal depression across waves 1–4
      exposure_early = if (
        all(is.na(mat_dep[WAVE %in% 1:4]))
      ) {
        NA_real_
      } else {
        mean(mat_dep[WAVE %in% 1:4], na.rm = TRUE)
      },
      
      # Fixed race indicator
      RACE_early = as.integer(any(RACEAI == 1, na.rm = TRUE)),
      
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
    ) |>
    mutate(
      exposure_early_z = as.numeric(scale(exposure_early))
    )
  
  # Outcomes: psychiatric outcomes assessed at waves 5–18
  analysis_dat <- dat |>
    filter(WAVE %in% 5:18) |>
    select(
      GSMSID,
      WAVE,
      int,
      ext_comp,
      any_suicidality,
      any_suicidality_broad,
      WT2
    ) |>
    left_join(person_vars, by = "GSMSID")
  
  return(analysis_dat)
}

## Apply to all imputed datasets
analysis_list <- lapply(imp_list, make_analysis_dat)

# Standardise outcome vars
analysis_list <- lapply(analysis_list, function(dat) {
  dat %>%
    mutate(
      int_z = as.numeric(scale(int)),
      ext_comp_z = as.numeric(scale(ext_comp))
    )
})

# Check analytic sample and missingness ----

## Check missingness in internalising problems for every wave post-cash-transfer
missing_by_wave <- map_dfr(
  seq_along(analysis_list),
  function(i) {
    analysis_list[[i]] %>%
      group_by(WAVE) %>%
      summarise(
        n_total = n(),
        n_missing = sum(is.na(int)),
        pct_missing = mean(is.na(int)) * 100,
        .groups = "drop"
      ) %>%
      mutate(imputation = i, .before = 1)
  }
)

print(missing_by_wave, n=320) # ranges from 1154 (wave 18) to 308 (wave 14)
table(analysis_list[[1]]$WAVE)
length(unique(analysis_list[[1]]$GSMSID)) #1348 participants are represented

# Check the number of individuals post-cash transfer who had ANY internalising assessment (1348)
sapply(analysis_list, function(dat) {
  dat %>%
    filter(!is.na(int)) %>%
    summarise(n = n_distinct(GSMSID)) %>%
    pull(n)
})

## Check missingness in externalising problems for every wave post-cash-transfer
missing_by_wave <- map_dfr(
  seq_along(analysis_list),
  function(i) {
    analysis_list[[i]] %>%
      group_by(WAVE) %>%
      summarise(
        n_total = n(),
        n_missing = sum(is.na(ext_comp)),
        pct_missing = mean(is.na(ext_comp)) * 100,
        .groups = "drop"
      ) %>%
      mutate(imputation = i, .before = 1)
  }
)

print(missing_by_wave, n=320) # ranges from 1154 (wave 18) to 308 (wave 14)
table(analysis_list[[1]]$WAVE)
length(unique(analysis_list[[1]]$GSMSID)) #1348 participants are represented

# Check the number of individuals post-cash transfer who had ANY externalising assessment (1348)
sapply(analysis_list, function(dat) {
  dat %>%
    filter(!is.na(ext_comp)) %>%
    summarise(n = n_distinct(GSMSID)) %>%
    pull(n)
})

# Check the number of individuals post-cash transfer who had ANY suicidality assessment (1338)
sapply(analysis_list, function(dat) {
  dat %>%
    filter(!is.na(any_suicidality_broad)) %>%
    summarise(n = n_distinct(GSMSID)) %>%
    pull(n)
})

# Check number of individuals in imputed dataset who had info on cash transfer (1420)
sapply(imp_list, function(dat) {
  dat %>%
    filter(!is.na(RACEAI)) %>%
    summarise(n = n_distinct(GSMSID)) %>%
    pull(n)
})

# RQ1: Maternal depression and offspring psychopathology ----

## Fit the model to internalising problems
fit_int <- lapply(analysis_list, function(dat) {
  lmer(
    int_z ~ exposure_early_z + 
      WAVE + 
      sex +
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
    weights = WT2,   
    REML = FALSE
  )
})

pool_int <- pool(as.mira(fit_int))
summary(pool_int, conf.int = TRUE)

sapply(fit_int, function(m) length(unique(model.frame(m)$GSMSID)))

### Externalising
fit_ext <- lapply(analysis_list, function(dat) {
  lmer(
    ext_comp_z ~ exposure_early_z + 
      WAVE + 
      sex +
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
    weights = WT2,  
    REML = FALSE
  )
})

pool_ext <- pool(as.mira(fit_ext))
summary(pool_ext, conf.int = TRUE)
sapply(fit_ext, function(m) length(unique(model.frame(m)$GSMSID)))

## Suicidality
fit_sui <- lapply(analysis_list, function(dat) {
  glmer(
    any_suicidality_broad ~ exposure_early_z + 
      WAVE + 
      sex +
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
    weights = WT2,   # <<< ADDED
    family = binomial(link = "logit"),
    nAGQ = 0,
    control = glmerControl(optimizer = "bobyqa")
  )
})
pool_sui <- pool(as.mira(fit_sui))
summary(pool_sui, conf.int = TRUE, exponentiate = TRUE)
sapply(fit_sui, function(m) length(unique(model.frame(m)$GSMSID)))

# RQ1b: Differences between developmental periods ----

make_analysis_dat <- function(dat) {
  
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
analysis_list <- lapply(imp_list, make_analysis_dat)

# Check developmental-period classification
table(analysis_list[[1]]$period, useNA = "ifany")

# RQ1b: Fit developmental-period models ----

# Internalising + externalising 

analysis_list_adolref <- lapply(analysis_list, function(dat) {
  dat %>% mutate(period = factor(period, levels = c("adolescence", "young_adulthood")))
})

analysis_list_yaref <- lapply(analysis_list, function(dat) {
  dat %>% mutate(period = factor(period, levels = c("young_adulthood", "adolescence")))
})

# Internalising
fit_int_period_adolref <- lapply(analysis_list_adolref, function(dat) {
  lmer(
    scale(int) ~ exposure_early_z * period + WAVE + sex + AGEGRUP +
      singleP_baseline + income_baseline + lowSES_baseline +
      low_bw + bmi_baseline + int_baseline +
      ext_baseline + (1 | GSMSID),
    data = dat, weights = WT2, REML = FALSE
  )
})

fit_int_period_yaref <- lapply(analysis_list_yaref, function(dat) {
  lmer(
    scale(int) ~ exposure_early_z * period + WAVE + sex + AGEGRUP +
      singleP_baseline + income_baseline + lowSES_baseline +
      low_bw + bmi_baseline + int_baseline +
      ext_baseline + (1 | GSMSID),
    data = dat, weights = WT2, REML = FALSE
  )
})


# Externalising
fit_ext_period_adolref <- lapply(analysis_list_adolref, function(dat) {
  lmer(
    scale(ext_comp) ~ exposure_early_z * period + WAVE + sex + AGEGRUP +
      singleP_baseline + income_baseline + lowSES_baseline +
      low_bw + bmi_baseline + int_baseline +
      ext_baseline + (1 | GSMSID),
    data = dat, weights = WT2, REML = FALSE
  )
})

fit_ext_period_yaref <- lapply(analysis_list_yaref, function(dat) {
  lmer(
    scale(ext_comp) ~ exposure_early_z * period + WAVE + sex + AGEGRUP +
      singleP_baseline + income_baseline + lowSES_baseline +
      low_bw + bmi_baseline + int_baseline +
      ext_baseline + (1 | GSMSID),
    data = dat, weights = WT2, REML = FALSE
  )
})

#  Suicidality (developmental change) 

fit_sui_period_adolref <- lapply(analysis_list_adolref, function(dat) {
  glmer(
    any_suicidality_broad ~ exposure_early_z * period +
      WAVE + sex + AGEGRUP + singleP_baseline + income_baseline +
      lowSES_baseline + low_bw + bmi_baseline + int_baseline +
      ext_baseline +
      (1 | GSMSID),
    data = dat,
    weights = WT2,
    family = binomial(link = "logit"),
    nAGQ = 0,
    control = glmerControl(optimizer = "bobyqa")
  )
})

fit_sui_period_yaref <- lapply(analysis_list_yaref, function(dat) {
  glmer(
    any_suicidality_broad ~ exposure_early_z * period +
      WAVE + sex + AGEGRUP + singleP_baseline + income_baseline +
      lowSES_baseline + low_bw + bmi_baseline + int_baseline +
      ext_baseline +
      (1 | GSMSID),
    data = dat,
    weights = WT2,
    family = binomial(link = "logit"),
    nAGQ = 0,
    control = glmerControl(optimizer = "bobyqa")
  )
})

# RQ1b: Extract results and create Table 2 ----

# HELPER FUNCTIONS

extract_exposure_row <- function(model, outcome, period, exponentiate = FALSE) {
  summary(model, conf.int = TRUE, exponentiate = exponentiate) %>%
    as.data.frame() %>%
    filter(term == "exposure_early_z") %>%
    mutate(
      outcome = outcome,
      period = period
    ) %>%
    transmute(
      outcome,
      period,
      estimate,
      conf.low = `2.5 %`,
      conf.high = `97.5 %`,
      p.value
    )
}

extract_interaction_row <- function(model, outcome, exponentiate = FALSE) {
  out <- summary(model, conf.int = TRUE, exponentiate = exponentiate) %>%
    as.data.frame()
  
  interaction_row <- out %>%
    filter(
      str_detect(term, "^exposure_early_z:period") |
        str_detect(term, "^period.*:exposure_early_z$")
    )
  
  if (nrow(interaction_row) != 1) {
    stop(paste("Could not uniquely identify interaction term for", outcome))
  }
  
  interaction_row %>%
    mutate(outcome = outcome) %>%
    transmute(
      outcome,
      estimate,
      conf.low = `2.5 %`,
      conf.high = `97.5 %`,
      p.value
    )
}

fmt_est_ci <- function(est, low, high, digits = 3) {
  paste0(
    sprintf(paste0("%.", digits, "f"), est),
    " (",
    sprintf(paste0("%.", digits, "f"), low),
    ", ",
    sprintf(paste0("%.", digits, "f"), high),
    ")"
  )
}

fmt_p <- function(p) {
  ifelse(p < .001, "<.001", sprintf("%.3f", p))
}

# Internalising

rq1b_int_adol <- extract_exposure_row(
  pool(as.mira(fit_int_period_adolref)),
  outcome = "Internalising problems",
  period = "Adolescence"
)

rq1b_int_ya <- extract_exposure_row(
  pool(as.mira(fit_int_period_yaref)),
  outcome = "Internalising problems",
  period = "Young adulthood"
)

rq1b_int_interaction <- extract_interaction_row(
  pool(as.mira(fit_int_period_adolref)),
  outcome = "Internalising problems"
)

# Externalising

rq1b_ext_adol <- extract_exposure_row(
  pool(as.mira(fit_ext_period_adolref)),
  outcome = "Externalising problems",
  period = "Adolescence"
)

rq1b_ext_ya <- extract_exposure_row(
  pool(as.mira(fit_ext_period_yaref)),
  outcome = "Externalising problems",
  period = "Young adulthood"
)

rq1b_ext_interaction <- extract_interaction_row(
  pool(as.mira(fit_ext_period_adolref)),
  outcome = "Externalising problems"
)

# Suicidality

rq1b_sui_adol <- extract_exposure_row(
  pool(as.mira(fit_sui_period_adolref)),
  outcome = "Suicidality",
  period = "Adolescence",
  exponentiate = FALSE
)

rq1b_sui_ya <- extract_exposure_row(
  pool(as.mira(fit_sui_period_yaref)),
  outcome = "Suicidality",
  period = "Young adulthood",
  exponentiate = FALSE
)

rq1b_sui_interaction <- extract_interaction_row(
  pool(as.mira(fit_sui_period_adolref)),
  outcome = "Suicidality",
  exponentiate = FALSE
)

# Simple slopes

rq1b_slopes <- bind_rows(
  rq1b_int_adol, rq1b_int_ya,
  rq1b_ext_adol, rq1b_ext_ya,
  rq1b_sui_adol, rq1b_sui_ya
) %>%
  mutate(
    outcome = factor(outcome, levels = c("Internalising problems", "Externalising problems", "Suicidality")),
    period = recode(period,
                    "Adolescence" = "adol",
                    "Young adulthood" = "ya"),
    effect_ci = fmt_est_ci(estimate, conf.low, conf.high),
    p.value = fmt_p(p.value)
  ) %>%
  select(outcome, period, effect_ci, p.value) %>%
  pivot_wider(
    names_from = period,
    values_from = c(effect_ci, p.value)
  )

# Interaction terms

rq1b_interactions <- bind_rows(
  rq1b_int_interaction,
  rq1b_ext_interaction,
  rq1b_sui_interaction
) %>%
  mutate(
    outcome = factor(outcome, levels = c("Internalising problems", "Externalising problems", "Suicidality")),
    interaction_effect_ci = fmt_est_ci(estimate, conf.low, conf.high),
    interaction_p = fmt_p(p.value)
  ) %>%
  select(outcome, interaction_effect_ci, interaction_p)

# Final table
table_rq1b <- rq1b_slopes %>%
  left_join(rq1b_interactions, by = "outcome") %>%
  rename(
    Outcome = outcome,
    `Adolescence coefficient (95% CI)` = effect_ci_adol,
    `Adolescence p` = p.value_adol,
    `Adulthood coefficient (95% CI)` = effect_ci_ya,
    `Adulthood p` = p.value_ya,
    `Interaction coefficient (95% CI)` = interaction_effect_ci,
    `Interaction p` = interaction_p
  ) %>%
  select(
    Outcome,
    `Adolescence coefficient (95% CI)`,
    `Adolescence p`,
    `Adulthood coefficient (95% CI)`,
    `Adulthood p`,
    `Interaction coefficient (95% CI)`,
    `Interaction p`
  )

ft <- table_rq1b %>%
  flextable() %>%
  autofit() %>%
  theme_booktabs() %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  bold(part = "header") %>%
  fontsize(size = 10, part = "all")

ft

doc <- read_docx() |>
  body_add_flextable(ft)
print(doc, target = file.path(tables_dir, "Table2.docx"))

# RQ2: Moderation by cash-transfer receipt ----

make_analysis_dat_rq2 <- function(dat) {
  
  # Derive person-level early/baseline variables
  person_vars <- dat |>
    arrange(GSMSID, WAVE) |>
    group_by(GSMSID) |>
    
    # First available pre-transfer assessment
    mutate(
      baseline_wave = min(WAVE[WAVE %in% 1:3], na.rm = TRUE)
    ) |>
    
    summarise(
      
      # exposure = mean maternal depression score across waves 1-4
      exposure_early = mean(mat_dep[WAVE %in% 1:4], na.rm = TRUE),
      
      # fixed race indicator
      RACE_early = as.integer(any(RACEAI == 1, na.rm = TRUE)),
      
      # baseline covariates
      singleP_baseline = DRISK4[WAVE == baseline_wave][1],
      income_baseline  = INCOME[WAVE == baseline_wave][1],
      lowSES_baseline  = DRISK1[WAVE == baseline_wave][1],
      bmi_baseline     = bmi[WAVE == baseline_wave][1],
      int_baseline = int[WAVE == baseline_wave][1],
      ext_baseline = ext_comp[WAVE == baseline_wave][1],
      
      # fixed covariates
      sex     = SEX[1],
      low_bw  = DRISK38[1],
      AGEGRUP = AGEGRUP[1],
      WT2     = WT2[1],
      AIpar   = AIpar[1],
      
      .groups = "drop"
    ) |>
    mutate(
      exposure_early = ifelse(is.nan(exposure_early), NA, exposure_early),
      exposure_early_z = as.numeric(scale(exposure_early))
    )
  
  analysis_dat <- dat |>
    filter(WAVE %in% 5:18) |> 
    select(GSMSID, WAVE, COHAGE, int, ext_comp, any_suicidality, any_suicidality_broad) |>
    left_join(person_vars, by = "GSMSID") |>
    mutate(
      RACE_early = factor(RACE_early, levels = c(0, 1)),
      AGEGRUP = factor(AGEGRUP, levels = c(1, 2, 3)),
      AIpar = factor(AIpar, levels = c(0, 1, 2))
    )
  
  return(analysis_dat)
}

analysis_list_rq2 <- lapply(imp_list, make_analysis_dat_rq2)

## Create restricted datasets within each imputed dataset and standardise outcomes
analysis_split <- lapply(analysis_list_rq2, function(dat) {
  list(
    adol = dat |>
      filter(COHAGE <= 21) |>
      mutate(
        age_c = COHAGE - 16,
        int_z = as.numeric(scale(int)),
        ext_comp_z = as.numeric(scale(ext_comp))
      ),
    
    young_adult = dat |>
      filter(COHAGE >= 25, COHAGE <= 30) |>
      mutate(
        age_c = COHAGE - 26,
        int_z = as.numeric(scale(int)),
        ext_comp_z = as.numeric(scale(ext_comp))
      )
  )
})

adol_list <- lapply(analysis_split, `[[`, "adol")
young_adult_list <- lapply(analysis_split, `[[`, "young_adult")

# Relevel RACE_early so that 1 is the reference category
adol_list_race1 <- lapply(adol_list, function(dat) {
  dat |>
    mutate(RACE_early = relevel(RACE_early, ref = "1"))
})

young_adult_list_race1 <- lapply(young_adult_list, function(dat) {
  dat |>
    mutate(RACE_early = relevel(RACE_early, ref = "1"))
})

# Overall analytic sample in each imputation
sapply(analysis_list_rq2, function(dat) {
  dplyr::n_distinct(dat$GSMSID)
})

# Adolescence sample before fitting the model
sapply(adol_list, function(dat) {
  dplyr::n_distinct(dat$GSMSID)
})

# RQ2: Adolescence models ----

## Internalising 
fit_int_adol <- lapply(adol_list, function(dat) {
  lmer(
    int_z ~ exposure_early_z * RACE_early + WAVE + 
      sex +
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
    weights=WT2,
    REML = FALSE
  )
})

pool_int_adol <- pool(as.mira(fit_int_adol))
summary(pool_int_adol, conf.int = TRUE)
sapply(fit_int_adol, function(m) length(unique(model.frame(m)$GSMSID)))

## Externalising 
fit_ext_adol <- lapply(adol_list, function(dat) {
  lmer(
    ext_comp_z ~ exposure_early_z * RACE_early + WAVE + 
      sex +
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
    weights=WT2,
    REML = FALSE
  )
})

pool_ext_adol <- pool(as.mira(fit_ext_adol))
summary(pool_ext_adol, conf.int = TRUE)

## Suicidality 
fit_sui_adol <- lapply(adol_list, function(dat) {
  glmer(
    any_suicidality_broad ~ exposure_early_z * RACE_early + WAVE + 
      sex +
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
    weights=WT2,
    family = binomial(link = "logit"),
    nAGQ = 0,
    control = glmerControl(optimizer = "bobyqa")
  )
})

pool_sui_adol <- pool(as.mira(fit_sui_adol))
summary(pool_sui_adol, conf.int = TRUE, exponentiate = TRUE)

# RQ2: Young adulthood models ----
# Internalising 
fit_int_ya <- lapply(young_adult_list, function(dat) {
  lmer(
    int_z ~ exposure_early_z * RACE_early + WAVE + 
      sex +
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
    weights=WT2,
    REML = FALSE
  )
})

pool_int_ya <- pool(as.mira(fit_int_ya))
summary(pool_int_ya, conf.int = TRUE)

## Externalising 
fit_ext_ya <- lapply(young_adult_list, function(dat) {
  lmer(
    ext_comp_z ~ exposure_early_z * RACE_early + WAVE + 
      sex +
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
    weights = WT2,
    REML = FALSE
  )
})

pool_ext_ya <- pool(as.mira(fit_ext_ya))
summary(pool_ext_ya, conf.int = TRUE)

## Suicidality 
fit_sui_ya <- lapply(young_adult_list, function(dat) {
  glmer(
    any_suicidality_broad ~ exposure_early_z * RACE_early + WAVE + 
      sex +
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
    weights = WT2,
    family = binomial(link = "logit"),
    nAGQ = 0,
    control = glmerControl(optimizer = "bobyqa")
  )
})

pool_sui_ya <- pool(as.mira(fit_sui_ya))
summary(pool_sui_ya, conf.int = TRUE, exponentiate = TRUE)

# Check Ns in the model
sapply(fit_int_adol, function(m) length(unique(m@frame$GSMSID)))
sapply(fit_int_ya, function(m) length(unique(m@frame$GSMSID)))

# RQ2: Relevelled models for group-specific estimates ----

# Internalising: adolescence
fit_int_adol_race1 <- lapply(adol_list_race1, function(dat) {
  lmer(
    int_z ~ exposure_early_z * RACE_early + WAVE +
      sex +
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
    weights = WT2,
    REML = FALSE
  )
})

pool_int_adol_race1 <- pool(as.mira(fit_int_adol_race1))

# Internalising: young adulthood
fit_int_ya_race1 <- lapply(young_adult_list_race1, function(dat) {
  lmer(
    int_z ~ exposure_early_z * RACE_early + WAVE +
      sex +
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
    weights = WT2,
    REML = FALSE
  )
})

pool_int_ya_race1 <- pool(as.mira(fit_int_ya_race1))

# Externalising: adolescence
fit_ext_adol_race1 <- lapply(adol_list_race1, function(dat) {
  lmer(
    ext_comp_z ~ exposure_early_z * RACE_early + WAVE +
      sex +
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
    weights = WT2,
    REML = FALSE
  )
})

pool_ext_adol_race1 <- pool(as.mira(fit_ext_adol_race1))

# Externalising: young adulthood
fit_ext_ya_race1 <- lapply(young_adult_list_race1, function(dat) {
  lmer(
    ext_comp_z ~ exposure_early_z * RACE_early + WAVE +
      sex +
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
    weights = WT2,
    REML = FALSE
  )
})

pool_ext_ya_race1 <- pool(as.mira(fit_ext_ya_race1))

# Suicidality: adolescence
fit_sui_adol_race1 <- lapply(adol_list_race1, function(dat) {
  glmer(
    any_suicidality_broad ~ exposure_early_z * RACE_early + WAVE +
      sex +
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
    weights = WT2,
    family = binomial(link = "logit"),
    nAGQ = 0,
    control = glmerControl(optimizer = "bobyqa")
  )
})

pool_sui_adol_race1 <- pool(as.mira(fit_sui_adol_race1))

# Suicidality: young adulthood
fit_sui_ya_race1 <- lapply(young_adult_list_race1, function(dat) {
  glmer(
    any_suicidality_broad ~ exposure_early_z * RACE_early + WAVE +
      sex +
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
    weights = WT2,
    family = binomial(link = "logit"),
    nAGQ = 0,
    control = glmerControl(optimizer = "bobyqa")
  )
})

pool_sui_ya_race1 <- pool(as.mira(fit_sui_ya_race1))


# RQ2: Calculate BIC-approximated Bayes factors ----

# Compare models with and without the cash-transfer interaction
calculate_bf <- function(dat, outcome, family = "gaussian") {
  
  # Same covariates and random intercept in both models
  covariates <- paste(
    "WAVE + sex + AGEGRUP + singleP_baseline +",
    "income_baseline + lowSES_baseline + low_bw +",
    "bmi_baseline + int_baseline + ext_baseline +",
    "(1 | GSMSID)"
  )
  
  formula_null <- as.formula(paste(
    outcome, "~ exposure_early_z + RACE_early +", covariates
  ))
  
  formula_full <- as.formula(paste(
    outcome, "~ exposure_early_z * RACE_early +", covariates
  ))
  
  # Fit models using original estimation settings
  if (family == "gaussian") {
    
    m0 <- lmer(
      formula_null,
      data = dat,
      weights = WT2,
      REML = FALSE
    )
    
    m1 <- lmer(
      formula_full,
      data = dat,
      weights = WT2,
      REML = FALSE
    )
    
  } else if (family == "binomial") {
    
    m0 <- glmer(
      formula_null,
      data = dat,
      weights = WT2,
      family = binomial(link = "logit"),
      nAGQ = 0,
      control = glmerControl(optimizer = "bobyqa")
    )
    
    m1 <- glmer(
      formula_full,
      data = dat,
      weights = WT2,
      family = binomial(link = "logit"),
      nAGQ = 0,
      control = glmerControl(optimizer = "bobyqa")
    )
    
  } else {
    stop("Unknown model family")
  }
  
  # BIC approximation: BF01 favours no interaction
  bic0 <- BIC(m0)
  bic1 <- BIC(m1)
  
  bf01 <- exp((bic1 - bic0) / 2)
  
  data.frame(
    bic_null = bic0,
    bic_full = bic1,
    bf01 = bf01,
    bf10 = 1 / bf01
  )
}

#### Calculate Bayes factors across imputations ####

# Adolescence
bf_int_adol_df <- do.call(rbind, lapply(
  adol_list, calculate_bf,
  outcome = "int_z"
))

bf_ext_adol_df <- do.call(rbind, lapply(
  adol_list, calculate_bf,
  outcome = "ext_comp_z"
))

bf_sui_adol_df <- do.call(rbind, lapply(
  adol_list, calculate_bf,
  outcome = "any_suicidality_broad",
  family = "binomial"
))

# Young adulthood
bf_int_ya_df <- do.call(rbind, lapply(
  young_adult_list, calculate_bf,
  outcome = "int_z"
))

bf_ext_ya_df <- do.call(rbind, lapply(
  young_adult_list, calculate_bf,
  outcome = "ext_comp_z"
))

bf_sui_ya_df <- do.call(rbind, lapply(
  young_adult_list, calculate_bf,
  outcome = "any_suicidality_broad",
  family = "binomial"
))

#### Summarise Bayes factors ####

bf_datasets <- list(
  "Adolescent internalising" = bf_int_adol_df,
  "Adolescent externalising" = bf_ext_adol_df,
  "Adolescent suicidality" = bf_sui_adol_df,
  "Adult internalising" = bf_int_ya_df,
  "Adult externalising" = bf_ext_ya_df,
  "Adult suicidality" = bf_sui_ya_df
)

bf_summary <- do.call(rbind, lapply(names(bf_datasets), function(nm) {
  x <- bf_datasets[[nm]]$bf01
  
  data.frame(
    outcome = nm,
    median_BF01 = median(x),
    min_BF01 = min(x),
    max_BF01 = max(x)
  )
}))

rownames(bf_summary) <- NULL

bf_summary

# RQ2: Create and export Table 3 ----

get_exposure_row <- function(
    pool_obj,
    outcome,
    timepoint,
    group,
    exponentiate = FALSE
) {
  
  s <- summary(
    pool_obj,
    conf.int = TRUE,
    exponentiate = exponentiate
  )
  
  exposure_row <- s |>
    filter(term == "exposure_early_z")
  
  if (nrow(exposure_row) != 1) {
    stop(
      paste(
        "Could not identify exposure_early_z for",
        outcome, timepoint, group
      )
    )
  }
  
  tibble(
    Outcome = outcome,
    Timepoint = timepoint,
    Group = group,
    Estimate = exposure_row$estimate,
    CI_low = exposure_row$`2.5 %`,
    CI_high = exposure_row$`97.5 %`
  )
}

get_interaction_row <- function(
    pool_obj,
    outcome,
    timepoint,
    exponentiate = FALSE
) {
  
  s <- summary(
    pool_obj,
    conf.int = TRUE,
    exponentiate = exponentiate
  )
  
  int_row <- s |>
    filter(
      grepl("exposure_early_z:RACE_early", term) |
        grepl("RACE_early.*:exposure_early_z", term)
    )
  
  if (nrow(int_row) != 1) {
    stop(
      paste(
        "Could not uniquely identify interaction term for",
        outcome, timepoint
      )
    )
  }
  
  tibble(
    Outcome = outcome,
    Timepoint = timepoint,
    Estimate = int_row$estimate,
    CI_low = int_row$`2.5 %`,
    CI_high = int_row$`97.5 %`,
    p_value = int_row$p.value
  )
}

# Specify the pooled models for each outcome and developmental period
model_specs <- tibble(
  Outcome = rep(
    c("Internalising symptoms", "Externalising symptoms", "Suicidality"),
    each = 2
  ),
  Timepoint = rep(c("Adolescence", "Adulthood"), 3),
  pool_noCT = list(
    pool_int_adol, pool_int_ya,
    pool_ext_adol, pool_ext_ya,
    pool_sui_adol, pool_sui_ya
  ),
  pool_CT = list(
    pool_int_adol_race1, pool_int_ya_race1,
    pool_ext_adol_race1, pool_ext_ya_race1,
    pool_sui_adol_race1, pool_sui_ya_race1
  )
)

# Extract group-specific estimates and interaction effects
group_results_long <- bind_rows(lapply(seq_len(nrow(model_specs)), function(i) {
  x <- model_specs[i, ]
  
  bind_rows(
    get_exposure_row(
      x$pool_noCT[[1]], x$Outcome, x$Timepoint,
      "No cash transfer"
    ),
    get_exposure_row(
      x$pool_CT[[1]], x$Outcome, x$Timepoint,
      "Cash transfer"
    )
  )
}))

interaction_results <- bind_rows(lapply(seq_len(nrow(model_specs)), function(i) {
  x <- model_specs[i, ]
  
  get_interaction_row(
    x$pool_noCT[[1]], x$Outcome, x$Timepoint
  )
})) |>
  rename(
    Interaction_estimate = Estimate,
    Interaction_low = CI_low,
    Interaction_high = CI_high,
    Interaction_p = p_value
  )


# Convert group results into table columns

group_results_wide <- group_results_long |>
  mutate(
    Effect = sprintf(
      "%.3f (%.3f, %.3f)",
      Estimate,
      CI_low,
      CI_high
    )
  ) |>
  select(Outcome, Timepoint, Group, Effect) |>
  pivot_wider(
    names_from = Group,
    values_from = Effect
  )

# Bayes factor summary matched to same rows
bf_results <- tibble(
  Outcome = c(
    "Internalising symptoms", "Internalising symptoms",
    "Externalising symptoms", "Externalising symptoms",
    "Suicidality", "Suicidality"
  ),
  Timepoint = c(
    "Adolescence", "Adulthood",
    "Adolescence", "Adulthood",
    "Adolescence", "Adulthood"
  ),
  BF01 = c(
    median(bf_int_adol_df$bf01),
    median(bf_int_ya_df$bf01),
    median(bf_ext_adol_df$bf01),
    median(bf_ext_ya_df$bf01),
    median(bf_sui_adol_df$bf01),
    median(bf_sui_ya_df$bf01)
  )
)

# Combine group-specific effects, interaction p-values and Bayes factors
table_main_groups <- group_results_wide |>
  left_join(
    interaction_results,
    by = c("Outcome", "Timepoint")
  ) |>
  left_join(
    bf_results,
    by = c("Outcome", "Timepoint")
  ) |>
  mutate(
    Interaction = sprintf(
      "%.3f (%.3f, %.3f)",
      Interaction_estimate,
      Interaction_low,
      Interaction_high
    ),
    
    # Format interaction p-values
    P_value = ifelse(
      Interaction_p < 0.001,
      "<0.001",
      sprintf("%.3f", Interaction_p)
    ),
    
    BF01 = round(BF01, 1)
  ) |>
  select(
    Outcome,
    Timepoint,
    `No cash transfer`,
    `Cash transfer`,
    Interaction,
    P_value,
    BF01
  )

table_main_groups

# Combine into flextable

ft <- flextable(table_main_groups) |>
  set_header_labels(
    Outcome = "Outcome",
    Timepoint = "Developmental period",
    `No cash transfer` = "No cash transfer: association (95% CI)",
    `Cash transfer` = "Cash transfer: association (95% CI)",
    Interaction = "Interaction effect (95% CI)",
    P_value = "p-value",
    BF01 = "Bayes factor (null vs interaction)"
  ) |>
  set_table_properties(layout = "fixed") |>
  width(j = "Outcome", width = 2.3) |>
  width(j = "Timepoint", width = 1.7) |>
  width(j = "No cash transfer", width = 2.5) |>
  width(j = "Cash transfer", width = 2.5) |>
  width(j = "Interaction", width = 2.5) |>
  width(j = "P_value", width = 1.0) |>
  width(j = "BF01", width = 1.5) |>
  align(
    j = c("Outcome", "Timepoint"),
    align = "left"
  ) |>
  align(
    j = c(
      "No cash transfer",
      "Cash transfer",
      "P_value",
      "Interaction",
      "BF01"
    ),
    align = "center"
  ) |>
  theme_booktabs()

save_as_docx(
  "Table 3" = ft,
  path = file.path(tables_dir, "Table3.docx")
)

# Heterogeneity in cash-transfer effects: cash amount and duration ----


adol_list <- lapply(adol_list, function(dat) {
  dat %>%
    mutate(
      AIpar_num = as.numeric(as.character(AIpar)),
      duration_num = case_when(
        AGEGRUP == 1 ~ 2,  # youngest: longest exposure
        AGEGRUP == 2 ~ 1,  # middle
        AGEGRUP == 3 ~ 0   # oldest: shortest exposure
      )
    )
})

young_adult_list <- lapply(young_adult_list, function(dat) {
  dat %>%
    mutate(
      AIpar_num = as.numeric(as.character(AIpar)),
      duration_num = case_when(
        AGEGRUP == 1 ~ 2,  # youngest: longest exposure
        AGEGRUP == 2 ~ 1,  # middle
        AGEGRUP == 3 ~ 0   # oldest: shortest exposure
      )
    )
})


# Fit models across imputations and pool estimates
fit_heterogeneity <- function(data_list, outcome, moderator,
                              binary = FALSE) {
  
  # Cash amount models adjust for cohort;
  # duration models include cohort in a three-way interaction
  if (moderator == "AIpar_num") {
    interaction <- "exposure_early_z * AIpar_num"
    adjustment <- "AGEGRUP"
  } else if (moderator == "duration_num") {
    interaction <- "exposure_early_z * RACE_early * duration_num"
    adjustment <- NULL
  } else {
    stop("Unknown moderator")
  }
  
  covariates <- paste(
    "WAVE + sex +",
    if (!is.null(adjustment)) paste0(adjustment, " +") else "",
    "singleP_baseline + income_baseline + lowSES_baseline +",
    "low_bw + bmi_baseline + int_baseline + ext_baseline +",
    "(1 | GSMSID)"
  )
  
  model_formula <- as.formula(
    paste(outcome, "~", interaction, "+", covariates)
  )
  
  fits <- lapply(data_list, function(dat) {
    
    if (binary) {
      glmer(
        model_formula,
        data = dat,
        weights = WT2,
        family = binomial(link = "logit"),
        nAGQ = 0,
        control = glmerControl(optimizer = "bobyqa")
      )
    } else {
      lmer(
        model_formula,
        data = dat,
        weights = WT2,
        REML = FALSE
      )
    }
  })
  
  pool(as.mira(fits))
}


# Heterogeneity: Cash-transfer amount ----

# Adolescence
pool_int_adol_cashdose <- fit_heterogeneity(
  adol_list, "int_z", "AIpar_num"
)

pool_ext_adol_cashdose <- fit_heterogeneity(
  adol_list, "ext_comp_z", "AIpar_num"
)

pool_sui_adol_cashdose <- fit_heterogeneity(
  adol_list, "any_suicidality_broad", "AIpar_num",
  binary = TRUE
)

# Young adulthood
pool_int_ya_cashdose <- fit_heterogeneity(
  young_adult_list, "int_z", "AIpar_num"
)

pool_ext_ya_cashdose <- fit_heterogeneity(
  young_adult_list, "ext_comp_z", "AIpar_num"
)

pool_sui_ya_cashdose <- fit_heterogeneity(
  young_adult_list, "any_suicidality_broad", "AIpar_num",
  binary = TRUE
)

# Heterogeneity: Duration of exposure ----

# Adolescence
pool_int_adol_duration <- fit_heterogeneity(
  adol_list, "int_z", "duration_num"
)

pool_ext_adol_duration <- fit_heterogeneity(
  adol_list, "ext_comp_z", "duration_num"
)

pool_sui_adol_duration <- fit_heterogeneity(
  adol_list, "any_suicidality_broad", "duration_num",
  binary = TRUE
)

# Young adulthood
pool_int_ya_duration <- fit_heterogeneity(
  young_adult_list, "int_z", "duration_num"
)

pool_ext_ya_duration <- fit_heterogeneity(
  young_adult_list, "ext_comp_z", "duration_num"
)

pool_sui_ya_duration <- fit_heterogeneity(
  young_adult_list, "any_suicidality_broad", "duration_num",
  binary = TRUE
)

# Heterogeneity: Create and export Table ----

# Extract one interaction term per model

get_one_interaction <- function(pool_obj, outcome, timepoint, analysis,
                                term_name, exponentiate = FALSE) {
  
  s <- summary(pool_obj, conf.int = TRUE, exponentiate = exponentiate)
  
  out <- s %>%
    filter(term == term_name)
  
  if (nrow(out) != 1) {
    stop(paste("Could not uniquely identify interaction term for", outcome, timepoint, analysis))
  }
  
  tibble(
    Outcome = outcome,
    Timepoint = timepoint,
    Analysis = analysis,
    Estimate = out$estimate,
    CI_low = out$`2.5 %`,
    CI_high = out$`97.5 %`
  )
}


# Specify models and interaction terms for Table 4
robustness_specs <- tibble(
  Outcome = rep(
    c("Internalising symptoms", "Externalising symptoms", "Suicidality"),
    each = 4
  ),
  Timepoint = rep(
    c("Adolescence", "Adulthood", "Adolescence", "Adulthood"),
    3
  ),
  Analysis = rep(
    c("Cash amount", "Cash amount",
      "Duration of cash transfer", "Duration of cash transfer"),
    3
  ),
  Model = list(
    pool_int_adol_cashdose, pool_int_ya_cashdose,
    pool_int_adol_duration, pool_int_ya_duration,
    pool_ext_adol_cashdose, pool_ext_ya_cashdose,
    pool_ext_adol_duration, pool_ext_ya_duration,
    pool_sui_adol_cashdose, pool_sui_ya_cashdose,
    pool_sui_adol_duration, pool_sui_ya_duration
  )
)

# Extract the relevant interaction coefficient from each model
table_robustness_gsms <- bind_rows(
  lapply(seq_len(nrow(robustness_specs)), function(i) {
    
    x <- robustness_specs[i, ]
    
    term_name <- if (x$Analysis == "Cash amount") {
      "exposure_early_z:AIpar_num"
    } else {
      "exposure_early_z:RACE_early1:duration_num"
    }
    
    get_one_interaction(
      pool_obj = x$Model[[1]],
      outcome = x$Outcome,
      timepoint = x$Timepoint,
      analysis = x$Analysis,
      term_name = term_name,
      exponentiate = (x$Outcome == "Suicidality")
    )
  })
) |>
  mutate(
    Effect = if_else(
      Outcome == "Suicidality",
      sprintf("%.2f (%.2f to %.2f)", Estimate, CI_low, CI_high),
      sprintf("%.3f (%.3f to %.3f)", Estimate, CI_low, CI_high)
    )
  ) |>
  select(Outcome, Timepoint, Analysis, Effect)

table_robustness_gsms

ft_robust <- flextable(table_robustness_gsms)

ft_robust <- ft_robust |>
  set_header_labels(
    Outcome = "Outcome",
    Timepoint = "Developmental period",
    Analysis = "Robustness analysis",
    Effect = "Interaction effect (95% CI)"
  ) |>
  set_table_properties(layout = "fixed") |>
  width(j = "Outcome", width = 2.7) |>
  width(j = "Timepoint", width = 1.8) |>
  width(j = "Analysis", width = 2.3) |>
  width(j = "Effect", width = 2.7) |>
  align(j = c("Outcome", "Timepoint", "Analysis"), align = "left") |>
  align(j = "Effect", align = "center") |>
  theme_booktabs()

ft_robust
save_as_docx(
  "Table 4" = ft_robust,
  path = file.path(tables_dir, "Table4.docx")
)
