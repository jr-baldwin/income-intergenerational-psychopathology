# ============================================================
# 06: Baby's First Years (BFY) analyses
# ============================================================
#
# This script conducts analyses of the Baby's First Years
# randomised cash-transfer study, examining:
#
# 1. Associations between baseline maternal depressive
#    symptoms and child internalising and externalising symptoms.
#
# 2. Moderation of these associations by cash-transfer group.
#
# 3. Sensitivity analyses examining differential attrition
#    and inverse-probability weighting.
#
# 4. Exploratory analyses of concurrent maternal depressive
#    symptoms and subsequent maternal depressive symptoms.
#
# The script also produces descriptive and results tables.
#
# ============================================================

# Load required packages ----
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(geepack)
library(broom)
library(broom.mixed)
library(flextable)
library(officer)

# Set up file paths ----
source("config.R")

bfy_data_dir <- file.path(gsms_dir, "restricted_data", "BFY")
tables_dir <- file.path(gsms_dir, "results", "tables")

# Load age-4 BFY data ----
load(file.path(bfy_data_dir, "37871-0009-Data.rda"))

babyfy4 <- da37871.0009
rm(da37871.0009)

# Merge with age 3 data
load(file.path(bfy_data_dir, "37871-0004-Data.rda"))

babyfy3 <- da37871.0004
rm(da37871.0004)

babyfy_43merged <- merge(babyfy4, babyfy3,by="PUBLICSAMPLEID")

# Make internalising and externalising outcomes
# Age 3
babyfy_43merged$internalizing_raw_a3 <-
  babyfy_43merged$CCBCLANXDEPRAWA3 +
  babyfy_43merged$CCBCLEMOREACTIVERAWA3

babyfy_43merged$externalizing_raw_a3 <-
  babyfy_43merged$CCBCLAGGBEHAVIORRAWA3 +
  babyfy_43merged$CCBCLATTENTIONRAWA3

# Age 4
babyfy_43merged$internalizing_raw_a4 <-
  babyfy_43merged$CCBCLANXDEPRAWA4 +
  babyfy_43merged$CCBCLEMOREACTIVERAWA4

babyfy_43merged$externalizing_raw_a4 <-
  babyfy_43merged$CCBCLAGGBEHAVIORRAWA4 +
  babyfy_43merged$CCBCLATTENTIONRAWA4

# Merge with age 1 data

load(file.path(bfy_data_dir, "37871-0002-Data.rda"))

babyfy1 <- da37871.0002
rm(da37871.0002)

babyfy_merged <- merge(babyfy_43merged, babyfy1,by="PUBLICSAMPLEID")

# Download and merge with baseline

load(file.path(bfy_data_dir, "37871-0001-Data.rda"))
babyfy0 <- da37871.0001
rm(da37871.0001)

babyfy_merged <- merge(babyfy_merged, babyfy0,by="PUBLICSAMPLEID")

# Remove intermediate datasets no longer needed
rm(babyfy4, babyfy3, babyfy1, babyfy0, babyfy_43merged)

# Reshape into a long dataset
babyfy_outcomes_long <- babyfy_merged %>%
  select(
    PUBLICSAMPLEID,
    MDEPRESSIONA0, # maternal depression baseline
    MPHQ8A3, # maternal depression age 3
    MPHQ8A4, # maternal depression age 4
    # Covariates
    HHREVISEDINCOMEA0,
    MEDLEVELA0,
    HHBIODADA0,
    CFEMALEA0, 
    CWEIGHTOZA0,
    MPHQ8A1,
    TREAT, # Cash-gift group
    SITE, # site
    internalizing_raw_a3,
    internalizing_raw_a4,
    externalizing_raw_a3,
    externalizing_raw_a4
  ) %>%
  # Rename maternal depression variables to match outcome naming
  rename(
    maternal_depression_a3 = MPHQ8A3,
    maternal_depression_a4 = MPHQ8A4
  ) %>%
  
  pivot_longer(
    cols = c(
      internalizing_raw_a3,
      internalizing_raw_a4,
      externalizing_raw_a3,
      externalizing_raw_a4,
      maternal_depression_a3,
      maternal_depression_a4
    ),
    names_to = c(".value", "age_wave"),
    names_pattern = "(internalizing_raw|externalizing_raw|maternal_depression)_a(3|4)"
  ) %>%
  
  mutate(
    age_wave = factor(
      age_wave,
      levels = c("3", "4"),
      labels = c("Age 3", "Age 4")
    )
  )

head(babyfy_outcomes_long)

# Main analyses: maternal depression and child psychopathology ----

## No covariates
# Internalising problems
rq1_int <- lmer(scale(internalizing_raw) ~ scale(MDEPRESSIONA0) + TREAT + 
                  age_wave + factor(SITE) + (1 | PUBLICSAMPLEID), data=babyfy_outcomes_long)
broom.mixed::tidy(rq1_int, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
  filter(grepl("MDEPRESSIONA0", term)) %>%
  select(estimate, conf.low, conf.high, p.value) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ round(.x, 2)), p.value = round(p.value, 3))

# Externalising problems
rq1_ext <- lmer(scale(externalizing_raw) ~ scale(MDEPRESSIONA0) + TREAT + 
                  age_wave + factor(SITE) + (1 | PUBLICSAMPLEID), data=babyfy_outcomes_long)
broom.mixed::tidy(rq1_ext, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
  filter(grepl("MDEPRESSIONA0", term)) %>%
  select(estimate, conf.low, conf.high, p.value) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ round(.x, 2)), p.value = round(p.value, 3))

## Including covariates 
babyfy_outcomes_long <- babyfy_outcomes_long %>%
  mutate(
    income_z = as.numeric(scale(HHREVISEDINCOMEA0)),
    birthweight_z = as.numeric(scale(CWEIGHTOZA0))
  )
# Internalising problems
rq1_int_cov <- lmer(scale(internalizing_raw) ~ scale(MDEPRESSIONA0) + TREAT + 
                      income_z + MEDLEVELA0 + HHBIODADA0 + CFEMALEA0 + birthweight_z +
                      age_wave + factor(SITE) + (1 | PUBLICSAMPLEID), data=babyfy_outcomes_long)
broom.mixed::tidy(rq1_int_cov, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
  filter(grepl("MDEPRESSIONA0", term)) %>%
  select(estimate, conf.low, conf.high, p.value) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ round(.x, 2)), p.value = round(p.value, 3))

# Externalising problems
rq1_ext_cov <- lmer(scale(externalizing_raw) ~ scale(MDEPRESSIONA0) + TREAT + 
                      income_z + MEDLEVELA0 + HHBIODADA0 + CFEMALEA0 + birthweight_z +
                      age_wave + factor(SITE) + (1 | PUBLICSAMPLEID), data=babyfy_outcomes_long)
broom.mixed::tidy(rq1_ext_cov, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
  filter(grepl("MDEPRESSIONA0", term)) %>%
  select(estimate, conf.low, conf.high, p.value) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ round(.x, 2)), p.value = round(p.value, 3))

# Cash-transfer moderation analyses ----

### Standardise variables ------------------------------------------------

# Standardise variables ----

babyfy_outcomes_long <- babyfy_outcomes_long %>%
  mutate(
    depression_z = as.numeric(scale(MDEPRESSIONA0)),
    income_z = as.numeric(scale(HHREVISEDINCOMEA0)),
    birthweight_z = as.numeric(scale(CWEIGHTOZA0))
  ) %>%
  group_by(age_wave) %>%
  mutate(
    internalizing_z = as.numeric(scale(internalizing_raw)),
    externalizing_z = as.numeric(scale(externalizing_raw))
  ) %>%
  ungroup()

# Fit primary and covariate-adjusted moderation models ----

# Primary models
rq2_int <- lmer(
  internalizing_z ~ depression_z * TREAT +
    age_wave + factor(SITE) + (1 | PUBLICSAMPLEID),
  data = babyfy_outcomes_long
)

rq2_ext <- lmer(
  externalizing_z ~ depression_z * TREAT +
    age_wave + factor(SITE) + (1 | PUBLICSAMPLEID),
  data = babyfy_outcomes_long
)

# Covariate-adjusted models
rq2_int_cov <- lmer(
  internalizing_z ~ depression_z * TREAT +
    income_z + MEDLEVELA0 + HHBIODADA0 +
    CFEMALEA0 + birthweight_z +
    age_wave + factor(SITE) + (1 | PUBLICSAMPLEID),
  data = babyfy_outcomes_long
)

rq2_ext_cov <- lmer(
  externalizing_z ~ depression_z * TREAT +
    income_z + MEDLEVELA0 + HHBIODADA0 +
    CFEMALEA0 + birthweight_z +
    age_wave + factor(SITE) + (1 | PUBLICSAMPLEID),
  data = babyfy_outcomes_long
)

summary(rq2_int)
summary(rq2_ext)
summary(rq2_int_cov)
summary(rq2_ext_cov)

# Calculate BIC-approximated Bayes factors ----
get_BF01 <- function(outcome, covariates = FALSE) {
  
  covs <- if (covariates) {
    " + HHREVISEDINCOMEA0 + MEDLEVELA0 + HHBIODADA0 + CFEMALEA0 + CWEIGHTOZA0"
  } else {
    ""
  }
  
  formula0 <- as.formula(
    paste0(
      outcome,
      " ~ depression_z + TREAT",
      covs,
      " + age_wave + factor(SITE) + (1 | PUBLICSAMPLEID)"
    )
  )
  
  formula1 <- as.formula(
    paste0(
      outcome,
      " ~ depression_z * TREAT",
      covs,
      " + age_wave + factor(SITE) + (1 | PUBLICSAMPLEID)"
    )
  )
  
  m0 <- lmer(
    formula0,
    data = babyfy_outcomes_long,
    REML = FALSE
  )
  
  m1 <- lmer(
    formula1,
    data = babyfy_outcomes_long,
    REML = FALSE
  )
  
  exp((BIC(m1) - BIC(m0)) / 2)
}

BF01_int     <- get_BF01("internalizing_z")
BF01_int_cov <- get_BF01("internalizing_z", covariates = TRUE)

BF01_ext     <- get_BF01("externalizing_z")
BF01_ext_cov <- get_BF01("externalizing_z", covariates = TRUE)

# Create cash-transfer moderation results table ----

make_mod_row <- function(model, outcome, model_name, BF = NA) {
  
  # Maternal depression association within each treatment group
  trends <- emtrends(
    model,
    ~ TREAT,
    var = "depression_z"
  ) %>%
    confint() %>%
    as.data.frame()
  
  # Interaction estimate, CI and p-value
  coefs <- coef(summary(model))
  
  interaction_name <- grep(
    "depression_z:TREAT",
    rownames(coefs),
    value = TRUE
  )
  
  interaction_ci <- confint(
    model,
    parm = interaction_name,
    method = "Wald"
  )
  
  interaction_est <- coefs[interaction_name, "Estimate"]
  interaction_p   <- coefs[interaction_name, "Pr(>|t|)"]
  
  tibble(
    Outcome = outcome,
    Model = model_name,
    
    `Low-cash group: association (95% CI)` =
      sprintf(
        "%.3f (%.3f, %.3f)",
        trends$depression_z.trend[1],
        trends$lower.CL[1],
        trends$upper.CL[1]
      ),
    
    `High-cash group: association (95% CI)` =
      sprintf(
        "%.3f (%.3f, %.3f)",
        trends$depression_z.trend[2],
        trends$lower.CL[2],
        trends$upper.CL[2]
      ),
    
    `Interaction effect (95% CI)` =
      sprintf(
        "%.3f (%.3f, %.3f)",
        interaction_est,
        interaction_ci[1],
        interaction_ci[2]
      ),
    
    `Interaction P value` =
      ifelse(
        interaction_p < .001,
        "<.001",
        sprintf("%.3f", interaction_p)
      ),
    
    `BF01` =
      ifelse(
        is.na(BF),
        "—",
        sprintf("%.1f", BF)
      )
  )
}

## Create table
table_bfy <- bind_rows(
  make_mod_row(
    rq2_int,
    "Internalising symptoms",
    "Primary",
    BF01_int
  ),
  make_mod_row(
    rq2_int_cov,
    "Internalising symptoms",
    "Covariate-adjusted",
    BF01_int_cov
  ),
  make_mod_row(
    rq2_ext,
    "Externalising symptoms",
    "Primary",
    BF01_ext
  ),
  make_mod_row(
    rq2_ext_cov,
    "Externalising symptoms",
    "Covariate-adjusted",
    BF01_ext_cov
  )
)

table_bfy <- table_bfy %>%
  rename(
    `Low-cash group, β (95% CI)` =
      `Low-cash group: association (95% CI)`,
    `High-cash group, β (95% CI)` =
      `High-cash group: association (95% CI)`,
    `Interaction, β (95% CI)` =
      `Interaction effect (95% CI)`,
    `P value` = `Interaction P value`
  )
table_bfy

# Format and export moderation results table ----

ft <- flextable(table_bfy) %>%
  autofit() %>%
  theme_booktabs() %>%
  align(align = "center", part = "all") %>%
  align(j = c("Outcome", "Model"), align = "left", part = "all")

doc <- read_docx() %>%
  body_add_flextable(ft)

print(
  doc,
  target = file.path(tables_dir, "BFY_moderation_table.docx")
)

# Sensitivity analyses: attrition and inverse-probability weighting ----

# Create outcome-availability indicators
babyfy_merged <- babyfy_merged %>%
  mutate(
    depression_z = as.numeric(scale(MDEPRESSIONA0)),
    
    observed_int_a3 = as.integer(!is.na(internalizing_raw_a3)),
    observed_int_a4 = as.integer(!is.na(internalizing_raw_a4)),
    observed_ext_a3 = as.integer(!is.na(externalizing_raw_a3)),
    observed_ext_a4 = as.integer(!is.na(externalizing_raw_a4))
  )


# Fit observation models ------------------------------------------------

attr_int_a3 <- glm(
  observed_int_a3 ~ depression_z * TREAT + factor(SITE),
  family = binomial,
  data = babyfy_merged
)

attr_int_a4 <- glm(
  observed_int_a4 ~ depression_z * TREAT + factor(SITE),
  family = binomial,
  data = babyfy_merged
)

attr_ext_a3 <- glm(
  observed_ext_a3 ~ depression_z * TREAT + factor(SITE),
  family = binomial,
  data = babyfy_merged
)

attr_ext_a4 <- glm(
  observed_ext_a4 ~ depression_z * TREAT + factor(SITE),
  family = binomial,
  data = babyfy_merged
)

summary(attr_int_a3)
summary(attr_int_a4)
summary(attr_ext_a3)
summary(attr_ext_a4)

# Predicted probability of observation ---------------------------------

babyfy_merged <- babyfy_merged %>%
  mutate(
    p_int_a3 = predict(attr_int_a3, newdata = babyfy_merged, type = "response"),
    p_int_a4 = predict(attr_int_a4, newdata = babyfy_merged, type = "response"),
    p_ext_a3 = predict(attr_ext_a3, newdata = babyfy_merged, type = "response"),
    p_ext_a4 = predict(attr_ext_a4, newdata = babyfy_merged, type = "response")
  )


# Marginal probabilities of observation --------------------------------

eligible <- !is.na(babyfy_merged$MDEPRESSIONA0)

pbar_int_a3 <- mean(babyfy_merged$observed_int_a3[eligible])
pbar_int_a4 <- mean(babyfy_merged$observed_int_a4[eligible])

pbar_ext_a3 <- mean(babyfy_merged$observed_ext_a3[eligible])
pbar_ext_a4 <- mean(babyfy_merged$observed_ext_a4[eligible])


# Stabilized inverse-probability-of-observation weights -----------------

babyfy_merged <- babyfy_merged %>%
  mutate(
    sw_int_a3 = ifelse(
      observed_int_a3 == 1,
      pbar_int_a3 / p_int_a3,
      NA
    ),
    
    sw_int_a4 = ifelse(
      observed_int_a4 == 1,
      pbar_int_a4 / p_int_a4,
      NA
    ),
    
    sw_ext_a3 = ifelse(
      observed_ext_a3 == 1,
      pbar_ext_a3 / p_ext_a3,
      NA
    ),
    
    sw_ext_a4 = ifelse(
      observed_ext_a4 == 1,
      pbar_ext_a4 / p_ext_a4,
      NA
    )
  )

# Create separate long dataset for IPW analyses ----
babyfy_ipw_long <- babyfy_merged %>%
  select(
    PUBLICSAMPLEID,
    MDEPRESSIONA0,
    HHREVISEDINCOMEA0,
    MEDLEVELA0,
    HHBIODADA0,
    CFEMALEA0, 
    CWEIGHTOZA0,
    MPHQ8A1,
    TREAT,
    SITE,
    sw_int_a3,
    sw_int_a4,
    sw_ext_a3,
    sw_ext_a4,
    internalizing_raw_a3,
    internalizing_raw_a4,
    externalizing_raw_a3,
    externalizing_raw_a4
  ) %>%
  pivot_longer(
    cols = c(
      internalizing_raw_a3,
      internalizing_raw_a4,
      externalizing_raw_a3,
      externalizing_raw_a4
    ),
    names_to = c(".value", "age_wave"),
    names_pattern = "(internalizing_raw|externalizing_raw)_a(3|4)"
  ) %>%
  mutate(
    age_wave = factor(
      age_wave,
      levels = c("3", "4"),
      labels = c("Age 3", "Age 4")
    ),
    
    # Appropriate weight for outcome and wave
    ipw_int = case_when(
      age_wave == "Age 3" ~ sw_int_a3,
      age_wave == "Age 4" ~ sw_int_a4
    ),
    
    ipw_ext = case_when(
      age_wave == "Age 3" ~ sw_ext_a3,
      age_wave == "Age 4" ~ sw_ext_a4
    )
  ) %>%
  
  # Standardise each outcome separately within wave
  group_by(age_wave) %>%
  mutate(
    internalizing_z = as.numeric(scale(internalizing_raw)),
    externalizing_z = as.numeric(scale(externalizing_raw))
  ) %>%
  ungroup() %>%
  
  # Standardise baseline depression once across participants
  mutate(
    depression_z = as.numeric(scale(MDEPRESSIONA0))
  )

### Check all four weights
babyfy_merged %>%
  summarise(
    across(
      c(sw_int_a3, sw_int_a4, sw_ext_a3, sw_ext_a4),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        min  = ~min(.x, na.rm = TRUE),
        max  = ~max(.x, na.rm = TRUE)
      )
    )
  )


### Fit the two IPW models

ipw_int <- babyfy_ipw_long %>%
  filter(
    !is.na(internalizing_z),
    !is.na(depression_z),
    !is.na(ipw_int)
  )

ipw_ext <- babyfy_ipw_long %>%
  filter(
    !is.na(externalizing_z),
    !is.na(depression_z),
    !is.na(ipw_ext)
  )

ipw_model_int <- geeglm(
  internalizing_z ~ depression_z * TREAT +
    age_wave + factor(SITE),
  id = PUBLICSAMPLEID,  
  weights = ipw_int,
  corstr = "exchangeable",
  family = gaussian,
  data = ipw_int
)

ipw_model_ext <- geeglm(
  externalizing_z ~ depression_z * TREAT +
    age_wave + factor(SITE),
  id = PUBLICSAMPLEID,
  weights = ipw_ext,
  corstr = "exchangeable",
  family = gaussian,
  data = ipw_ext
)

summary(ipw_model_int)
summary(ipw_model_ext)

# Create attrition and IPW sensitivity results table ----

get_attrition_result <- function(model, outcome, age) {
  
  x <- broom::tidy(model, conf.int = TRUE) %>%
    filter(grepl("depression_z.*TREAT|TREAT.*depression_z", term))
  
  tibble(
    Outcome = outcome,
    Assessment = age,
    N = nobs(model),
    `Interaction, β (95% CI)` = sprintf(
      "%.3f (%.3f, %.3f)",
      x$estimate,
      x$conf.low,
      x$conf.high
    ),
    `P value` = ifelse(
      x$p.value < .001,
      "<.001",
      sprintf("%.3f", x$p.value)
    )
  )
}

panel_a <- bind_rows(
  get_attrition_result(attr_int_a3, "Internalising", "Age 3"),
  get_attrition_result(attr_int_a4, "Internalising", "Age 4"),
  get_attrition_result(attr_ext_a3, "Externalising", "Age 3"),
  get_attrition_result(attr_ext_a4, "Externalising", "Age 4")
)

## Add N to panel B

get_primary_result <- function(model, outcome) {
  
  coefs <- coef(summary(model))
  
  term <- grep(
    "depression_z:TREAT",
    rownames(coefs),
    value = TRUE
  )
  
  est <- coefs[term, "Estimate"]
  se  <- coefs[term, "Std. Error"]
  p   <- coefs[term, "Pr(>|t|)"]
  
  tibble(
    Outcome = outcome,
    Analysis = "Primary",
    N = length(unique(model.frame(model)$PUBLICSAMPLEID)),
    `Interaction, β (95% CI)` = sprintf(
      "%.3f (%.3f, %.3f)",
      est,
      est - 1.96 * se,
      est + 1.96 * se
    ),
    `P value` = ifelse(
      p < .001, "<.001", sprintf("%.3f", p)
    )
  )
}

## For GEE
get_ipw_result <- function(model, outcome) {
  
  coefs <- summary(model)$coefficients
  
  term <- grep(
    "depression_z:TREAT",
    rownames(coefs),
    value = TRUE
  )
  
  est <- coefs[term, "Estimate"]
  se  <- coefs[term, "Std.err"]
  p   <- coefs[term, "Pr(>|W|)"]
  
  tibble(
    Outcome = outcome,
    Analysis = "IPW",
    N = length(unique(model$id)),
    `Interaction, β (95% CI)` = sprintf(
      "%.3f (%.3f, %.3f)",
      est,
      est - 1.96 * se,
      est + 1.96 * se
    ),
    `P value` = ifelse(
      p < .001, "<.001", sprintf("%.3f", p)
    )
  )
}

panel_b <- bind_rows(
  get_primary_result(rq2_int, "Internalising"),
  get_ipw_result(ipw_model_int, "Internalising"),
  get_primary_result(rq2_ext, "Externalising"),
  get_ipw_result(ipw_model_ext, "Externalising")
)

# Format and export attrition/IPW supplementary table ----

ft_a <- flextable(panel_a) %>%
  set_header_labels(
    Outcome = "Outcome",
    Assessment = "Assessment",
    N = "N",
    `Interaction, β (95% CI)` = "Interaction, β (95% CI)",
    `P value` = "P value"
  ) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  align(j = c("Outcome", "Assessment"), align = "left", part = "all") %>%
  align(j = c("N", "Interaction, β (95% CI)", "P value"),
        align = "center", part = "all") %>%
  autofit()

ft_b <- flextable(panel_b) %>%
  set_header_labels(
    Outcome = "Outcome",
    Analysis = "Analysis",
    N = "N",
    `Interaction, β (95% CI)` = "Interaction, β (95% CI)",
    `P value` = "P value"
  ) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  align(j = c("Outcome", "Analysis"), align = "left", part = "all") %>%
  align(j = c("N", "Interaction, β (95% CI)", "P value"),
        align = "center", part = "all") %>%
  autofit()

doc <- read_docx() %>%
  
  body_add_par(
    "Supplementary Table X. Sensitivity analyses examining differential outcome attrition in Baby’s First Years",
    style = "Normal"
  ) %>%
  
  body_add_par(
    "Panel A. Predictors of outcome availability",
    style = "Normal"
  ) %>%
  body_add_flextable(ft_a) %>%
  
  body_add_par("") %>%
  
  body_add_par(
    "Panel B. Inverse-probability-of-observation-weighted sensitivity analyses",
    style = "Normal"
  ) %>%
  body_add_flextable(ft_b) %>%
  
  body_add_par("") %>%
  
  body_add_par(
    paste0(
      "Note. Panel A presents logistic regression models predicting outcome availability ",
      "from baseline maternal depressive symptoms, randomized cash-transfer group, their ",
      "interaction, and study site; β coefficients are on the log-odds scale. ",
      "Panel B compares the primary moderation analyses with post-hoc sensitivity analyses ",
      "using stabilized inverse-probability-of-observation weights estimated separately ",
      "for each outcome and assessment wave. β coefficients in Panel B are standardized. ",
      "CI = confidence interval; IPW = inverse-probability weighting."
    ),
    style = "Normal"
  )

print(
  doc,
  target = file.path(tables_dir, "BFY_attrition_IPW_table.docx")
)

# Descriptive statistics ----

### Get analytic sample size to report descriptively
# Participant IDs used in each primary moderation model
analytic_int <- babyfy_outcomes_long %>%
  filter(
    !is.na(internalizing_z),
    !is.na(depression_z),
    !is.na(TREAT),
    !is.na(age_wave),
    !is.na(SITE),
    !is.na(PUBLICSAMPLEID)
  ) %>%
  distinct(PUBLICSAMPLEID, TREAT)

nrow(analytic_int)

### Get N low and high cash
analytic_int %>%
  count(TREAT) %>%
  mutate(
    percent = 100 * n / sum(n)
  )

### Initial step: subset to 937 included in main models ###
# Rows actually used in the model
model_rows <- as.numeric(rownames(model.frame(rq2_int)))

# IDs of participants actually included
analytic_ids <- babyfy_outcomes_long %>%
  slice(model_rows) %>%
  distinct(PUBLICSAMPLEID)

nrow(analytic_ids)

## Create the descriptive dataset
babyfy_analytic <- babyfy_merged %>%
  semi_join(
    analytic_ids,
    by = "PUBLICSAMPLEID"
  )

nrow(babyfy_analytic)

## Check treatment allocation
babyfy_analytic %>%
  count(TREAT) %>%
  mutate(
    percent = 100 * n / sum(n)
  )

# ============================================================
# 1. Create clean variables for the descriptive table
# ============================================================

babyfy_analytic <- babyfy_analytic %>%
  mutate(
    
    # Treatment: 0 = low cash, 1 = high cash
    treat_num = case_when(
      TREAT == "(0) Low cash gift group"  ~ 0,
      TREAT == "(1) High cash gift group" ~ 1,
      TRUE ~ NA_real_
    ),
    
    # Biological father in household: 0 = no, 1 = yes
    biodad_num = case_when(
      HHBIODADA0 == "(0) No"  ~ 0,
      HHBIODADA0 == "(1) Yes" ~ 1,
      TRUE ~ NA_real_
    ),
    
    # Child sex: 0 = male, 1 = female
    female_num = case_when(
      CFEMALEA0 == "(0) No"  ~ 0,
      CFEMALEA0 == "(1) Yes" ~ 1,
      TRUE ~ NA_real_
    ),
    
    # Convert birth weight from ounces to grams
    birthweight_g = CWEIGHTOZA0 * 28.3495
  )


# ============================================================
# 2. Helper functions
# ============================================================

# Mean (SD), with N
mean_sd_n <- function(x, digits = 2) {
  
  fmt <- paste0("%.", digits, "f (%.", digits, "f); n=%d")
  
  sprintf(
    fmt,
    mean(x, na.rm = TRUE),
    sd(x, na.rm = TRUE),
    sum(!is.na(x))
  )
}


# n (%)
n_pct <- function(x, value = 1) {
  
  denom <- sum(!is.na(x))
  n <- sum(x == value, na.rm = TRUE)
  
  sprintf(
    "%d (%.1f%%)",
    n,
    100 * n / denom
  )
}


# Format P values
format_p <- function(p) {
  
  ifelse(
    p < .001,
    "<.001",
    sprintf("%.3f", p)
  )
}


# ============================================================
# 3. Function for continuous variables
# ============================================================

continuous_row <- function(data, variable, label, digits = 2) {
  
  p <- t.test(
    data[[variable]] ~ data$treat_num
  )$p.value
  
  tibble(
    Characteristic = label,
    
    Overall = mean_sd_n(
      data[[variable]],
      digits
    ),
    
    `Low cash` = mean_sd_n(
      data[[variable]][data$treat_num == 0],
      digits
    ),
    
    `High cash` = mean_sd_n(
      data[[variable]][data$treat_num == 1],
      digits
    ),
    
    `P value` = format_p(p)
  )
}


# ============================================================
# 4. Function for binary variables
# ============================================================

binary_row <- function(data, variable, label) {
  
  p <- chisq.test(
    table(
      data[[variable]],
      data$treat_num
    )
  )$p.value
  
  tibble(
    Characteristic = label,
    
    Overall = n_pct(
      data[[variable]],
      1
    ),
    
    `Low cash` = n_pct(
      data[[variable]][data$treat_num == 0],
      1
    ),
    
    `High cash` = n_pct(
      data[[variable]][data$treat_num == 1],
      1
    ),
    
    `P value` = format_p(p)
  )
}


# ============================================================
# 5. Maternal education
# ============================================================

education_p <- chisq.test(
  table(
    babyfy_analytic$MEDLEVELA0,
    babyfy_analytic$treat_num
  )
)$p.value


education_levels <- levels(
  babyfy_analytic$MEDLEVELA0
)


education_rows <- bind_rows(
  
  lapply(seq_along(education_levels), function(i) {
    
    lev <- education_levels[i]
    
    # Remove "(1)", "(2)", etc. from displayed labels
    clean_label <- sub(
      "^\\([0-9]+\\) ",
      "",
      lev
    )
    
    tibble(
      
      Characteristic = ifelse(
        i == 1,
        paste0("Maternal education: ", clean_label),
        paste0("    ", clean_label)
      ),
      
      Overall = n_pct(
        babyfy_analytic$MEDLEVELA0,
        lev
      ),
      
      `Low cash` = n_pct(
        babyfy_analytic$MEDLEVELA0[
          babyfy_analytic$treat_num == 0
        ],
        lev
      ),
      
      `High cash` = n_pct(
        babyfy_analytic$MEDLEVELA0[
          babyfy_analytic$treat_num == 1
        ],
        lev
      ),
      
      # One overall P value for maternal education
      `P value` = ifelse(
        i == 1,
        format_p(education_p),
        ""
      )
    )
  })
)


# ============================================================
# 6. Assemble final descriptive table
# ============================================================

table_bfy_desc <- bind_rows(
  
  # --------------------------
  # Baseline characteristics
  # --------------------------
  
  continuous_row(
    babyfy_analytic,
    "MDEPRESSIONA0",
    "Maternal depressive symptoms"
  ),
  
  continuous_row(
    babyfy_analytic,
    "HHREVISEDINCOMEA0",
    "Household income, US$",
    digits = 0
  ),
  
  education_rows,
  
  binary_row(
    babyfy_analytic,
    "biodad_num",
    "Biological father living in household"
  ),
  
  binary_row(
    babyfy_analytic,
    "female_num",
    "Female child"
  ),
  
  continuous_row(
    babyfy_analytic,
    "birthweight_g",
    "Birth weight, g",
    digits = 0
  ),
  
  
  # --------------------------
  # Child psychopathology
  # --------------------------
  
  continuous_row(
    babyfy_analytic,
    "internalizing_raw_a3",
    "Internalising symptoms, Age 3"
  ),
  
  continuous_row(
    babyfy_analytic,
    "internalizing_raw_a4",
    "Internalising symptoms, Age 4"
  ),
  
  continuous_row(
    babyfy_analytic,
    "externalizing_raw_a3",
    "Externalising symptoms, Age 3"
  ),
  
  continuous_row(
    babyfy_analytic,
    "externalizing_raw_a4",
    "Externalising symptoms, Age 4"
  )
)

### Rename table headings with Ns and %
# Get analytic sample Ns
n_overall <- n_distinct(babyfy_analytic$PUBLICSAMPLEID)

n_low <- babyfy_analytic %>%
  filter(treat_num == 0) %>%
  summarise(n = n_distinct(PUBLICSAMPLEID)) %>%
  pull(n)

n_high <- babyfy_analytic %>%
  filter(treat_num == 1) %>%
  summarise(n = n_distinct(PUBLICSAMPLEID)) %>%
  pull(n)

# Rename table headings
names(table_bfy_desc)[names(table_bfy_desc) == "Overall"] <-
  paste0("Overall (N = ", n_overall, ")")

names(table_bfy_desc)[names(table_bfy_desc) == "Low cash"] <-
  paste0("Low cash (N = ", n_low, ")")

names(table_bfy_desc)[names(table_bfy_desc) == "High cash"] <-
  paste0("High cash (N = ", n_high, ")")

table_bfy_desc

# ============================================================
# 7.Convert table to Word
# ============================================================

# Convert to formatted Word table
ft <- flextable(table_bfy_desc) %>%
  
  # Nice journal-style borders
  theme_booktabs() %>%
  
  # Bold column headings
  bold(part = "header") %>%
  
  # Left-align characteristic names
  align(
    j = "Characteristic",
    align = "left",
    part = "all"
  ) %>%
  
  # Centre numerical columns
  align(
    j = 2:5,
    align = "center",
    part = "all"
  ) %>%
  
  # Fit column widths to content
  autofit() %>%
  
  # Font and size
  font(fontname = "Arial", part = "all") %>%
  fontsize(size = 10, part = "all")

# Add table title
doc <- read_docx() %>%
  body_add_flextable(ft)

# Save Word document
print(
  doc,
  target = file.path(tables_dir, "BFY_descriptives.docx")
)

# Exploratory analyses: concurrent maternal depressive symptoms ----

# Standardise concurrent maternal depression and outcomes within age wave
babyfy_outcomes_long <- babyfy_outcomes_long %>%
  group_by(age_wave) %>%
  mutate(
    maternal_depression_z = as.numeric(scale(maternal_depression)),
    internalizing_z = as.numeric(scale(internalizing_raw)),
    externalizing_z = as.numeric(scale(externalizing_raw))
  ) %>%
  ungroup()

# Internalising
explor_concur_int <- lmer(
  internalizing_z ~ maternal_depression_z * TREAT +
    age_wave + factor(SITE) + (1 | PUBLICSAMPLEID),
  data = babyfy_outcomes_long
)

summary(explor_concur_int)

# Externalising
explor_concur_ext <- lmer(
  externalizing_z ~ maternal_depression_z * TREAT +
    age_wave + factor(SITE) + (1 | PUBLICSAMPLEID),
  data = babyfy_outcomes_long
)

summary(explor_concur_ext)

# Exploratory analyses: cash-transfer effects on maternal depression ----

explor_mat_dep <- lmer(
  maternal_depression_z ~ TREAT +
    age_wave + factor(SITE) + (1 | PUBLICSAMPLEID),
  data = babyfy_outcomes_long
)

summary(explor_mat_dep)

# Assemble BFY exploratory results panel ----

# Check coefficient names first
tidy(explor_concur_int, effects = "fixed")
tidy(explor_mat_dep, effects = "fixed")

## Extract results
# Function for formatting estimate and CI
format_est_ci <- function(est, low, high) {
  sprintf("%.3f (%.3f, %.3f)", est, low, high)
}

format_p <- function(p) {
  ifelse(p < .001, "<.001", sprintf("%.3f", p))
}

# Extract interaction from internalising model
int_result <- tidy(
  explor_concur_int,
  effects = "fixed",
  conf.int = TRUE,
  conf.method = "Wald"
) %>%
  filter(grepl("maternal_depression_z:TREAT", term))

# Extract interaction from externalising model
ext_result <- tidy(
  explor_concur_ext,
  effects = "fixed",
  conf.int = TRUE,
  conf.method = "Wald"
) %>%
  filter(grepl("maternal_depression_z:TREAT", term))

# Extract treatment effect from maternal-depression model
matdep_result <- tidy(
  explor_mat_dep,
  effects = "fixed",
  conf.int = TRUE,
  conf.method = "Wald"
) %>%
  filter(grepl("^TREAT", term))

## Construct table panel
panel_bfy <- bind_rows(
  
  tibble(
    Analysis = "Moderation of concurrent associations",
    Outcome = "Internalising symptoms",
    `Estimate (95% CI)` = format_est_ci(
      int_result$estimate,
      int_result$conf.low,
      int_result$conf.high
    ),
    `P value` = format_p(int_result$p.value)
  ),
  
  tibble(
    Analysis = "",
    Outcome = "Externalising symptoms",
    `Estimate (95% CI)` = format_est_ci(
      ext_result$estimate,
      ext_result$conf.low,
      ext_result$conf.high
    ),
    `P value` = format_p(ext_result$p.value)
  ),
  
  tibble(
    Analysis = "Subsequent maternal depressive symptoms",
    Outcome = "Maternal depressive symptoms",
    `Estimate (95% CI)` = format_est_ci(
      matdep_result$estimate,
      matdep_result$conf.low,
      matdep_result$conf.high
    ),
    `P value` = format_p(matdep_result$p.value)
  )
)

panel_bfy

# Combine GSMS and BFY exploratory results ----

# Read GSMS in
panel_gsms_exploratory <- readRDS(
  file.path(tables_dir, "panel_gsms_exploratory.rds")
)

# Add panel-heading rows
panel_A <- bind_rows(
  tibble(
    Analysis = "Panel A. Great Smoky Mountains Study",
    Outcome = "",
    `Estimate (95% CI)` = "",
    `P value` = ""
  ),
  panel_gsms_exploratory
)

panel_B <- bind_rows(
  tibble(
    Analysis = "Panel B. Baby's First Years",
    Outcome = "",
    `Estimate (95% CI)` = "",
    `P value` = ""
  ),
  panel_bfy
)

# Combine
combined_table <- bind_rows(panel_A, panel_B)

# Identify panel-heading rows
panel_rows <- which(
  combined_table$Analysis %in% c(
    "Panel A. Great Smoky Mountains Study",
    "Panel B. Baby's First Years"
  )
)

# Create flextable
ft <- flextable(combined_table)

# Column headings
ft <- set_header_labels(
  ft,
  Analysis = "Analysis",
  Outcome = "Outcome",
  `Estimate (95% CI)` = "Estimate (95% CI)",
  `P value` = "P value"
)

# Merge panel-heading rows across all columns
ft <- merge_at(
  ft,
  i = panel_rows[1],
  j = 1:4,
  part = "body"
)

ft <- merge_at(
  ft,
  i = panel_rows[2],
  j = 1:4,
  part = "body"
)

# Bold panel headings
ft <- bold(
  ft,
  i = panel_rows,
  bold = TRUE,
  part = "body"
)

# Alignment
ft <- align(
  ft,
  j = c("Analysis", "Outcome"),
  align = "left",
  part = "all"
)

ft <- align(
  ft,
  j = c("Estimate (95% CI)", "P value"),
  align = "center",
  part = "all"
)

ft <- valign(
  ft,
  valign = "center",
  part = "all"
)

ft <- autofit(ft)

ft

# Save combined exploratory results table ----
doc <- read_docx() %>%
  body_add_flextable(ft)

print(
  doc,
  target = file.path(tables_dir, "BFY_GSMS_exploratory_table.docx")
)
