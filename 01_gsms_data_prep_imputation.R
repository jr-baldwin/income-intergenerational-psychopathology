# ============================================================
# 01: Data preparation and multiple imputation
# ============================================================
#
# This script:
#
# 1. Loads the original GSMS dataset.
#
# 2. Prepares analytic, covariate, and auxiliary variables
#    for multiple imputation.
#
# 3. Performs multiple imputation of missing data.
#
# 4. Derives maternal depressive symptoms and offspring
#    psychiatric outcome variables in the imputed datasets.
#
# 5. Saves the completed imputed datasets for subsequent
#    analyses.
#
# This script requires access to the restricted GSMS data.
# The imputed datasets are saved locally and are not
# included in the public GitHub repository.
# ============================================================

# Load required packages ----

library(dplyr)
library(tidyr)
library(haven)
library(sjlabelled)
library(mice)

# Set up file paths ----

source("../config.R")

data_dir <- file.path(gsms_dir, "restricted_data")

# Location of restricted datasets
data_dir <- file.path(gsms_dir, "restricted_data")

# Load original GSMS dataset
gsms_long <- read_sas(file.path(data_dir, "Jessie.sas7bdat"))

# includes all of the historic GSMS variables including demographics, mental health, maternal depression. 
# The dataset is long, thus one observation per subject per wave. 
# It covers all obs. from age 9 to age 30. 
# The ID variable is GSMSID and the wave variable is wave. 
# The first 3 waves of the study are pre-casino and the remainder are post-casino

# Prepare variables for imputation ----

# Variable documentation ----
# See variable_dictionary.R for descriptions of the
# analytic, covariate, identifying and auxiliary variables.

#============================================================
# 1. HELPER FUNCTIONS
#============================================================

# Sum across variables, but retain NA when every component is missing
safe_rowSums <- function(x) {
  rs <- rowSums(x, na.rm = TRUE)
  all_na <- rowSums(!is.na(x)) == 0
  rs[all_na] <- NA
  rs
}

# Code:
# 1 = at least one component equals 1
# 0 = at least one component is observed, but none equals 1
# NA = all components are missing
safe_any01 <- function(x) {
  out <- as.integer(safe_rowSums(x == 1) > 0)
  all_na <- rowSums(!is.na(x)) == 0
  out[all_na] <- NA
  out
}


#============================================================
# 2. DERIVE AUXILIARY COMPOSITES IN gsms_long
#============================================================

# Pregnancy stress
pregstress_vars <- c(
  "DRISK33",
  "DRISK34"
)

gsms_long$preg_stress_any <-
  safe_any01(gsms_long[pregstress_vars])


# Pregnancy unwelcome
unwelcome_vars <- c(
  "DRISK36",
  "DRISK37"
)

gsms_long$preg_unwelcome_any <-
  safe_any01(gsms_long[unwelcome_vars])


# Parental mental health/substance service use
service_vars <- c(
  "DRISK121",
  "DRISK122",
  "DRISK123",
  "DRISK124"
)

gsms_long$parent_service_any <-
  safe_any01(gsms_long[service_vars])


# Parental arrest
arrest_vars <- c(
  "DRISK125",
  "DRISK126"
)

gsms_long$parent_arrest_any <-
  safe_any01(gsms_long[arrest_vars])


# Parenting quality/conflict:
# retain separately rather than creating one summary
parenting_vars <- c(
  "DRISK10",
  "DRISK11",
  "DRISK14",
  "DRISK24",
  "DRISK17",
  "DRISK19",
  "DRISK15"
)


# Peer/social risk
peer_vars <- c(
  "DRISK41",
  "DRISK42",
  "DRISK43",
  "DRISK44",
  "DRISK45"
)

gsms_long$peer_risk_any <-
  safe_any01(gsms_long[peer_vars])


# Recent events
violent_vars <- c(
  "DRISK50",
  "DRISK51",
  "DRISK52"
)

loss_vars <- c(
  "DRISK53",
  "DRISK77"
)

instab_vars <- c(
  "DRISK57",
  "DRISK61",
  "DRISK62",
  "DRISK65",
  "DRISK66",
  "DRISK67",
  "DRISK68",
  "DRISK72",
  "DRISK80",
  "DRISK89"
)

gsms_long$violent_events_n <-
  safe_rowSums(gsms_long[violent_vars] == 1)

gsms_long$loss_events_n <-
  safe_rowSums(gsms_long[loss_vars] == 1)

gsms_long$instability_n <-
  safe_rowSums(gsms_long[instab_vars] == 1)


#============================================================
# 3. DEFINE AUXILIARY VARIABLES FOR MI
#============================================================

aux_vars <- c(
  # Fixed or early-life variables
  "RACE3",
  "PSA0X01",
  "PSA0X02",
  "preg_stress_any",
  "preg_unwelcome_any",
  
  # Variables whose later measurements will subsequently
  # be masked within long_mi_dat
  "parent_service_any",
  "parent_arrest_any",
  parenting_vars,
  "peer_risk_any",
  "violent_events_n",
  "loss_events_n",
  "instability_n"
)


#============================================================
# 4. CREATE THE LONG MI DATASET
#============================================================

long_mi_dat <- gsms_long %>%
  dplyr::select(
    # Identifying/time variables
    GSMSID,
    WAVE,
    AGEGRUP,
    COHAGE,
    
    # Cash-transfer/race status
    RACEAI,
    AIpar,
    
    # Sampling weight
    WT2,
    
    #--------------------------------------------------------
    # Maternal depression exposure
    #--------------------------------------------------------
    DRISK22,
    SRISK22,
    SRISK22A,
    
    #--------------------------------------------------------
    # Child mental health outcomes
    #--------------------------------------------------------
    
    # Internalising
    j4nda,
    jtotanx,
    
    # Externalising
    j4ncd,
    j4nopp,
    c4nASP,
    p4nadh,
    
    # Substance use
    J4NUSE1,
    j4nsdp1,
    J4NUSE4,
    j4nsdp4,
    J4NUSE5,
    j4nsdp5,
    J4NUSE6,
    j4nsdp6,
    J4NUSE7,
    j4nsdp7,
    J4NUSE8,
    j4nsdp8,
    J4NUSE9,
    j4nsdp9,
    j5nthc,
    j5_thc,
    j5nalc,
    j5_alc,
    
    # Suicidality
    CDC9I01,
    PDC9I01,
    CDD0I01,
    PDD0I01,
    CDD1I01,
    PDD1I01,
    CDD2E01,
    PDD2E01,
    CDD2I01,
    PDD2I01,
    j4ndep9,
    
    #--------------------------------------------------------
    # Confounders
    #--------------------------------------------------------
    
    # Fixed/pre-existing confounders
    SEX,
    DRISK38,
    
    # Repeated confounder measures:
    # later values will be masked below
    DRISK4,
    POVERTY,
    INCOME,
    DRISK1,
    bmi,
    j4_core,
    
    #--------------------------------------------------------
    # Auxiliary variables
    #--------------------------------------------------------
    dplyr::all_of(aux_vars)
  )


#============================================================
# 5. RESTRICT CONFOUNDERS AND PSYCHOSOCIAL AUXILIARIES
#    TO PRE-CASH-TRANSFER WAVES 1-4
#============================================================

# These variables may vary over time and could plausibly be
# influenced by the cash-transfer programme after wave 4.
#
# Their wave 1-4 values remain available in the multilevel
# imputation model. Their wave 5+ values are masked only in
# long_mi_dat.

pre_transfer_predictors <- c(
  # Repeated confounders
  "DRISK4",
  "POVERTY",
  "INCOME",
  "DRISK1",
  "bmi",
  "j4_core",
  
  # Psychosocial auxiliaries
  parenting_vars,
  "peer_risk_any",
  "violent_events_n",
  "loss_events_n",
  "instability_n",
  "parent_service_any",
  "parent_arrest_any"
)

# Confirm all listed variables exist
setdiff(pre_transfer_predictors, names(long_mi_dat))

# This should return character(0).
# Stop with an informative error if any are absent.
missing_pre_vars <- setdiff(
  pre_transfer_predictors,
  names(long_mi_dat)
)

if (length(missing_pre_vars) > 0) {
  stop(
    "These pre-transfer predictor variables are missing from long_mi_dat: ",
    paste(missing_pre_vars, collapse = ", ")
  )
}

# Mask wave 5+ values only in the MI dataset
long_mi_dat <- long_mi_dat %>%
  mutate(
    across(
      all_of(pre_transfer_predictors),
      ~ replace(.x, WAVE > 4, NA)
    )
  )


# Confirm that no observed wave 5+ values remain
post_transfer_check <- long_mi_dat %>%
  filter(WAVE > 4) %>%
  summarise(
    across(
      all_of(pre_transfer_predictors),
      ~ sum(!is.na(.x))
    )
  )

post_transfer_check

# Every result in post_transfer_check should equal 0.


#============================================================
# 6. CARRY SAMPLING WEIGHT WITHIN PARTICIPANT
#============================================================

long_mi_dat <- long_mi_dat %>%
  group_by(GSMSID) %>%
  fill(WT2, .direction = "downup") %>%
  ungroup()

# Sampling-weight missingness check
sum(is.na(long_mi_dat$WT2))


#============================================================
# 7. PREPARE SUICIDALITY VARIABLES
#============================================================

suicide_vars_all <- c(
  "CDC9I01",
  "PDC9I01",
  "CDD0I01",
  "PDD0I01",
  "CDD1I01",
  "PDD1I01",
  "CDD2E01",
  "PDD2E01",
  "CDD2I01",
  "PDD2I01"
)

clean_suicide <- function(x) {
  x_char <- as.character(x)
  x_char[x_char %in% c("", "<NA>")] <- NA
  as.numeric(x_char)
}

long_mi_dat <- long_mi_dat %>%
  mutate(
    across(
      all_of(suicide_vars_all),
      clean_suicide
    )
  )

# Check classes and missingness
sapply(
  long_mi_dat[suicide_vars_all],
  function(x) {
    c(
      class = class(x)[1],
      n_na = sum(is.na(x))
    )
  }
)


#============================================================
# 8. GENERAL DIAGNOSTICS BEFORE MICE
#============================================================

# Variables that are completely missing
all_missing <-
  colSums(is.na(long_mi_dat)) == nrow(long_mi_dat)

names(all_missing[all_missing])


# Variables with no observed variation
no_var <- sapply(
  long_mi_dat,
  function(x) {
    length(unique(na.omit(x))) <= 1
  }
)

names(no_var[no_var])


# Variation by wave
wave_var_check <- long_mi_dat %>%
  group_by(WAVE) %>%
  summarise(
    across(
      everything(),
      ~ n_distinct(na.omit(.x))
    ),
    .groups = "drop"
  )

constant_within_wave <- names(wave_var_check)[
  sapply(
    wave_var_check,
    function(x) all(x <= 1)
  )
]

constant_within_wave


# Raw overall missingness in the MI dataset,
# including deliberately masked post-transfer cells
miss_prop_all <- sort(
  colMeans(is.na(long_mi_dat)),
  decreasing = TRUE
)
miss_prop_all

# Relevant observed missingness during waves 1-4
miss_prop_pre <- long_mi_dat %>%
  filter(WAVE <= 4) %>%
  summarise(
    across(
      all_of(pre_transfer_predictors),
      ~ mean(is.na(.x))
    )
  )
miss_prop_pre
# The apparent overall missingness of about 59% was therefore overwhelmingly caused 
# by deliberately masking waves 5 onward. It was not evidence that nearly 60% of 
# the pre-transfer income information was genuinely absent.

#============================================================
# 9. AUXILIARY-VARIABLE CORRELATIONS
#============================================================

aux_check <- aux_vars[
  aux_vars %in% names(long_mi_dat)
]

aux_numeric <- aux_check[
  sapply(long_mi_dat[aux_check], is.numeric)
]

cors_child_dep <- sapply(
  aux_numeric,
  function(v) {
    suppressWarnings(
      cor(
        long_mi_dat[[v]],
        long_mi_dat$j4nda,
        use = "pairwise.complete.obs"
      )
    )
  }
)

cors_exp <- sapply(
  aux_numeric,
  function(v) {
    suppressWarnings(
      cor(
        long_mi_dat[[v]],
        long_mi_dat$SRISK22,
        use = "pairwise.complete.obs"
      )
    )
  }
)

aux_corrs <- data.frame(
  variable = aux_numeric,
  cor_with_outcome = cors_child_dep,
  cor_with_exposure = cors_exp,
  abs_cor_with_outcome = abs(cors_child_dep),
  abs_cor_with_exposure = abs(cors_exp)
)

aux_corrs_outcome <- aux_corrs[
  order(-aux_corrs$abs_cor_with_outcome),
]

aux_corrs_exposure <- aux_corrs[
  order(-aux_corrs$abs_cor_with_exposure),
]

aux_corrs_outcome
aux_corrs_exposure

ncol(long_mi_dat)


#============================================================
# 10. PARENT-INFORMANT CHECKS
#============================================================

parent_vars <- c(
  "PDC9I01",
  "PDD0I01",
  "PDD1I01",
  "PDD2E01",
  "PDD2I01"
)

parent_vars <- intersect(
  parent_vars,
  names(long_mi_dat)
)

parent_wave_summary <- long_mi_dat %>%
  group_by(WAVE) %>%
  summarise(
    parent_any = any(
      rowSums(
        !is.na(
          across(all_of(parent_vars))
        )
      ) > 0
    ),
    .groups = "drop"
  )

parent_wave_summary


parent_missing_by_wave <- long_mi_dat %>%
  group_by(WAVE) %>%
  summarise(
    across(
      all_of(parent_vars),
      ~ mean(is.na(.x))
    ),
    .groups = "drop"
  )

parent_missing_by_wave
# Parent reports of suicidality are completely missing from wave 9 onwards
# But partially missing from wave 7, as cohort 3 already turned 19 and did not 
# include parent reports at this age

#============================================================
# 11. CONVERT SUICIDALITY VARIABLES TO FACTORS
#============================================================

suicide_vars_present <- intersect(
  suicide_vars_all,
  names(long_mi_dat)
)

long_mi_dat <- long_mi_dat %>%
  mutate(
    across(
      all_of(suicide_vars_present),
      factor
    )
  )


#============================================================
# 12. INITIAL MICE SETUP
#============================================================

ini <- mice(
  long_mi_dat,
  maxit = 0,
  printFlag = FALSE
)

meth <- ini$method
pred <- ini$predictorMatrix
where <- is.na(long_mi_dat)


#============================================================
# 13. IDENTIFIERS AND TIME VARIABLES
#============================================================

# Participant ID
meth["GSMSID"] <- ""

pred[, "GSMSID"] <- 0
pred["GSMSID", ] <- 0


# Wave and cohort variables
meth["WAVE"] <- ""
meth["AGEGRUP"] <- ""

pred[, "WAVE"] <- 0
pred["WAVE", ] <- 0

pred[, "AGEGRUP"] <- 0
pred["AGEGRUP", ] <- 0


#============================================================
# 14. REDUNDANT EXPOSURE VARIABLE
#============================================================

# SRISK22 is excluded because SRISK22A is retained
meth["SRISK22"] <- ""

pred[, "SRISK22"] <- 0
pred["SRISK22", ] <- 0


#============================================================
# 15. PREDICTOR-ONLY VARIABLES
#============================================================

# AI parental-status variable
meth["AIpar"] <- ""

pred["AIpar", ] <- 0
pred[, "AIpar"] <- 1
pred["AIpar", "AIpar"] <- 0


# Sampling weight
meth["WT2"] <- ""

pred["WT2", ] <- 0
pred[, "WT2"] <- 1
pred["WT2", "WT2"] <- 0


#============================================================
# 16. STRUCTURAL MISSINGNESS
#============================================================

# Parent-report variables were not assessed from wave 7 (cohort 3), 9 (cohort 2), 11 (cohort 1) onward
# This corresponds to when participants turned 19 
where[
  (
    (long_mi_dat$AGEGRUP == 3 & long_mi_dat$WAVE >= 7) |
      (long_mi_dat$AGEGRUP == 2 & long_mi_dat$WAVE >= 9) |
      (long_mi_dat$AGEGRUP == 1 & long_mi_dat$WAVE >= 11)
  ),
  parent_vars
] <- FALSE

# Important:
# Do not impute the deliberately masked post-transfer values
# of the confounders and auxiliary variables.
#
# Otherwise MICE would refill the wave 5+ values that were
# deliberately removed above.

where[
  long_mi_dat$WAVE > 4,
  pre_transfer_predictors
] <- FALSE


#============================================================
# 17. SUICIDALITY IMPUTATION METHODS
#============================================================

meth[suicide_vars_present] <- "polyreg"


# Sparse heroin variables
meth["J4NUSE7"] <- "cart"
meth["j4nsdp7"] <- "cart"


# Do not use suicidality variables as predictors
pred[, suicide_vars_present] <- 0
pred[suicide_vars_present, suicide_vars_present] <- 0


#============================================================
# 18. RE-ENFORCE PREDICTOR-ONLY SETTINGS
#============================================================

pred[, "WT2"] <- 1
pred["WT2", ] <- 0
pred["WT2", "WT2"] <- 0

pred[, "AIpar"] <- 1
pred["AIpar", ] <- 0
pred["AIpar", "AIpar"] <- 0


#============================================================
# 19. REMOVE SPARSE HEROIN VARIABLES AS PREDICTORS
#============================================================

pred[, "J4NUSE7"] <- 0
pred[, "j4nsdp7"] <- 0

meth[c("J4NUSE7", "j4nsdp7")]

pred[
  1:5,
  c("J4NUSE7", "j4nsdp7")
]


#============================================================
# 20. CHECK THE WAVE RESTRICTIONS IN `where`
#============================================================

# These should all be FALSE: MICE should not impute the
# deliberately masked post-transfer values.
table(
  where[
    long_mi_dat$WAVE > 4,
    "INCOME"
  ]
)

sapply(
  pre_transfer_predictors,
  function(v) {
    any(
      where[
        long_mi_dat$WAVE > 4,
        v
      ]
    )
  }
)

# Every value returned above should be FALSE.

# Do not impute ADHD after childhood due to structural missingness
where[
  (
    (long_mi_dat$AGEGRUP == 3 & long_mi_dat$WAVE > 4) |
      (long_mi_dat$AGEGRUP == 2 & long_mi_dat$WAVE > 6) |
      (long_mi_dat$AGEGRUP == 1 & long_mi_dat$WAVE > 8)
  ),
  "p4nadh"
] <- FALSE

# See below for justification (missingness by cohort and wave)
long_mi_dat %>%
  group_by(AGEGRUP, WAVE) %>%
  summarise(
    n_observed = sum(!is.na(p4nadh)),
    .groups = "drop"
  ) %>%
  print(n = Inf)

## Do not impute c4nASP (ASPD) in adolescence due to structural missingness
# c4nASP first assessed at different waves by cohort
where[
  (
    (long_mi_dat$AGEGRUP == 3 & long_mi_dat$WAVE < 7) |
      (long_mi_dat$AGEGRUP == 2 & long_mi_dat$WAVE < 9) |
      (long_mi_dat$AGEGRUP == 1 & long_mi_dat$WAVE < 11)
  ),
  "c4nASP"
] <- FALSE

## Remove pre-transfer-only variables as predictors of c4nASP (assessed only in adolescence):
# Start with no predictors for c4nASP
pred["c4nASP", ] <- 0

# Candidate predictors likely to be available at adolescent/adult waves
c4_candidate_preds <- c(
  "RACEAI",
  "AIpar",
  "WT2",
  "SEX",
  "COHAGE",
  
  # Other child MH measures
  "j4nda",
  "jtotanx",
  "j4ncd",
  "j4nopp",
  
  # Substance-use measures
  "j5nthc",
  "j5_thc",
  "j5nalc",
  "j5_alc"
)

c4_candidate_preds <- intersect(
  c4_candidate_preds,
  colnames(pred)
)

pred["c4nASP", c4_candidate_preds] <- 1

# Do not use the target to predict itself
pred["c4nASP", "c4nASP"] <- 0

c4_preds <- names(which(pred["c4nASP", ] != 0))

sum(
  !is.na(long_mi_dat$c4nASP) &
    complete.cases(
      long_mi_dat[, c4_preds, drop = FALSE]
    )
)

data.frame(
  predictor = c4_preds,
  n_jointly_observed = sapply(
    c4_preds,
    function(v) {
      sum(
        !is.na(long_mi_dat$c4nASP) &
          !is.na(long_mi_dat[[v]])
      )
    }
  )
) %>%
  arrange(n_jointly_observed)

# Remove heroin vars as predictors
heroin_vars <- c(
  "J4NUSE4","J4NUSE5","J4NUSE6","J4NUSE7","J4NUSE8","J4NUSE9",
  "j4nsdp4","j4nsdp5","j4nsdp6","j4nsdp7","j4nsdp8","j4nsdp9"
)

pred[, heroin_vars] <- 0

# Pre-transfer confounders and auxiliary variables were imputed
# using a reduced predictor set restricted to variables available
# during the pre-transfer period. This avoids age-specific variables
# (e.g. ASPD) preventing imputations because of structural missingness.

# Variables that are reliable predictors during waves 1-4
core_pre_predictors <- c(
  "WAVE",
  "COHAGE",
  "RACEAI",
  "AIpar",
  "WT2",
  "SEX",
  "DRISK38"
)

core_pre_predictors <- intersect(
  core_pre_predictors,
  colnames(pred)
)

# Pre-transfer variables can also predict one another
pre_vars_available <- intersect(
  pre_transfer_predictors,
  colnames(pred)
)

for (target in pre_vars_available) {
  
  # Skip variables that are not imputed
  if (meth[target] == "") next
  
  # Reset this target's predictor row
  pred[target, ] <- 0
  
  # Add stable demographic/time predictors
  pred[target, core_pre_predictors] <- 1
  
  # Allow the other pre-transfer variables to predict it
  other_pre_vars <- setdiff(pre_vars_available, target)
  pred[target, other_pre_vars] <- 1
  
  # Never predict a variable from itself
  pred[target, target] <- 0
}

meth["WAVE"] <- ""
pred["WAVE", ] <- 0

# ------------------------------------------------------------
# Predictor model for maternal depressive symptoms
# ------------------------------------------------------------
eligible_dep_rows <-
  long_mi_dat$WAVE <= 4 &
  is.na(long_mi_dat$SRISK22A) &
  where[, "SRISK22A"]

dep_candidates <- c(
  "WAVE",
  "COHAGE",
  "RACEAI",
  "AIpar",
  "WT2",
  "SEX",
  "DRISK38",
  "DRISK22",
  "j4nda",
  "jtotanx",
  pre_transfer_predictors
)

dep_candidates <- intersect(
  dep_candidates,
  names(long_mi_dat)
)

dep_availability <- data.frame(
  predictor = dep_candidates,
  n_missing_in_eligible_rows = vapply(
    dep_candidates,
    function(v) {
      sum(is.na(long_mi_dat[[v]][eligible_dep_rows]))
    },
    integer(1)
  )
) |>
  arrange(n_missing_in_eligible_rows)

dep_availability

srisk22a_predictors <- dep_availability |>
  filter(n_missing_in_eligible_rows == 0) |>
  pull(predictor)

pred["SRISK22A", ] <- 0
pred["SRISK22A", srisk22a_predictors] <- 1
pred["SRISK22A", "SRISK22A"] <- 0

srisk22a_predictors

# ------------------------------------------------------------
# Predictor model for low birth weight
# ------------------------------------------------------------

# DRISK38 is a fixed early-life characteristic.
# Use a reduced set of stable, complete predictors to avoid
# structurally missing age-specific variables preventing imputation.

pred["DRISK38", ] <- 0

bw_predictors <- c(
  "COHAGE",
  "RACEAI",
  "AIpar",
  "WT2",
  "SEX",
  "RACE3"
)

bw_predictors <- intersect(
  bw_predictors,
  colnames(pred)
)

pred["DRISK38", bw_predictors] <- 1
pred["DRISK38", "DRISK38"] <- 0

# ------------------------------------------------------------
# Predictor model for ADHD
# ------------------------------------------------------------

pred["p4nadh", ] <- 0

adhd_predictors <- c(
  # Stable participant characteristics
  "COHAGE",
  "RACEAI",
  "AIpar",
  "WT2",
  "SEX",
  "DRISK38",
  "RACE3",
  
  # Concurrent child psychiatric measures
  "j4nda",
  "jtotanx",
  "j4ncd",
  "j4nopp",
  
  # Substance indicators that are available
  "J4NUSE1",
  "j5nthc",
  "j5_thc",
  "j5nalc",
  "j5_alc"
)

adhd_predictors <- intersect(
  adhd_predictors,
  colnames(pred)
)

pred["p4nadh", adhd_predictors] <- 1
pred["p4nadh", "p4nadh"] <- 0

# ------------------------------------------------------------
# Predictor model for nicotine dependence
# ------------------------------------------------------------

pred["j4nsdp1", ] <- 0

nicdep_predictors <- c(
  "COHAGE",
  "RACEAI",
  "AIpar",
  "WT2",
  "SEX",
  "DRISK38",
  "RACE3",
  
  "J4NUSE1",
  
  "j4nda",
  "jtotanx",
  "j4ncd",
  
  "j5nthc",
  "j5_thc",
  "j5nalc",
  "j5_alc"
)

pred["j4nsdp1", ] <- 0
pred["j4nsdp1", nicdep_predictors] <- 1
pred["j4nsdp1", "j4nsdp1"] <- 0

# ------------------------------------------------------------
# Predictor model for nicotine use (J4NUSE1)
# ------------------------------------------------------------

pred["J4NUSE1", ] <- 0

use_predictors <- c(
  # Stable participant characteristics
  "COHAGE",
  "RACEAI",
  "AIpar",
  "WT2",
  "SEX",
  "DRISK38",
  "RACE3",
  
  # Concurrent psychiatric measures that are available
  "j4nda",
  "jtotanx",
  
  # Other substance-use measures that are complete on the failed rows
  "j5nthc",
  "j5_thc",
  "j5nalc",
  "j5_alc"
)

use_predictors <- intersect(
  use_predictors,
  colnames(pred)
)

pred["J4NUSE1", use_predictors] <- 1
pred["J4NUSE1", "J4NUSE1"] <- 0

# ------------------------------------------------------------
# Predictor model for conduct disorder symptoms (j4ncd)
# ------------------------------------------------------------

pred["j4ncd", ] <- 0

cd_predictors <- c(
  # Stable participant characteristics
  "COHAGE",
  "RACEAI",
  "AIpar",
  "WT2",
  "SEX",
  "DRISK38",
  "RACE3",
  
  # Psychiatric measures available on the failed rows
  "j4nda",
  "jtotanx",
  
  # Substance-use measures available on the failed rows
  "J4NUSE1",
  "j5nthc",
  "j5_thc",
  "j5nalc",
  "j5_alc"
)

cd_predictors <- intersect(
  cd_predictors,
  colnames(pred)
)

pred["j4ncd", cd_predictors] <- 1
pred["j4ncd", "j4ncd"] <- 0

# ------------------------------------------------------------
# Wave 17 assessment pattern
# ------------------------------------------------------------

wave17_structural <- 
  long_mi_dat$WAVE == 17 &
  is.na(long_mi_dat$j4ncd) &
  is.na(long_mi_dat$j4nopp) &
  is.na(long_mi_dat$c4nASP) &
  is.na(long_mi_dat$p4nadh) &
  is.na(long_mi_dat$J4NUSE1) &
  is.na(long_mi_dat$j4nsdp1) &
  !is.na(long_mi_dat$j4nda) &
  !is.na(long_mi_dat$jtotanx)

sum(wave17_structural)
# should be 11

wave17_structural_vars <- c(
  "j4ncd",
  "j4nopp",
  "c4nASP",
  "p4nadh",
  "J4NUSE1",
  "j4nsdp1"
)

where[
  wave17_structural,
  wave17_structural_vars
] <- FALSE

# ------------------------------------------------------------
# Predictor model for ODD symptoms (j4nopp)
# ------------------------------------------------------------

pred["j4nopp", ] <- 0

opp_predictors <- c(
  "COHAGE",
  "RACEAI",
  "WT2",
  "SEX",
  "DRISK38",
  "RACE3",
  
  "j4nda",
  "jtotanx",
  "j4ncd",
  "c4nASP",
  
  "J4NUSE1",
  "j4nsdp1",
  "j5nthc",
  "j5_thc",
  "j5nalc",
  "j5_alc"
)

opp_predictors <- intersect(
  opp_predictors,
  colnames(pred)
)

pred["j4nopp", opp_predictors] <- 1
pred["j4nopp", "j4nopp"] <- 0

# ------------------------------------------------------------
# Predictor model for binary maternal depression (DRISK22)
# ------------------------------------------------------------

pred["DRISK22", ] <- 0

drisk22_predictors <- c(
  # Stable participant characteristics
  "COHAGE",
  "RACEAI",
  "AIpar",
  "WT2",
  "SEX",
  "DRISK38",
  "RACE3",
  
  # Maternal depressive symptoms
  "SRISK22A",
  
  # Psychiatric measures
  "j4nda",
  "jtotanx",
  "j4ncd",
  "j4nopp",
  "p4nadh",
  "J4NUSE1",
  "j4nsdp1",
  "j5nthc",
  "j5_thc",
  "j5nalc",
  "j5_alc",
  "j4ndep9",
  
  # Pre-transfer family/individual characteristics
  "DRISK4",
  "POVERTY",
  "INCOME",
  "DRISK1",
  "bmi",
  "j4_core",
  
  # Auxiliaries available on these rows
  "parent_service_any",
  "parent_arrest_any",
  "DRISK10",
  "DRISK11",
  "DRISK14",
  "DRISK24",
  "DRISK17",
  "DRISK19",
  "DRISK15",
  "peer_risk_any",
  "violent_events_n",
  "loss_events_n",
  "instability_n"
)

drisk22_predictors <- intersect(
  drisk22_predictors,
  colnames(pred)
)

pred["DRISK22", drisk22_predictors] <- 1
pred["DRISK22", "DRISK22"] <- 0

# ------------------------------------------------------------
# Reduced predictor models for child suicidality components
# ------------------------------------------------------------

# Treat anomalous child attempt code as missing BEFORE imputation
long_mi_dat$CDD2E01[long_mi_dat$CDD2E01 == 1] <- NA

# Child-reported suicidality components
child_suic_vars <- c(
  "CDC9I01",
  "CDD0I01",
  "CDD1I01",
  "CDD2E01",
  "CDD2I01"
)

# Identify rows where the entire child suicidality block is absent
child_suic_block_missing <- rowSums(
  !is.na(long_mi_dat[, child_suic_vars])
) == 0

sum(child_suic_block_missing)
# Expected: 276

# Do not impute these structurally absent cells
where[
  child_suic_block_missing,
  child_suic_vars
] <- FALSE



# Give the four problematic child items reduced predictor models
suic_targets <- c(
  "CDC9I01",
  "CDD0I01",
  "CDD2E01",
  "CDD2I01"
)

suic_predictors <- c(
  # Stable characteristics
  "COHAGE",
  "RACEAI",
  "AIpar",
  "WT2",
  "SEX",
  "DRISK38",
  "RACE3",
  
  # Psychiatric symptoms
  "j4nda",
  "jtotanx",
  "j4ncd",
  "j4nopp",
  
  # Substance measures
  "J4NUSE1",
  "j4nsdp1",
  "j5nthc",
  "j5_thc",
  "j5nalc",
  "j5_alc"
)

suic_predictors <- intersect(
  suic_predictors,
  colnames(pred)
)

for (v in suic_targets) {
  pred[v, ] <- 0
  pred[v, suic_predictors] <- 1
  
  # A variable cannot predict itself
  pred[v, v] <- 0
}


# Multiple imputation ----

imp_test <- mice(
  data = long_mi_dat,
  m = 1,
  maxit = 5,
  method = meth,
  predictorMatrix = pred,
  where = where,
  seed = 20260316,
  printFlag = TRUE
)

test_complete <- complete(imp_test, 1)

saveRDS(
  list(
    meth = meth,
    pred = pred,
    where = where
  ),
  file.path(data_dir, "mi_settings_20260818.rds")
)

saveRDS(
  long_mi_dat,
  file.path(data_dir, "long_mi_dat_20260818.rds")
)

## Run imputation
imp_0826 <- mice(
  data = long_mi_dat,
  m = 20,
  maxit = 20,
  method = meth,
  predictorMatrix = pred,
  where = where,
  seed = 20260316,
  printFlag = TRUE
)

imp_0826$loggedEvents

## Check variables were imputed correctly:
imp_list_check <- complete(imp_0826, action = "all")

# Save mids object
saveRDS(
  list(
    meth = meth,
    pred = pred,
    where = where
  ),
  file.path(data_dir, "mi_settings_20260818.rds")
)

# Check residual missingness in key variables targeted for imputation
# Structural missingness is excluded using the original where matrix

final_check_vars <- c(
  "DRISK38",
  "DRISK4",
  "POVERTY",
  "INCOME",
  "DRISK1",
  "bmi",
  "j4_core",
  "j4nda",
  "jtotanx",
  "j4ncd",
  "j4nopp",
  "c4nASP",
  "p4nadh",
  "J4NUSE1",
  "j4nsdp1",
  "CDC9I01",
  "CDD0I01",
  "CDD1I01",
  "CDD2E01",
  "CDD2I01"
)

# Count eligible cells that remain missing in each imputation
final_audit <- purrr::map_dfr(
  seq_len(imp_0826$m),
  function(i) {
    
    dat_i <- complete(imp_0826, i)
    
    data.frame(
      imputation = i,
      variable = final_check_vars,
      still_missing = sapply(
        final_check_vars,
        function(v) {
          sum(
            where[, v] &
              is.na(dat_i[[v]])
          )
        }
      )
    )
  }
)

# Summarise across all 20 imputations
final_audit_summary <- final_audit |>
  dplyr::group_by(variable) |>
  dplyr::summarise(
    min_missing = min(still_missing),
    max_missing = max(still_missing),
    .groups = "drop"
  )

final_audit_summary

# All 20 imputations had zero residual eligible missingness
# for the variables included in final_check_vars.

##### Safe to proceed to the analysis stage

# Derive variables in imputed datasets ----

imp_list <- complete(imp_0826, action = "all")

#### Derive maternal depression at each wave ####
imp_list <- lapply(imp_list, function(d) {
  
  # Maternal depression at each wave:
  # use biological mother if present, otherwise resident mother
  d <- d %>%
    mutate(
      mat_dep = dplyr::coalesce(SRISK22, SRISK22A)
    )
  
  d
})

# Check the derived variable in the first imputed dataset
summary(imp_list[[1]]$mat_dep)

# Check for any missingness in this variable 
sapply(imp_list, function(d) sum(is.na(d$mat_dep)))

## Derive internalising problems at each wave ####
imp_list <- lapply(imp_list, function(d) {
  
  d$int <- NA_real_
  
  cells <- unique(d[, c("WAVE", "AGEGRUP")])
  
  for (i in seq_len(nrow(cells))) {
    
    w <- cells$WAVE[i]
    a <- cells$AGEGRUP[i]
    
    idx <- d$WAVE == w & d$AGEGRUP == a
    
    Z <- scale(d[idx, c("j4nda", "jtotanx")]) # standardise by wave and cohort
    
    # Mean composite of depression and anxiety
    d$int[idx] <- rowMeans(Z, na.rm = TRUE)
  }
  
  d$int[is.nan(d$int)] <- NA_real_
  
  d
})

# Check the distribution of the derived composite
summary(imp_list[[1]]$int)

# Check residual missingness across imputations
sapply(imp_list, function(d) sum(is.na(d$int)))

#### Derive substance abuse factor ####
## derive the substance severity variables inside each imputed dataset

# Create a make severity function
make_severity <- function(use, dep) {
  ifelse(use == 1 & dep == 1, 2,
         ifelse((use == 1 & dep == 0) | (use == 0 & dep == 1), 1,
                ifelse(use == 0 & dep == 0, 0, NA)))
}

# Apply it to each imputed dataset
imp_list <- lapply(imp_list, function(d) {
  
  ## 0/1/2 severity for substances with use + dependence indicators
  d$sev_nicotine <- make_severity(d$J4NUSE1, d$j4nsdp1)
  d$sev_cocaine  <- make_severity(d$J4NUSE4, d$j4nsdp4)
  d$sev_amphet   <- make_severity(d$J4NUSE5, d$j4nsdp5)
  d$sev_inhal    <- make_severity(d$J4NUSE6, d$j4nsdp6)
  d$sev_heroin   <- make_severity(d$J4NUSE7, d$j4nsdp7)
  d$sev_halluc   <- make_severity(d$J4NUSE8, d$j4nsdp8)
  d$sev_sed      <- make_severity(d$J4NUSE9, d$j4nsdp9)
  
  ## Cannabis severity
  d$sev_cannabis <- NA_real_
  d$sev_cannabis[d$j5nthc == 0] <- 0
  d$sev_cannabis[d$j5nthc %in% 1:2] <- 1
  d$sev_cannabis[d$j5nthc >= 3] <- 2
  
  ## Alcohol severity
  d$sev_alcohol <- NA_real_
  d$sev_alcohol[d$j5nalc == 0] <- 0
  d$sev_alcohol[d$j5nalc %in% 1:2] <- 1
  d$sev_alcohol[d$j5nalc >= 3] <- 2
  
  d
})

# Fit a one-factor CFA to nine ordered substance-use severity indicators
# separately in each imputed dataset (WLSMV estimator).
# Extract participant-wave factor scores for the externalising composite

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

ordered_vars <- c(
  "sev_alcohol",
  "sev_cannabis",
  "sev_nicotine",
  "sev_cocaine",
  "sev_amphet",
  "sev_inhal",
  "sev_heroin",
  "sev_halluc",
  "sev_sed"
)

imp_list <- lapply(imp_list, function(d) {
  
  ## make ordered categorical
  d[ordered_vars] <- lapply(d[ordered_vars], function(x) {
    ordered(x, levels = c(0, 1, 2))
  })
  
  ## fit CFA
  substance_fit <- cfa(
    substance_model,
    data = d,
    ordered = ordered_vars,
    estimator = "WLSMV"
  )
  
  ## extract factor scores
  fs  <- lavPredict(substance_fit)[, 1]
  idx <- lavInspect(substance_fit, "case.idx")
  
  d$substance_score <- NA_real_
  d$substance_score[idx] <- fs
  
  d
})

# Inspect the distribution of CFA factor scores
summary(imp_list[[1]]$substance_score)

#### Derive externalising problems composite #####
imp_list <- lapply(imp_list, function(d) {
  
  ext_vars <- c(
    "j4ncd",          # conduct problems
    "j4nopp",         # ODD
    "c4nASP",         # ASPD
    "p4nadh",         # ADHD (continuous)
    "substance_score" # CFA factor score; standardised below within wave × age group
  )
  
  d$ext_comp <- NA_real_
  
  cells <- unique(d[, c("WAVE", "AGEGRUP")])
  
  for (i in seq_len(nrow(cells))) {
    
    w <- cells$WAVE[i]
    a <- cells$AGEGRUP[i]
    
    idx <- d$WAVE == w & d$AGEGRUP == a
    X <- d[idx, ext_vars, drop = FALSE]
    
    # keep indicators that vary within this wave x cohort cell
    keep <- sapply(X, function(x) {
      sdx <- sd(x, na.rm = TRUE)
      !is.na(sdx) && sdx > 0
    })
    
    # Need at least 2 usable indicators in the cell
    if (sum(keep) < 2) next
    
    # Standardise each usable indicator within the cell
    Z <- scale(X[, keep, drop = FALSE])
    
    # Number of observed indicators for each participant-wave
    n_observed <- rowSums(!is.na(Z))
    
    # Require at least two observed indicators per participant-wave
    d$ext_comp[idx] <- ifelse(
      n_observed >= 2,
      rowMeans(Z, na.rm = TRUE),
      NA_real_
    )
  }
  
  d$ext_comp[is.nan(d$ext_comp)] <- NA_real_
  
  d
})

# Inspect the distribution of externalising scores
summary(imp_list[[1]]$ext_comp)

# Check missingness across imputations
sapply(imp_list, function(d) sum(is.na(d$ext_comp)))

#### Derive suicidality variable #####

# First derive a narrower indicator of suicidal thoughts, plans and attempts.
# This is an intermediate variable used to construct the primary
# suicidality outcome, which additionally includes passive death wishes. 

## Make any suicidality var = suicidal ideation OR suicide plan OR suicide attempt
make_suicidality <- function(dat) {
  
  # Treat anomalous child attempt code as missing
  dat$CDD2E01[dat$CDD2E01 == 1] <- NA
  
  # Child report
  child_thoughts <- dat$CDD0I01 %in% c(2, 3) # suicidal thoughts (2=at least sometimes, 3=usually)
  child_plans    <- dat$CDD1I01 %in% c(2, 3) # suicidal plans, (2=a specific plan considered >1 occasion, 3 plan with action taken)
  child_attempt1 <- dat$CDD2E01 %in% 2 # Suicidal behavior (2=present)
  child_attempt2 <- dat$CDD2I01 %in% 2 # Suicidal behavior in past 3 months (2=present)
  
  # Parent report
  parent_thoughts <- dat$PDD0I01 %in% c(2, 3) # suicidal thoughts (2=at least sometimes, 3=usually)
  parent_plans    <- dat$PDD1I01 %in% c(2, 3) # suicidal plans, (2=a specific plan considered >1 occasion, 3 plan with action taken)
  parent_attempt1 <- dat$PDD2E01 %in% 2 # Suicidal behavior (2=present)
  parent_attempt2 <- dat$PDD2I01 %in% 2 # Suicidal behavior in past 3 months (2=present)
  
  # Whether any relevant item was observed
  child_observed <- rowSums(!is.na(dat[, c(
    "CDD0I01", "CDD1I01", "CDD2E01", "CDD2I01"
  )])) > 0
  
  parent_observed <- rowSums(!is.na(dat[, c(
    "PDD0I01", "PDD1I01", "PDD2E01", "PDD2I01"
  )])) > 0
  
  # Code child positive if they reported any evidence of thoughts, plans, or suicide attempts
  child_positive <-
    child_thoughts |
    child_plans |
    child_attempt1 |
    child_attempt2
  
  # Code parent positive if they reported any evidence of thoughts, plans, or suicide attempts
  parent_positive <-
    parent_thoughts |
    parent_plans |
    parent_attempt1 |
    parent_attempt2
  
  dat$any_suicidality <- ifelse(
    dat$WAVE < 9,
    
    # Waves 1–8: child or parent report
    ifelse(
      child_positive | parent_positive,
      1L,
      ifelse(
        child_observed | parent_observed,
        0L,
        NA_integer_
      )
    ),
    
    # Wave 9 onwards: child report only
    ifelse(
      child_positive,
      1L,
      ifelse(
        child_observed,
        0L,
        NA_integer_
      )
    )
  )
  
  dat
}

imp_list <- lapply(imp_list, make_suicidality)
prop.table(table(imp_list[[1]]$any_suicidality, useNA = "ifany")) #3%

##------------------------------------------------------------
## Derive broad binary suicidality outcome
##------------------------------------------------------------
## Outcome includes:
##   - wanting to die (passive death wishes)
##   - suicidal thoughts
##   - suicidal plans
##   - suicide attempts/behaviour
##
## Informants:
##   - Waves 1–8: child OR parent report
##   - Waves 9 onwards: child report only (parent interview no longer available)
##
## A participant is coded as having suicidality (1) if any of the above
## symptoms are endorsed by the relevant informant(s) at that assessment.

add_broad_suicidality <- function(dat) {
  
  # Wanting-to-die indicators
  ## Codes 2 or 3 indicate endorsement
  child_death_wish  <- dat$CDC9I01 %in% c(2, 3)
  parent_death_wish <- dat$PDC9I01 %in% c(2, 3)
  
  # Whether wanting-to-die information was observed
  child_death_wish_observed  <- !is.na(dat$CDC9I01)
  parent_death_wish_observed <- !is.na(dat$PDC9I01)
  
  # Positive indicator for the narrower measure.
  # %in% returns FALSE, rather than NA, when any_suicidality is missing.
  narrower_positive <- dat$any_suicidality %in% 1
  
  dat$any_suicidality_broad <- ifelse(
    dat$WAVE < 9,
    
    # Waves 1–8: combine child and parent reports
    ifelse(
      !is.na(dat$any_suicidality) | # if at least 1 source of suicidality info was observed
        child_death_wish_observed |
        parent_death_wish_observed,
      
      as.integer(
        narrower_positive | # Code 1 if thoughts/plans/attempts were present, or
          child_death_wish | # child reported wanting to die, or
          parent_death_wish # parent reported that the child wanted to die
      ),
      
      # All relevant information was missing
      NA_integer_
    ),
    
    # Wave 9 onwards: child report only
    ifelse(
      !is.na(dat$any_suicidality) | # if at least one child self-report was observed
        child_death_wish_observed,
      
      # Code 1 if thoughts/plans/attempts or wanting to die were present
      as.integer(
        narrower_positive |
          child_death_wish
      ),
      
      # All relevant child-report information was missing
      NA_integer_
    )
  )
  
  dat
}

## Apply to each imputed dataset
imp_list <- lapply(imp_list, add_broad_suicidality)

# Check prevalence and missingness of the broad suicidality outcome
prop.table(table(
  imp_list[[1]]$any_suicidality_broad,
  useNA = "ifany"
))

# Save final imputed datasets ----
saveRDS(
  imp_list,
  file.path(data_dir, "imp_list_derivedVars_20260818.rds")
)