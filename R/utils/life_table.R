# Age- and sex-specific background mortality (analysis_plan.md §7.1 item 7: "US life tables,
# age- and sex-specific; no CD excess mortality (per Aliyev / TA456)") -- the sourced input a
# lifetime horizon needs and that R/05_deterministic_results.R's own module header flagged as
# missing ("no life-table data has been sourced anywhere in this repository. A true lifetime run
# is a real, flagged gap").
#
# LIFETIME HORIZON, 2026-08-05: this file plus R/utils/transition_matrix.R's age_adjust_matrix()
# and the death_prob_schedule wiring through R/02/R/03/R/05 are what close that gap. See
# R/05_deterministic_results.R's module header for the end-to-end picture of what changed.
#
# ---- Source -----------------------------------------------------------------------------------
#
# NCHS "United States Life Tables, 2021" (National Vital Statistics Reports, Vol. 72, No. 12,
# Nov 7 2023), Table 2 (life table for males, pp. 16-17) and Table 3 (life table for females,
# pp. 18-19) -- data/raw/nchs_us_life_tables_2021.csv, transcribed directly from the published
# PDF's table images (https://www.cdc.gov/nchs/data/nvsr/nvsr72/nvsr72-12.pdf; ssa.gov's own
# equivalent Actuarial Life Table page, table4c6.html, returns HTTP 403 to automated fetches
# from this environment -- the NCHS report was fetchable and used instead, both being valid
# candidates per the resolutions memo). Only the qx column (probability of dying between exact
# age x and x+1) is used; lx/Lx/Tx/ex are the table's own downstream columns, not independent
# inputs, and are not transcribed. Row "100 and older" is the table's own terminal bucket
# (qx = 1 by the table's own construction, not a modelling choice) -- 101 rows, ages 0-100.
#
# 2021 is the most recent NCHS single-year US life table published as of this work. It reflects
# elevated COVID-era mortality (life expectancy at birth 73.5 (male) / 79.3 (female) in this
# vintage, several years below the pre-pandemic trend) -- a documented, not hidden, property of
# using the most current available data, and a candidate limitation/sensitivity check (a
# pre-pandemic, e.g. 2019, vintage) for a future pass, not implemented here.
#
# ---- Why REPLACE Aliyev's own embedded mortality, not add to it -------------------------------
#
# Aliyev's own maintenance matrices (Table 4, data/raw/aliyev2019_appendixS2_table4_...) already
# carry a flat 0.00006-per-2-week (~0.00156/year) all-cause death probability from every live
# state. That figure is realistic for his trial's actual (young, mean age 35) cohort over a
# short follow-up, but held CONSTANT for a lifetime run it would silently imply a healthy
# 35-year-old's mortality risk never rises with age -- structurally wrong once the horizon
# extends into old age, which is the entire point of extending it. analysis_plan.md's own
# parameter table (item 7) already anticipates this: background mortality should be
# "age- and sex-specific... no CD excess mortality" -- i.e. total mortality = general-population
# background only, not background-plus-Aliyev's-trial-figure. death_prob_schedule() below feeds
# age_adjust_matrix() (R/utils/transition_matrix.R), which REPLACES each row's Death probability
# with the life-table figure for the cohort's current attained age rather than adding it on top
# -- see that function's own header for the exact mechanics. This is a documented modelling
# assumption (no CD-specific excess mortality, same assumption TA456/NICE makes for this
# disease area), not a finding, and should be named as such in Methods, same as the mixture-cure
# structure's own assumptions (analysis_plan.md §15's closing list).
#
# ---- Why two sub-cohorts (male/female), not one blended curve ---------------------------------
#
# This study's population is Aliyev's own cohort: mean age 35, 50% male (analysis_plan.md §5).
# Rather than average male/female qx into one sex-neutral curve (which would misstate both
# sexes' actual risk and can't be labelled "age- and sex-specific" honestly),
# R/02_markov_engine.R's run_maintenance_arm_with_mortality() and
# R/03_cure_fraction_module.R's run_treg_arm_with_mortality() run the SAME transition structure
# twice, once per sex at half the cohort's initial mass, each against that sex's own
# death_prob_schedule(), and sum the two occupancy traces. This is exact (not an approximation):
# cohort-Markov mass is additively separable across independent sub-cohorts run through
# structurally identical matrices that differ only in their Death-column entries, and cost/
# utility attachment (R/04_costs_utilities.R) depends only on health state, not sex, so a summed
# trace feeds the existing attach_*() pipeline completely unchanged.

if (!exists("age_adjust_matrix")) source("R/utils/transition_matrix.R")

#' Load the sourced NCHS 2021 life table (age 0-100, qx by sex; age 100 is the table's own
#' terminal "100 and older" bucket, qx = 1). One row per exact age. Only checks the file's own
#' shape (ages 0..max present in order, terminal qx = 1 for both sexes) -- the qx values
#' themselves are a verbatim published transcription, not independently re-derivable here, same
#' provenance model as R/00_derive_transition_probs.R's own Aliyev tables.
load_life_table <- function(raw_dir = "data/raw") {
  df <- utils::read.csv(file.path(raw_dir, "nchs_us_life_tables_2021.csv"), stringsAsFactors = FALSE)
  stopifnot(
    identical(df$age, 0:(nrow(df) - 1)),
    df$qx_male[nrow(df)] == 1, df$qx_female[nrow(df)] == 1,
    all(df$qx_male >= 0 & df$qx_male <= 1), all(df$qx_female >= 0 & df$qx_female <= 1)
  )
  df
}

#' Convert an annual death probability (qx, the life table's own unit) to the equivalent
#' probability over a shorter period (`cycle_weeks`), under the standard actuarial assumption
#' that deaths are uniformly distributed within the year of age (UDD -- e.g. Bowers et al.,
#' *Actuarial Mathematics*) -- the same closed-form annual-to-period conversion in kind as
#' R/02_markov_engine.R's own discount_factor(), just for a probability instead of a discount
#' rate.
annual_qx_to_cycle_prob <- function(qx, cycle_weeks = 2) {
  stopifnot(qx >= 0, qx <= 1, cycle_weeks > 0)
  1 - (1 - qx)^(cycle_weeks / 52)
}

#' The per-cycle background-mortality probability for each of `n_cycles` transitions, for a
#' cohort starting at `baseline_age` (this study's own cohort mean, analysis_plan.md §5: 35) and
#' a given `sex`. Element t is the probability applied to the transition FROM cycle t-1 TO cycle
#' t, using the life table's qx at the cohort's attained age at the START of that cycle --
#' `floor(baseline_age + (t - 1) * cycle_weeks / 52)`, the standard "age last birthday"
#' actuarial convention. Ages at or beyond the table's terminal row (100) are held at that row's
#' qx = 1 (nobody survives past it, matching the table's own closing convention) -- relevant
#' only once a run reaches that age, e.g. HORIZON_CYCLES_LIFETIME's own derivation
#' (R/05_deterministic_results.R) is sized so the LAST cycle of a lifetime run lands exactly
#' there, fully exhausting the cohort rather than leaving residual undiscounted mass alive
#' forever.
death_prob_schedule <- function(baseline_age, sex, n_cycles, cycle_weeks = 2, life_table = NULL,
                                 raw_dir = "data/raw") {
  sex <- match.arg(sex, c("male", "female"))
  stopifnot(baseline_age >= 0, n_cycles >= 0)
  if (is.null(life_table)) life_table <- load_life_table(raw_dir)
  qx_col <- paste0("qx_", sex)
  max_age <- max(life_table$age)

  vapply(seq_len(n_cycles), function(t) {
    attained_age <- floor(baseline_age + (t - 1) * cycle_weeks / 52)
    attained_age <- min(attained_age, max_age)
    qx <- life_table[[qx_col]][life_table$age == attained_age]
    annual_qx_to_cycle_prob(qx, cycle_weeks)
  }, numeric(1))
}
