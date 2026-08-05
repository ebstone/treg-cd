# US general-population utility norms for scenario S7 (analysis_plan.md §6.2/§7.1 item 21:
# "SDR utility... scenario = age-matched general population") -- the sourced input that
# scenario needs, and that R/04_costs_utilities.R's own docstring flagged as missing
# ("general-population-utility scenario is not implemented here").
#
# ---- Source -----------------------------------------------------------------------------------
#
# Jiang R, Janssen MFB, Pickard AS. "US population norms for the EQ-5D-5L and comparison of norms
# from face-to-face and online samples." Qual Life Res. 2021;30(3):803-816.
# data/raw/jiang2020_us_eq5d5l_population_norms.csv -- transcribed from the paper's face-to-face
# sample (its primary/recommended norms, not the online comparison sample), by age band and sex,
# using the US value set (composite time trade-off). Same role for this scenario as
# nchs_us_life_tables_2021.csv plays for lifetime-horizon mortality: a real, independently
# published, age-stratified US population dataset, not a single invented number.
#
# ---- A stated limitation, not hidden: EQ-5D-5L norms applied to (likely) EQ-5D-3L-anchored
# ---- disease-state utilities ------------------------------------------------------------------
#
# This project's own health-state utilities (load_health_state_utilities(), R/04) trace to
# Aliyev/NICE TA352 (analysis_plan.md §7.1 items 19-20), which are EQ-5D-3L-based, not 5L. Jiang
# et al. 2021 is the best available AGE-STRATIFIED US general-population norm set found, but it
# uses the newer EQ-5D-5L descriptive system and US value set, a genuine (if minor, and widely
# tolerated in practice) version mismatch between the disease-state utilities and this scenario's
# comparator -- the same class of limitation as, e.g., mixing ASP and NADAC pricing bases across
# drugs in this project's own cost data. Stated here, not smoothed over.
#
# ---- Why sex-averaged, not sex-specific --------------------------------------------------------
#
# Jiang et al. report separate male/female norms. This project's cohort is a stated 50/50 male/
# female split (analysis_plan.md §5) but R/03_cure_fraction_module.R's run_treg_arm_with_mortality()
# already sums the two sexes' SDR occupancy into one combined trace before R/04 ever sees it
# (that function's own header) -- splitting SDR utility attachment by sex again downstream would
# require re-threading sex through a boundary this project has deliberately kept clean (R/03
# simulates occupancy, R/04 attaches values). general_population_utility_at_age() below instead
# uses the unweighted mean of the male/female figures for each age band, consistent with the
# cohort's own stated 50/50 assumption -- not sample-size-weighted by Jiang's own (unrelated)
# survey Ns, which would misrepresent this project's cohort composition as Jiang's.

#' Load the age/sex-stratified US EQ-5D-5L population-norm table.
load_population_utility_norms <- function(raw_dir = "data/raw") {
  df <- utils::read.csv(
    file.path(raw_dir, "jiang2020_us_eq5d5l_population_norms.csv"),
    stringsAsFactors = FALSE
  )
  stopifnot(all(df$mean_utility >= -1, df$mean_utility <= 1), all(df$sex %in% c("male", "female")))
  df
}

#' Sex-averaged (unweighted mean of the male/female rows -- see this file's own header for why)
#' general-population utility at a given age, from the age-band table `norms` (default: load via
#' load_population_utility_norms()). `age` is floored to a whole year before the band lookup --
#' the paper's own bands are defined in completed years (e.g. "<25" = ages 0-24 inclusive), the
#' same convention death_prob_schedule() (R/utils/life_table.R) already uses for exactly the same
#' reason: a continuous cycle-by-cycle age (e.g. 24.5, from baseline_age + t*cycle_weeks/52) has
#' to land in a whole-year band, not fall through a gap between two integer-only cutoffs. Ages
#' above the table's own terminal open-ended band (age_high = 999, mirroring
#' nchs_us_life_tables_2021.csv's "100 and older" convention) fall into that band, not an error.
general_population_utility_at_age <- function(age, norms = NULL, raw_dir = "data/raw") {
  stopifnot(age >= 0)
  if (is.null(norms)) norms <- load_population_utility_norms(raw_dir)
  age <- floor(age)
  band <- norms[norms$age_low <= age & age <= norms$age_high, ]
  stopifnot(nrow(band) == 2)  # exactly one male + one female row for the matching band
  mean(band$mean_utility)
}

#' Per-cycle general-population utility schedule for cycles 0..n_cycles (length n_cycles + 1,
#' matching attach_sdr_costs_utilities()'s own on_sdr trace length -- NOT death_prob_schedule()'s
#' length-n_cycles convention, which indexes cycle STEPS 1..n_cycles rather than cycle STATES
#' 0..n_cycles; the two schedules serve different callers with different indexing needs, not an
#' inconsistency). `baseline_age = NULL` (default) returns a length-1 vector at
#' ASSUMED_PATIENT_AGE_YEARS (R/04_costs_utilities.R) -- R's own recycling then applies that one
#' reference-age value to every cycle uniformly, the correct behaviour for the 6.15-year/10-year
#' horizons, which don't track attained age at all (R/05's own module header on baseline_age).
#' Pass an explicit `baseline_age` (the lifetime-horizon caller's own ASSUMED_PATIENT_AGE_YEARS)
#' for a genuinely age-varying, cycle-by-cycle schedule instead.
general_population_utility_schedule <- function(n_cycles, cycle_weeks = 2, baseline_age = NULL,
                                                  norms = NULL, raw_dir = "data/raw") {
  stopifnot(n_cycles >= 0)
  if (is.null(norms)) norms <- load_population_utility_norms(raw_dir)
  if (is.null(baseline_age)) {
    if (!exists("ASSUMED_PATIENT_AGE_YEARS")) stop("ASSUMED_PATIENT_AGE_YEARS not found -- source R/04_costs_utilities.R first")
    return(general_population_utility_at_age(ASSUMED_PATIENT_AGE_YEARS, norms))
  }
  vapply(0:n_cycles, function(c) {
    attained_age <- baseline_age + c * cycle_weeks / 52
    general_population_utility_at_age(attained_age, norms)
  }, numeric(1))
}
