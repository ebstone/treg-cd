# Cost and utility attachment by health state and arm (analysis_plan.md §7.1, §8). Combines a
# per-cycle occupancy trace (R/02_markov_engine.R / R/03_cure_fraction_module.R's on_biologic/
# on_ct/on_sdr matrices) with per-cycle costs and utilities to produce discounted QALYs and
# non-drug costs. 2025 USD throughout; Aliyev's native 2-week cycle throughout (no cycle-length
# conversion anywhere in this pipeline -- R/00_derive_transition_probs.R's third-revision header
# comment).
#
# FIRST PASS, 2026-08-04. Implements utility attachment and the arm-independent per-cycle
# health-state monitoring cost in full. UST/IFX drug acquisition + administration costs added
# 2026-08-04 (second pass) -- see "UST/IFX drug acquisition/administration costs" below. ADA and
# Treg's own ten Ham-derived dose cost are still not implemented -- see "Not yet implemented"
# below for why. Non-drug-cost outputs from the first pass ("non_drug_cost_by_cycle" etc.) are
# kept as separately named fields even now that drug costs exist elsewhere in this module --
# don't conflate the two.
#
# ---- Why this doesn't read data/processed/model_health_state_costs.csv ------------------------
#
# That file's non-Surgery rows (Moderate-Severe/M-SR/Mild/Remission) are extracted as-is from the
# v6 Excel workbook (data_dictionary.md: "extracted directly from IBD_CEA_v6_PSA.xlsm") and are
# NOT simply "Aliyev's per-cycle health-state cost" -- they are workbook-era composites that bake
# in an implicit drug-acquisition assumption for the biologic arms. E.g. UST's $14,312 M-SR
# figure decomposes as $282.86 (Aliyev's own state cost, see below) + $14,029.14 -- and
# $14,029.14 = $155.883/mg (data/raw/cms_asp_and_hcup_cost_sources.csv) x 90mg, i.e. a full q8w
# UST maintenance dose, matching the "$14,029 per cycle for ustekinumab" analysis_plan.md §6.4
# cites. That composite was computed under the workbook's original 8-week-cycle design, where one
# cycle *was* one dosing interval. This project now runs a native 2-week cycle: naively applying
# that composite every 2-week cycle would charge a full q8w dose four times as often as it's
# actually given -- a ~4x overstatement of UST/IFX drug cost, not a rounding issue. Same class of
# bug as A9 (Surgery) and the stale cost_per_8wk_cycle_usd_2025 label (both resolved 2026-08-04,
# data/data_dictionary.md), just not caught until now because nothing had consumed this file yet.
#
# The CT track has the mirror-image problem: A8 (docs/model_audit_v6.md) already documents that
# the workbook's flat $88 CT-track figure is the CT *drug* cost only ($67 2017 USD x ~1.306 =
# $87.53), not a health-state cost -- the workbook omits the state-specific monitoring cost A8
# says should be added on top. Reusing $88 as if it were a state cost would carry A8's bug
# forward, not fix it.
#
# So this module treats "cost" as two independently-sourced layers, matching how Aliyev actually
# costed his own model (Appendix S2; analysis_plan.md §8's reconciliation table):
#   1. A monitoring/management cost that varies by health state but is constant across every arm
#      (Aliyev Costs Assumption #1 -- already established for Surgery in the A9 fix). Sourced
#      directly from Aliyev's own reported 2017 USD per-cycle figures below, not from the
#      workbook's arm-specific composites.
#   2. An arm-specific drug-acquisition/administration cost, layered on top only for arms and
#      cycles where a dose is actually given. Implemented for UST/IFX -- see below; ADA and
#      Treg's own dose cost are not.
#
# ---- UST/IFX drug acquisition/administration costs (added 2026-08-04) -------------------------
#
# Sourced from data/raw/biologic_dosing_schedule.csv (dose size, route, and frequency, cited to
# the STELARA/REMICADE FDA-approved Prescribing Information) and
# data/raw/cms_asp_and_hcup_cost_sources.csv ($/mg ASP-based payment limits + the flat IV
# administration fee). Two things this pairing needed that neither file resolves on its own:
#
# 1. **Patient weight**, for UST's induction tier lookup and IFX's mg/kg dosing (both weight-
#    based). Uses ASSUMED_PATIENT_WEIGHT_KG below, Aliyev's own cohort mean (analysis_plan.md
#    §5: "mean age 35, 71 kg, 50% male") -- not manuscript_supplement1_costs_and_discounting.csv's
#    70kg worked example, which is just a rounded illustration of that same figure, not an
#    independent parameter.
# 2. **Vial rounding for IFX**: 71kg x 5mg/kg = 355mg doesn't divide evenly into IFX's 100mg
#    vials (data/raw/aliyev2019_appendixS1_table2_parameters.csv). Rounded UP to the nearest whole
#    vial (400mg), per Aliyev's own Induction Assumption #8 ("No vial sharing... Excess
#    medication is assumed to be wasted," data/raw/aliyev2019_appendixS1_table1_assumptions.csv)
#    -- applied to every IFX dosing event (induction AND maintenance), not just induction, since
#    Remicade vials are single-use/non-preserved and every administration independently wastes
#    any partial vial. UST's induction doses (260/390/520mg) are already exact vial multiples
#    (130mg vial) by design -- no rounding needed there. UST's maintenance dose (90mg) is a fixed
#    pre-filled syringe/pen, not assembled from vials -- also no rounding.
#
# Cost is split the same way the module header already explains: a one-time INDUCTION cost
# (induction_drug_cost(), decision-tree-level, not part of any per-cycle trace -- analysis_plan.md
# §6.1's 8-week induction window) and a recurring MAINTENANCE cost
# (maintenance_drug_cost_by_cycle(), charged only on cycles a real dose is actually given --
# every 4th cycle at this project's native 2-week cycle, i.e. every 8 weeks -- and only against
# mass still on that arm's own biologic track, so patients who've switched to CT or hit the
# 2-year cap stop being charged automatically, with no separate cap-awareness needed here).
#
# **Not applied to Treg's non-cured track.** attach_treg_costs_utilities() reuses UST's own
# maintenance matrix for Treg's non-cured patients (efficacy-equivalent to UST, analysis_plan.md
# §6.2) but that is an efficacy-only equivalence, not a claim that non-cured Treg patients are
# actually receiving ustekinumab. Confirmed with E. Stone, 2026-08-04: non-cured Treg patients
# are NOT charged UST's drug cost -- they received the Treg product (its own one-time/two-dose
# cost, still not wired in, see below), and giving them UST's ongoing drug cost on top would
# double up an acquisition cost they never actually incurred. This is why
# attach_maintenance_costs_utilities() takes an explicit `therapy` argument that must be supplied
# for the drug-cost layer to activate at all -- attach_treg_costs_utilities() deliberately never
# passes one for its Markov portion.
#
# ---- Not yet implemented ------------------------------------------------------------------
#
# ADA's dosing schedule (160mg wk0 + 80mg wk2 induction, 40mg SC every-other-week maintenance --
# note: EOW is a materially different cadence from UST/IFX's every-8-weeks, every native 2-week
# cycle rather than every 4th) is now SOURCED in data/raw/biologic_dosing_schedule.csv (added
# 2026-08-04, same provenance standard as UST/IFX: HUMIRA Prescribing Information via
# humirapro.com, cross-checked against medicalnewstoday.com, and against Aliyev's own 40mg
# syringe unit cost already in this repo). ADA's drug-COST layer is NOT wired into this module
# yet, though -- unlike UST/IFX, ADA is entirely self-administered SC (no Part-B-style ASP price
# applies to any of its doses, including induction; Aliyev's own paper prices it via an FSS unit
# cost instead, data/raw/aliyev2019_appendixS1_table2_parameters.csv). A current $/mg price for
# ADA (NADAC or WAC benchmark, not CMS ASP) still needs sourcing before induction_drug_cost()/
# maintenance_drug_cost_by_cycle() can be extended to it -- flagged here rather than reusing the
# stale 2017 FSS figure or guessing at a current price.
#
# Treg's own acquisition cost (ten Ham-derived, data/processed/model_dose_costs_and_psa_ranges.csv)
# is a decision-tree-level one-time/two-dose event, not a recurring per-cycle charge (R/00's
# third-revision note) -- it belongs in R/01_decision_tree.R's output or a one-time-cost step in
# R/05_deterministic_results.R, not in this per-cycle attachment function, and is still not wired
# in.

if (!exists("MAINTENANCE_STATES")) source("R/utils/transition_matrix.R")
if (!exists("discount_factor")) source("R/02_markov_engine.R")

# ---- Health-state monitoring/management cost (Aliyev Appendix S2; analysis_plan.md §8) --------

#' Aliyev's own per-cycle health-state costs, 2017 USD, at his native 2-week cycle -- constant
#' across every arm (Costs Assumption #1, Aliyev Suppl. Table 1; already established for Surgery
#' in the A9 fix, data/data_dictionary.md, resolved 2026-08-04). M-SR uses the same figure as
#' Moderate-Severe ("assumed equal to M-S", analysis_plan.md §8's reconciliation table) -- Aliyev
#' does not report M-SR separately. Death is always 0.
ALIYEV_MONITORING_COST_2017_USD <- c(
  "Moderate-Severe" = 217, "Moderate-Severe Responder" = 217, "Mild" = 91,
  "Remission" = 10, "Surgery" = 884, "Death" = 0
)

#' 2017 -> 2025 USD inflation factor. Reproduces the workbook's published M-S/Mild/Remission
#' figures exactly (analysis_plan.md §8) but the underlying index is still uncited -- Decision
#' 5's first required fix (analysis_plan.md §15), not resolved by this module.
INFLATION_2017_TO_2025 <- 1.3035

#' Named vector of per-cycle monitoring cost (2025 USD) by health state, applied identically to
#' every arm. Using this uniformly for TREG as well as UST/IFX/CT also resolves A4/Decision 5's
#' second required fix (Treg's health-state costs were un-inflated 2017 dollars) as a side
#' effect of sourcing from Aliyev directly rather than from the arm-specific workbook rows.
health_state_monitoring_costs <- function() {
  ALIYEV_MONITORING_COST_2017_USD * INFLATION_2017_TO_2025
}

# ---- CT-track drug cost (Aliyev Appendix S2 item 18; A8 fix) ----------------------------------

#' Flat per-cycle CT (conventional therapy) drug cost, added on top of the state-specific
#' monitoring cost for every living CT-track state (A8, docs/model_audit_v6.md: the workbook's
#' $88 flat figure IS this drug cost, mislabeled as the whole health-state cost, with the
#' monitoring component missing entirely). $67 (2017 USD) inflated by the same factor as every
#' other line in this module, for internal consistency -- analysis_plan.md §8 reports this one
#' line's factor as "~1.306" (giving $87.53) rather than 1.3035 ($87.33); both are the workbook's
#' own published rounding of the same underlying number, not a discrepancy this module needs to
#' resolve.
CT_DRUG_COST_2017_USD <- 67
ct_drug_cost <- function() CT_DRUG_COST_2017_USD * INFLATION_2017_TO_2025

# ---- UST/IFX drug acquisition + administration cost (data/raw/biologic_dosing_schedule.csv) ----

#' Assumed patient body weight for weight-based dosing (UST induction tier selection, IFX mg/kg
#' conversion). Aliyev's own cohort mean (analysis_plan.md §5: "mean age 35, 71 kg, 50% male"),
#' not the 70kg figure in manuscript_supplement1_costs_and_discounting.csv's worked example,
#' which is a rounded illustration of the same underlying population, not an independent source.
ASSUMED_PATIENT_WEIGHT_KG <- 71

#' IFX vial size (data/raw/aliyev2019_appendixS1_table2_parameters.csv: "FSS IFX Unit Cost (100
#' mg Vial)"). Used for vial-rounding -- see module header.
IFX_VIAL_SIZE_MG <- 100

#' UST/IFX/ADA dosing schedule (data/raw/biologic_dosing_schedule.csv, sourced 2026-08-04:
#' STELARA/REMICADE/HUMIRA Prescribing Information). ADA rows are present but not yet consumed by
#' any cost function in this module -- see "Not yet implemented" in the module header (no current
#' ADA drug price sourced yet).
load_dosing_schedule <- function(raw_dir = "data/raw") {
  utils::read.csv(file.path(raw_dir, "biologic_dosing_schedule.csv"), stringsAsFactors = FALSE)
}

#' $/mg (or $/mg-equivalent) drug prices and the flat IV administration fee
#' (data/raw/cms_asp_and_hcup_cost_sources.csv), keyed by HCPCS code so a code typo fails loudly
#' rather than silently returning NA.
load_drug_prices <- function(raw_dir = "data/raw") {
  df <- utils::read.csv(file.path(raw_dir, "cms_asp_and_hcup_cost_sources.csv"), stringsAsFactors = FALSE)
  price_for <- function(code) {
    val <- df$payment_limit_usd_per_unit[df$hcpcs_code == code]
    stopifnot(length(val) == 1, !is.na(val))
    val
  }
  list(
    ust_induction_usd_per_mg = price_for("J3358"),
    ust_maintenance_usd_per_mg = price_for("J3357"),
    ifx_usd_per_mg = price_for("J1745"),
    iv_administration_usd = price_for("96365")
  )
}

#' UST's induction dose is a discrete weight-tiered label dose, already vial-exact (260/390/520mg
#' are exact 2x/3x/4x multiples of the 130mg vial Aliyev prices induction against) -- a tier
#' lookup against the sourced schedule, not a formula. Tier boundaries read from `schedule`
#' rather than hardcoded here, so drift between this function and the sourced CSV is impossible.
ust_induction_dose_mg <- function(weight_kg, schedule) {
  rows <- schedule[schedule$therapy == "UST" & schedule$phase == "Induction", ]
  lower_ok <- is.na(rows$weight_tier_min_kg) | weight_kg > rows$weight_tier_min_kg
  upper_ok <- is.na(rows$weight_tier_max_kg) | weight_kg <= rows$weight_tier_max_kg
  match <- rows[lower_ok & upper_ok, ]
  stopifnot(nrow(match) == 1)
  match$dose_amount
}

#' IFX's dose is mg/kg, not vial-exact -- round UP to the nearest whole vial (see module header:
#' Aliyev's own "no vial sharing" assumption). Applied identically to induction and maintenance
#' doses (both are single administration events).
ifx_dose_mg <- function(weight_kg, schedule, vial_size_mg = IFX_VIAL_SIZE_MG) {
  mg_per_kg <- unique(schedule$dose_amount[schedule$therapy == "IFX"])
  stopifnot(length(mg_per_kg) == 1)
  raw_mg <- weight_kg * mg_per_kg
  ceiling(raw_mg / vial_size_mg) * vial_size_mg
}

#' One-time induction-phase drug cost (analysis_plan.md §6.1's 8-week decision-tree induction) --
#' NOT part of any per-cycle Markov trace; the caller (eventually R/05_deterministic_results.R)
#' adds this once per patient/arm on top of the per-cycle results below. UST is a single IV dose;
#' IFX is n_doses (3, at weeks 0/2/6, all within the induction window) independently vial-rounded
#' IV doses, each with its own administration fee.
induction_drug_cost <- function(therapy, weight_kg, schedule, prices) {
  therapy <- match.arg(therapy, c("UST", "IFX"))
  if (therapy == "UST") {
    dose_mg <- ust_induction_dose_mg(weight_kg, schedule)
    dose_mg * prices$ust_induction_usd_per_mg + prices$iv_administration_usd
  } else {
    row <- schedule[schedule$therapy == "IFX" & schedule$phase == "Induction", ]
    stopifnot(nrow(row) == 1)
    dose_mg <- ifx_dose_mg(weight_kg, schedule)
    row$n_doses * (dose_mg * prices$ifx_usd_per_mg + prices$iv_administration_usd)
  }
}

#' Per-cycle maintenance drug cost, added on top of on_biologic occupancy at the real dosing
#' cadence read from `schedule` (every 4th cycle at this project's native 2-week cycle = every 8
#' weeks, first dose at maintenance-phase cycle 4 -- data/raw/biologic_dosing_schedule.csv), for as
#' long as a patient remains on that arm's own biologic track. `on_biologic_trace` mass naturally
#' goes to 0 for patients who've switched to CT or hit the 2-year cap (R/02_markov_engine.R), so
#' this needs no separate cap-awareness.
#'
#' UST maintenance is SC (self-injection): drug cost only, no IV administration fee -- Aliyev's
#' own Costs Assumption #3 (data/raw/aliyev2019_appendixS1_table1_assumptions.csv) covers SC
#' patients' visit costs via the monitoring/management cost, not a separate procedure fee. IFX
#' maintenance is IV: drug cost + administration fee every dose, same as induction.
maintenance_drug_cost_by_cycle <- function(on_biologic_trace, therapy, weight_kg, schedule, prices) {
  therapy <- match.arg(therapy, c("UST", "IFX"))
  row <- schedule[schedule$therapy == therapy & schedule$phase == "Maintenance", ]
  stopifnot(nrow(row) == 1)

  n_cycles <- nrow(on_biologic_trace) - 1
  living <- setdiff(colnames(on_biologic_trace), "Death")
  mass_on_biologic <- rowSums(on_biologic_trace[, living, drop = FALSE])

  is_dose_cycle <- rep(FALSE, n_cycles + 1)
  dose_cycles <- seq(row$first_maintenance_cycle_native_2wk, n_cycles,
                      by = row$maintenance_interval_cycles_native_2wk)
  is_dose_cycle[dose_cycles + 1] <- TRUE  # +1: trace row 1 is cycle 0

  if (therapy == "UST") {
    per_dose_cost <- row$dose_amount * prices$ust_maintenance_usd_per_mg  # SC, no admin fee
  } else {
    per_dose_cost <- ifx_dose_mg(weight_kg, schedule) * prices$ifx_usd_per_mg + prices$iv_administration_usd
  }

  mass_on_biologic * is_dose_cycle * per_dose_cost
}

# ---- Utilities (data/processed/model_health_utilities.csv) ------------------------------------

#' Named vector of per-state utility values. Source: model_health_utilities.csv (workbook
#' extraction). Unlike the cost figures above, utilities aren't cycle-length- or
#' drug-cost-entangled (a utility is "how good is this health state," not a per-dosing-interval
#' charge), so no re-derivation concern applies -- safe to use as extracted. Deterministic
#' base-case values only; PSA-time sampling from the underlying ratios
#' (model_utility_ratios_and_psa.csv, analysis_plan.md §7.1 item 20) is R/06_psa.R's job, not
#' this module's.
load_health_state_utilities <- function(proc_dir = "data/processed") {
  df <- utils::read.csv(file.path(proc_dir, "model_health_utilities.csv"), stringsAsFactors = FALSE)
  # Source file uses the workbook's abbreviated state names (Mod/Sev, Mod/Sev Resp); this
  # project's canonical names are MAINTENANCE_STATES (R/00_derive_transition_probs.R). Mapped
  # explicitly, not by row-order alignment, which would silently break if either file's row order
  # ever changed.
  name_map <- c(
    "Mod/Sev" = "Moderate-Severe", "Mod/Sev Resp" = "Moderate-Severe Responder",
    "Mild" = "Mild", "Remission" = "Remission", "Surgery" = "Surgery", "Death" = "Death"
  )
  stopifnot(all(df$health_state %in% names(name_map)))
  stats::setNames(df$utility_value, name_map[df$health_state])
}

# ---- Discounted per-trace QALYs and non-drug costs ---------------------------------------------

#' Per-cycle discount factors for a trace of length n_cycles+1 (cycles 0..n_cycles), reusing
#' R/02_markov_engine.R's discount_factor() rather than recomputing the formula here.
discount_factors_for_trace <- function(n_cycles, cycle_weeks = 2, annual_rate = 0.03) {
  vapply(0:n_cycles, discount_factor, numeric(1), cycle_weeks = cycle_weeks, annual_rate = annual_rate)
}

#' Discounted per-cycle QALYs for a full occupancy trace (rows = cycles 0..n, columns = states).
#' Undiscounted per-cycle utility mass is trace %*% utilities[states] (matches R/02's own
#' trace[t, ] %*% m style), multiplied by cycle length in years and the cycle's discount factor.
#' Half-cycle correction (A12, docs/model_audit_v6.md) is a known, not-yet-implemented gap -- not
#' silently applied here.
trace_qalys <- function(trace, utilities, cycle_weeks = 2, annual_rate = 0.03) {
  states <- colnames(trace)
  stopifnot(!is.null(states), all(states %in% names(utilities)))
  n_cycles <- nrow(trace) - 1
  utility_mass <- as.numeric(trace %*% utilities[states])
  utility_mass * cycle_weeks / 52 * discount_factors_for_trace(n_cycles, cycle_weeks, annual_rate)
}

#' Discounted per-cycle non-drug cost for a full occupancy trace. `add_ct_drug_cost = TRUE` adds
#' the flat CT-track drug cost (see module header) to every living state's occupancy each cycle --
#' pass this for a CT-track trace (e.g. run_maintenance_arm()'s on_ct), not a biologic-track one.
trace_costs <- function(trace, monitoring_costs, cycle_weeks = 2, annual_rate = 0.03,
                         add_ct_drug_cost = FALSE) {
  states <- colnames(trace)
  stopifnot(!is.null(states), all(states %in% names(monitoring_costs)))
  n_cycles <- nrow(trace) - 1
  cost_mass <- as.numeric(trace %*% monitoring_costs[states])
  if (add_ct_drug_cost) {
    living <- setdiff(states, "Death")
    cost_mass <- cost_mass + rowSums(trace[, living, drop = FALSE]) * ct_drug_cost()
  }
  cost_mass * discount_factors_for_trace(n_cycles, cycle_weeks, annual_rate)
}

# ---- Whole-arm attachment ------------------------------------------------------------------

#' Attach costs and QALYs to one run_maintenance_arm()-style result (on_biologic, on_ct, total;
#' R/02_markov_engine.R). Non-drug costs: on_ct gets the CT drug cost added, on_biologic doesn't
#' (its own drug cost, if any, is the separate `drug_cost_by_cycle` field below). QALYs use
#' on_biologic + on_ct (equivalent to `total` for a plain maintenance arm) since utility depends
#' only on health state, not which track a patient is on.
#'
#' `therapy` activates the UST/IFX drug-cost layer (module header) -- pass "UST" or "IFX" plus
#' `weight_kg`/`schedule`/`prices` to charge on_biologic's own maintenance drug cost;
#' `drug_cost_by_cycle` stays all-zero if `therapy` is NULL (the default), which is deliberate for
#' CT-only arms (CT's drug cost is already the separate ct_drug_cost() layer on on_ct, not this
#' one) and for Treg's non-cured track (see attach_treg_costs_utilities() -- never passes
#' `therapy` here, by design, not by omission).
attach_maintenance_costs_utilities <- function(arm_result, utilities, monitoring_costs,
                                                cycle_weeks = 2, annual_rate = 0.03,
                                                therapy = NULL, weight_kg = NULL,
                                                schedule = NULL, prices = NULL) {
  states_trace <- arm_result$on_biologic + arm_result$on_ct
  non_drug_cost_by_cycle <- (
    trace_costs(arm_result$on_biologic, monitoring_costs, cycle_weeks, annual_rate,
                add_ct_drug_cost = FALSE) +
    trace_costs(arm_result$on_ct, monitoring_costs, cycle_weeks, annual_rate,
                add_ct_drug_cost = TRUE)
  )

  n_cycles <- nrow(arm_result$on_biologic) - 1
  if (!is.null(therapy)) {
    stopifnot(!is.null(weight_kg), !is.null(schedule), !is.null(prices))
    drug_cost_by_cycle <- (
      maintenance_drug_cost_by_cycle(arm_result$on_biologic, therapy, weight_kg, schedule, prices) *
      discount_factors_for_trace(n_cycles, cycle_weeks, annual_rate)
    )
  } else {
    drug_cost_by_cycle <- rep(0, n_cycles + 1)
  }

  list(
    qalys_by_cycle = trace_qalys(states_trace, utilities, cycle_weeks, annual_rate),
    non_drug_cost_by_cycle = non_drug_cost_by_cycle,
    drug_cost_by_cycle = drug_cost_by_cycle,
    total_cost_by_cycle = non_drug_cost_by_cycle + drug_cost_by_cycle
  )
}

#' SDR-specific cost and QALYs (analysis_plan.md §6.2). `on_sdr` is a plain numeric vector
#' (cycles 0..n), not a states matrix -- R/03_cure_fraction_module.R's own representation (SDR
#' has no internal state structure). Utility = Remission's utility (base case only; the plan's
#' general-population-utility scenario is not implemented here). Cost = Remission's monitoring
#' cost, halved after `halve_after_cycle` (default the 2-year cap boundary, cycle 52 at this
#' project's native 2-week cycle) -- explicitly "recommend...flag as assumption" in the plan text
#' (§6.2), not an independently sourced figure. `drug_cost_by_cycle` is always all-zero: SDR
#' patients are off therapy by definition, not an omission to fill in later.
attach_sdr_costs_utilities <- function(on_sdr, utilities, monitoring_costs, cycle_weeks = 2,
                                        annual_rate = 0.03, halve_after_cycle = 52) {
  n_cycles <- length(on_sdr) - 1
  cycles <- 0:n_cycles
  discount <- discount_factors_for_trace(n_cycles, cycle_weeks, annual_rate)
  remission_utility <- utilities[["Remission"]]
  remission_cost <- monitoring_costs[["Remission"]]
  cost_rate <- ifelse(cycles > halve_after_cycle, remission_cost / 2, remission_cost)
  non_drug_cost_by_cycle <- on_sdr * cost_rate * discount

  list(
    qalys_by_cycle = on_sdr * remission_utility * cycle_weeks / 52 * discount,
    non_drug_cost_by_cycle = non_drug_cost_by_cycle,
    drug_cost_by_cycle = rep(0, n_cycles + 1),
    total_cost_by_cycle = non_drug_cost_by_cycle
  )
}

#' Attach costs and QALYs to a run_treg_arm()-style result (on_biologic, on_ct, on_sdr, total;
#' R/03_cure_fraction_module.R). Combines attach_maintenance_costs_utilities()'s logic for the
#' Markov portion with attach_sdr_costs_utilities()'s rule for the SDR portion -- arm_result$total
#' isn't used directly here since it's already collapsed to a headcount scalar per cycle
#' (rowSums(on_biologic) + rowSums(on_ct) + on_sdr) and can't be re-split by state.
#'
#' Deliberately never passes `therapy` to the Markov portion below -- see module header ("Not
#' applied to Treg's non-cured track"): non-cured Treg patients are efficacy-equivalent to UST
#' but were never actually given ustekinumab, so they are not charged its drug cost.
attach_treg_costs_utilities <- function(arm_result, utilities, monitoring_costs, cycle_weeks = 2,
                                         annual_rate = 0.03, halve_after_cycle = 52) {
  markov <- attach_maintenance_costs_utilities(
    list(on_biologic = arm_result$on_biologic, on_ct = arm_result$on_ct),
    utilities, monitoring_costs, cycle_weeks, annual_rate
  )
  sdr <- attach_sdr_costs_utilities(
    arm_result$on_sdr, utilities, monitoring_costs, cycle_weeks, annual_rate, halve_after_cycle
  )
  list(
    qalys_by_cycle = markov$qalys_by_cycle + sdr$qalys_by_cycle,
    non_drug_cost_by_cycle = markov$non_drug_cost_by_cycle + sdr$non_drug_cost_by_cycle,
    drug_cost_by_cycle = markov$drug_cost_by_cycle + sdr$drug_cost_by_cycle,
    total_cost_by_cycle = markov$total_cost_by_cycle + sdr$total_cost_by_cycle
  )
}

#' Lifetime (sum of discounted cycles) QALYs and cost for one arm's attach_*() output.
#' `induction_cost` (default 0) is added in undiscounted -- it's already a cycle-0 dollar amount
#' from induction_drug_cost(), and discount_factor(0) = 1 regardless -- for whichever arm's
#' one-time induction drug cost the caller has separately computed; this function only sums
#' per-cycle series, it doesn't call induction_drug_cost() itself.
summarise_arm <- function(attached, induction_cost = 0) {
  list(
    qalys = sum(attached$qalys_by_cycle),
    non_drug_cost = sum(attached$non_drug_cost_by_cycle),
    drug_cost = sum(attached$drug_cost_by_cycle) + induction_cost,
    total_cost = sum(attached$total_cost_by_cycle) + induction_cost
  )
}
