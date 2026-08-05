# Cost and utility attachment by health state and arm (analysis_plan.md §7.1, §8). Combines a
# per-cycle occupancy trace (R/02_markov_engine.R / R/03_cure_fraction_module.R's on_biologic/
# on_ct/on_sdr matrices) with per-cycle costs and utilities to produce discounted QALYs and
# non-drug costs. 2025 USD throughout; Aliyev's native 2-week cycle throughout (no cycle-length
# conversion anywhere in this pipeline -- R/00_derive_transition_probs.R's third-revision header
# comment).
#
# FIRST PASS, 2026-08-04. Implements utility attachment and the arm-independent per-cycle
# health-state monitoring cost in full. Deliberately does NOT YET implement drug acquisition/
# administration costs (UST/IFX/ADA induction+maintenance dosing, Treg's ten Ham-derived dose
# cost) -- see "Not yet implemented" below for why. Every cost this file returns is a NON-DRUG
# cost only ("non_drug_cost" is named that way deliberately throughout) -- do not treat it as the
# model's total cost until the drug-cost layer lands.
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
#      cycles where a dose is actually given. NOT YET IMPLEMENTED -- see below.
#
# ---- Not yet implemented: drug acquisition/administration costs -------------------------------
#
# UST/IFX/ADA maintenance and induction drug costs need a dose SIZE and administration FREQUENCY
# per therapy, not just the $/mg unit price already in data/raw/cms_asp_and_hcup_cost_sources.csv.
# That pairing (e.g. "90mg SC every 8 weeks" for UST maintenance) is not cited anywhere in
# data/raw/ -- it exists only implicitly, uncited, inside the legacy workbook composite this
# module deliberately does not reuse (see above). Treg's own acquisition cost (ten Ham-derived,
# data/processed/model_dose_costs_and_psa_ranges.csv) is a decision-tree-level one-time/two-dose
# event, not a recurring per-cycle charge (R/00's third-revision note) -- it belongs in
# R/01_decision_tree.R's output or a one-time-cost step in R/05_deterministic_results.R, not in
# this per-cycle attachment function, and is also not wired in yet.
#
# analysis_plan.md §7.1 items 12-15 already carry "Re-extract"/"Action required" status for
# exactly this reason. Sourcing dose size + frequency per therapy (FDA label or Aliyev's own
# dosing assumptions, cited) is the next step before this module is complete -- flag to
# E. Stone rather than guess at a dosing schedule from outside knowledge not present in this
# repository's data.

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

#' Attach non-drug costs and QALYs to one run_maintenance_arm()-style result (on_biologic, on_ct,
#' total; R/02_markov_engine.R). Costs from on_ct get the CT drug cost added; costs from
#' on_biologic don't (that arm's own drug cost is the separate, not-yet-implemented layer -- see
#' module header). QALYs use on_biologic + on_ct (equivalent to `total` for a plain maintenance
#' arm) since utility depends only on health state, not which track a patient is on.
attach_maintenance_costs_utilities <- function(arm_result, utilities, monitoring_costs,
                                                cycle_weeks = 2, annual_rate = 0.03) {
  states_trace <- arm_result$on_biologic + arm_result$on_ct
  list(
    qalys_by_cycle = trace_qalys(states_trace, utilities, cycle_weeks, annual_rate),
    non_drug_cost_by_cycle = (
      trace_costs(arm_result$on_biologic, monitoring_costs, cycle_weeks, annual_rate,
                  add_ct_drug_cost = FALSE) +
      trace_costs(arm_result$on_ct, monitoring_costs, cycle_weeks, annual_rate,
                  add_ct_drug_cost = TRUE)
    )
  )
}

#' SDR-specific non-drug cost and QALYs (analysis_plan.md §6.2). `on_sdr` is a plain numeric
#' vector (cycles 0..n), not a states matrix -- R/03_cure_fraction_module.R's own representation
#' (SDR has no internal state structure). Utility = Remission's utility (base case only; the
#' plan's general-population-utility scenario is not implemented here). Cost = Remission's
#' monitoring cost, halved after `halve_after_cycle` (default the 2-year cap boundary, cycle 52
#' at this project's native 2-week cycle) -- explicitly "recommend...flag as assumption" in the
#' plan text (§6.2), not an independently sourced figure. No drug cost, ever: SDR patients are
#' off therapy by definition.
attach_sdr_costs_utilities <- function(on_sdr, utilities, monitoring_costs, cycle_weeks = 2,
                                        annual_rate = 0.03, halve_after_cycle = 52) {
  n_cycles <- length(on_sdr) - 1
  cycles <- 0:n_cycles
  discount <- discount_factors_for_trace(n_cycles, cycle_weeks, annual_rate)
  remission_utility <- utilities[["Remission"]]
  remission_cost <- monitoring_costs[["Remission"]]
  cost_rate <- ifelse(cycles > halve_after_cycle, remission_cost / 2, remission_cost)

  list(
    qalys_by_cycle = on_sdr * remission_utility * cycle_weeks / 52 * discount,
    non_drug_cost_by_cycle = on_sdr * cost_rate * discount
  )
}

#' Attach non-drug costs and QALYs to a run_treg_arm()-style result (on_biologic, on_ct, on_sdr,
#' total; R/03_cure_fraction_module.R). Combines attach_maintenance_costs_utilities()'s logic for
#' the Markov portion with attach_sdr_costs_utilities()'s rule for the SDR portion -- arm_result$
#' total isn't used directly here since it's already collapsed to a headcount scalar per cycle
#' (rowSums(on_biologic) + rowSums(on_ct) + on_sdr) and can't be re-split by state.
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
    non_drug_cost_by_cycle = markov$non_drug_cost_by_cycle + sdr$non_drug_cost_by_cycle
  )
}

#' Lifetime (sum of discounted cycles) QALYs and non-drug cost for one arm's attach_*() output.
summarise_arm <- function(attached) {
  list(qalys = sum(attached$qalys_by_cycle), non_drug_cost = sum(attached$non_drug_cost_by_cycle))
}
