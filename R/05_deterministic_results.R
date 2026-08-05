# Deterministic base case and structural scenarios (analysis_plan.md §10.1, §10.3): tornado on
# incremental NMB, plus the two-way (pi, price) headroom frontier (Aim 4).
#
# FIRST PASS, 2026-08-05. Wires R/00-R/04 together end to end for all four arms (UST/IFX/ADA as
# comparators, TREG as the intervention) and implements the two genuinely load-bearing outputs
# for this paper: a base-case results table, and the (pi, price) headroom frontier (Aim 4) --
# deliberately prioritised over the one-way tornado, which this file does NOT yet implement (see
# "Not yet implemented" below). This mirrors how R/04 itself was built, in incremental passes
# each closing one part rather than the whole module at once.
#
# ---- Why there's no single numeric "Treg base case" ---------------------------------------------
#
# Decision 4 (analysis_plan.md §15) is explicit: no numeric value for the cure fraction (pi) or
# post-cure relapse hazard (h) is proposed anywhere in this project, on purpose -- elicitation
# hasn't run, and assigning a point estimate from an analog (PolTREG, transplant tolerance, Tr1
# CATS1) would be exactly the fabrication the whole redesign exists to avoid. The only pi value
# with an actual number behind it is the Null scenario, pi = 0 (analysis_plan.md §7.2's table:
# "Treg collapses to a one-time-cost UST-equivalent... the honest floor"). run_base_case() below
# therefore reports UST/IFX/ADA's real base-case results alongside Treg AT PI=0 ONLY, labelled as
# the Null floor scenario, not as "the" Treg base case -- do not read it as a claim about Treg's
# expected performance. The Optimistic/Moderate/Pessimistic named scenarios (§7.2) stay
# unimplemented here until elicitation (or the uniform-prior fallback) supplies actual numbers for
# pi and h; headroom_pi_star()/headroom_frontier() below are what carry the paper in the meantime,
# exactly as Decision 4 recommends.
#
# ---- Horizon: 6.15-year, not lifetime ------------------------------------------------------------
#
# analysis_plan.md §4.1 recommends a LIFETIME horizon (6.15-year and 10-year as scenarios). Lifetime
# requires US life-table mortality extrapolation beyond Aliyev's own trial-supported cycles, and no
# life-table data has been sourced anywhere in this repository (checked: no file under data/raw/ or
# data/processed/ contains one). Rather than silently truncate "lifetime" to whatever the transition
# matrices happen to do past their trial-supported range, this pass implements only the two horizons
# that need no unsourced extrapolation: HORIZON_CYCLES_6YR (Aliyev's own 40 x 8-week-cycle span,
# converted to this project's native 2-week cycle) and HORIZON_CYCLES_10YR. A true lifetime run is a
# real, flagged gap -- same class of "not silently done" omission as half-cycle correction (A12,
# docs/model_audit_v6.md) and cyclophosphamide dose (R/04's module header) -- not implemented here.
#
# ---- Aim 5 (external validation against Aliyev's own results): a scope limit worth recording ----
#
# Gate 2's closing criterion (analysis_plan.md §14) is reproducing Aliyev et al. (2019)'s own
# published results. What this repository has transcribed from Aliyev is his TRANSITION-PROBABILITY
# appendix (Supplementary Tables 3/4) only -- not his own cost or utility parameters, and not a
# results table (ICER/QALYs/costs) from his paper. This project's costs are deliberately sourced
# from current CMS pricing, not Aliyev's 2017 drug prices (module header, R/04), so a dollar-for-
# dollar reproduction of his reported ICER was never going to match by construction, independent of
# whether the engine is correct. What CAN be validated purely from data already in this repo is
# Aliyev's own EPIDEMIOLOGICAL trace (state occupancy / life-years by arm, undiscounted, no cost
# layer) -- run_comparator_arm_lifetime()'s qalys_by_cycle output before cost attachment is the
# right object for that comparison, if/when Aliyev's own reported life-years become available to
# check against. Not done in this pass; flagged so it isn't mistaken for already covered by the
# "matrices used unmodified" note in R/00_derive_transition_probs.R (that note is about the INPUTS
# matching Aliyev, not the OUTPUTS being checked against his paper).
#
# ---- Not yet implemented in this pass ------------------------------------------------------------
#
# - The full one-way tornado on every Section 7.1 parameter (§10.1) -- most of those parameters
#   (refractory multipliers, biosimilar re-pricing, half-cycle correction, AE costs) don't have
#   sourced ranges yet either, so a tornado built now would mostly be placeholder bars. A future
#   pass should add tornado_one_way() once more of Gate 3's parameterisation work lands.
# - Discount-rate scenarios (0%/1.5%/5%, §4.1) and the societal-perspective scenario -- the plumbing
#   (annual_rate is already a parameter throughout) supports them; no orchestration wrapper exists
#   yet to run and report them as named scenarios.
# - Refractory-population co-primary scenario (Decision 3) -- needs the sourced multipliers Gate 3
#   is responsible for; not available yet.
# - S1-S12 structural scenarios (§10.3, analysis/run_scenario_analyses.R) beyond the no-cap variant
#   that run_maintenance_arm(apply_cap = FALSE) already gives for free.

if (!exists("MAINTENANCE_STATES")) source("R/utils/transition_matrix.R")
if (!exists("load_published_induction")) source("R/00_derive_transition_probs.R")
if (!exists("run_decision_tree")) source("R/01_decision_tree.R")
if (!exists("run_maintenance_arm")) source("R/02_markov_engine.R")
if (!exists("run_treg_arm")) source("R/03_cure_fraction_module.R")
if (!exists("summarise_arm")) source("R/04_costs_utilities.R")

# ---- Horizons and thresholds ---------------------------------------------------------------------

#' Aliyev's own 40 x 8-week-cycle span (analysis_plan.md §4.1's "current draft" horizon, kept here
#' as a scenario rather than the current draft's mistaken use of it as the ONLY horizon), converted
#' to this project's native 2-week cycle: 40 * 8 / 2 = 160. 320 weeks / 52 = 6.1538 years.
HORIZON_CYCLES_6YR <- 160

#' 10-year horizon scenario (analysis_plan.md §4.1): 10 * 52 / 2 = 260 native 2-week cycles.
HORIZON_CYCLES_10YR <- 260

#' Willingness-to-pay thresholds to report (analysis_plan.md §4.1: "Report $50k, $100k, $150k").
WTP_THRESHOLDS_USD <- c(50000, 100000, 150000)

COMPARATOR_THERAPIES <- c("UST", "IFX", "ADA")

# ---- Transition matrices ---------------------------------------------------------------------

#' Build all four arms' maintenance matrices (UST/IFX/ADA/CT) in one call, so callers needing more
#' than one arm (run_base_case(), headroom_frontier()) don't each re-read and re-validate
#' data/raw/aliyev2019_appendixS2_table4_maintenance_transition_probabilities.csv independently.
build_all_transition_matrices <- function(raw_dir = "data/raw") {
  maintenance <- load_published_maintenance(raw_dir)
  therapies <- c(COMPARATOR_THERAPIES, "CT")
  stats::setNames(
    lapply(therapies, function(tx) {
      build_transition_matrix(maintenance[maintenance$therapy == tx, ], MAINTENANCE_STATES)
    }),
    therapies
  )
}

# ---- Comparator arms (UST/IFX/ADA) ---------------------------------------------------------------

#' Run one biologic comparator arm end to end (induction split -> maintenance Markov -> cost/
#' utility attachment -> lifetime summary), at the given horizon. `matrices` is
#' build_all_transition_matrices()'s output (or a compatible named list with `therapy` and "CT"
#' entries) -- passed in rather than rebuilt here so run_base_case() only pays the CSV-read/
#' validation cost once for however many arms it runs.
#'
#' `utilities = NULL` (default) loads the deterministic base-case values from `proc_dir`, as
#' before. R/06_psa.R passes its own per-draw sampled utility vector here instead -- utilities
#' apply identically to every arm (utility depends only on health state, not which arm a patient
#' is on), so a given PSA draw's utility sample must reach every arm's attachment call, not just
#' Treg's; this parameter is what makes that possible without R/06 reimplementing this function.
#'
#' `induction_data`/`schedule`/`prices = NULL` (default): load from `raw_dir`/`proc_dir`, as
#' before -- a single call pays the CSV-read cost three times (once each). Callers running this
#' function many times over the SAME data, only varying `utilities`/`therapy` (R/06_psa.R's PSA
#' loop, thousands of iterations) can load once and pass the same objects in on every call instead
#' -- purely a performance path, output is identical either way.
run_comparator_arm_lifetime <- function(therapy, n_cycles, matrices, weight_kg = ASSUMED_PATIENT_WEIGHT_KG,
                                         cycle_weeks = 2, annual_rate = 0.03, apply_cap = TRUE,
                                         cap_cycle = 52, raw_dir = "data/raw", proc_dir = "data/processed",
                                         utilities = NULL, induction_data = NULL, schedule = NULL,
                                         prices = NULL) {
  therapy <- match.arg(therapy, COMPARATOR_THERAPIES)
  stopifnot(all(c(therapy, "CT") %in% names(matrices)))

  if (is.null(induction_data)) induction_data <- load_published_induction(raw_dir)
  induction_row <- induction_data[induction_data$therapy == therapy, ]
  split <- run_decision_tree(induction_row)

  arm <- run_maintenance_arm(
    matrices[[therapy]], matrices[["CT"]], split$initial_on_biologic, n_cycles,
    cap_cycle = cap_cycle, apply_cap = apply_cap, initial_on_ct = split$initial_on_ct
  )

  if (is.null(utilities)) utilities <- load_health_state_utilities(proc_dir)
  monitoring_costs <- health_state_monitoring_costs()
  if (is.null(schedule)) schedule <- load_dosing_schedule(raw_dir)
  if (is.null(prices)) prices <- load_drug_prices(raw_dir)

  attached <- attach_maintenance_costs_utilities(
    arm, utilities, monitoring_costs, cycle_weeks, annual_rate,
    therapy = therapy, weight_kg = weight_kg, schedule = schedule, prices = prices
  )
  induction_cost <- induction_drug_cost(therapy, weight_kg, schedule, prices)

  summarise_arm(attached, induction_cost = induction_cost)
}

# ---- Treg arm (function of pi, h, and price) -----------------------------------------------------

#' Treg's one-time dose cost as a function of an EXPLICIT price, not the sourced acquisition figure
#' -- the object analysis_plan.md §9.1's EJP algebra needs (C_treg(P) = C_treg,non-drug + D-bar*P,
#' D-bar = 1 in the single-dose base case). Deliberately duplicates treg_dose_cost()'s arithmetic
#' (R/04_costs_utilities.R) rather than calling it: that function always reads the sourced
#' acquisition figure from data/processed and has no price-override parameter, by design (its own
#' header: it's the fixed, cited base-case cost). Sweeping price for the headroom/EJP frontier is a
#' different question ("what if the price were X") from reporting the sourced base case ("what the
#' price actually is") -- use treg_dose_cost() for the latter, this for the former. Same
#' cyclophosphamide-dose and observation-stay-cost opt-ins as treg_dose_cost(), same reasoning for
#' why both default to 0 (R/04's module header).
treg_price_dependent_dose_cost <- function(price_usd, cyclophosphamide_dose_mg = 0,
                                            observation_stay_cost_usd = 0, raw_dir = "data/raw",
                                            prices = NULL) {
  stopifnot(price_usd >= 0)
  if (is.null(prices)) prices <- load_drug_prices(raw_dir)
  cyclophosphamide_cost <- cyclophosphamide_dose_mg * prices$cyclophosphamide_usd_per_mg
  price_usd + prices$iv_administration_usd + cyclophosphamide_cost + observation_stay_cost_usd
}

#' Run the Treg arm end to end at a given (pi, h, price) triple. Uses UST's own induction split and
#' maintenance matrix for the pre-landmark and non-cured tracks, per R/03_cure_fraction_module.R's
#' module header ("efficacy-equivalent to ustekinumab... there is no separate Treg transition
#' matrix to source") -- not an independent parameter choice made here.
#'
#' `utilities = NULL` (default): see run_comparator_arm_lifetime()'s docstring -- same override
#' mechanism, same reason (R/06_psa.R needs one sampled utility vector applied consistently across
#' every arm within a single PSA draw).
#'
#' `induction_data`/`prices = NULL` (default): same caching mechanism as
#' run_comparator_arm_lifetime()'s equivalent parameters -- load once, pass in on every call, for
#' callers (R/06_psa.R) invoking this many times over unchanged source data.
run_treg_arm_lifetime <- function(n_cycles, pi_sdr, relapse_hazard_annual, price_usd, matrices,
                                   weight_kg = ASSUMED_PATIENT_WEIGHT_KG, cycle_weeks = 2,
                                   annual_rate = 0.03, landmark_cycle = 28, cap_cycle = 52,
                                   apply_cap = TRUE, cyclophosphamide_dose_mg = 0,
                                   observation_stay_cost_usd = 0, raw_dir = "data/raw",
                                   proc_dir = "data/processed", utilities = NULL,
                                   induction_data = NULL, prices = NULL) {
  stopifnot(all(c("UST", "CT") %in% names(matrices)))

  if (is.null(induction_data)) induction_data <- load_published_induction(raw_dir)
  induction_row <- induction_data[induction_data$therapy == "UST", ]
  split <- run_decision_tree(induction_row)

  arm <- run_treg_arm(
    matrices[["UST"]], matrices[["CT"]], split$initial_on_biologic, split$initial_on_ct, n_cycles,
    pi_sdr = pi_sdr, relapse_hazard_annual = relapse_hazard_annual, landmark_cycle = landmark_cycle,
    cap_cycle = cap_cycle, apply_cap = apply_cap
  )

  if (is.null(utilities)) utilities <- load_health_state_utilities(proc_dir)
  monitoring_costs <- health_state_monitoring_costs()
  attached <- attach_treg_costs_utilities(
    arm, utilities, monitoring_costs, cycle_weeks, annual_rate, halve_after_cycle = cap_cycle
  )

  if (is.null(prices)) prices <- load_drug_prices(raw_dir)
  dose_cost <- treg_price_dependent_dose_cost(
    price_usd, cyclophosphamide_dose_mg, observation_stay_cost_usd, prices = prices
  )

  summarise_arm(attached, induction_cost = dose_cost)
}

# ---- Net monetary benefit (analysis_plan.md §9.2's hygiene note: NMB scale, never averaged ICERs) --

#' NMB = WTP * QALYs - cost, for one arm summary (summarise_arm()'s output: needs $qalys/$total_cost).
net_monetary_benefit <- function(arm_summary, wtp_usd) {
  wtp_usd * arm_summary$qalys - arm_summary$total_cost
}

#' The "next-best non-dominated comparator" (analysis_plan.md §9.1) at a given WTP, represented via
#' NMB-maximisation rather than a separately-computed efficiency frontier: at a fixed lambda, the
#' comparator with the highest NMB IS the option a decision-maker would currently choose, which is
#' exactly what "next-best non-dominated" is trying to identify -- ranking by NMB already accounts
#' for (simple and extended) dominance without needing a separate frontier-construction step. This
#' is a documented simplification relative to reporting pairwise ICERs (where dominance has to be
#' identified explicitly before the ICER table makes sense), valid specifically because the
#' decision RULE only needs the max, not the full ICER ladder; R/08_ejp.R's own comparator
#' identification for solving P* should use this same function, not a separate implementation.
best_comparator_nmb <- function(comparator_summaries, wtp_usd) {
  nmbs <- vapply(comparator_summaries, net_monetary_benefit, numeric(1), wtp_usd = wtp_usd)
  list(comparator = names(nmbs)[which.max(nmbs)], nmb = max(nmbs))
}

# ---- Base case (Null scenario for Treg -- see module header) --------------------------------------

#' Deterministic results table: UST/IFX/ADA at their real (only) base case, plus TREG at pi=0 (the
#' Null floor scenario, analysis_plan.md §7.2 -- the only pi value with an actual number behind it;
#' see module header for why there is no other numeric "Treg base case" to report here). Treg's
#' price is the sourced ten-Ham acquisition cost (treg_dose_cost()'s default, module header, R/04)
#' -- the current best estimate of what the product would actually cost, not a price under test;
#' price-as-swept-variable is headroom_frontier()'s job, below.
#'
#' Returns one row per arm: qalys, total_cost, and NMB at each of WTP_THRESHOLDS_USD.
run_base_case <- function(n_cycles = HORIZON_CYCLES_6YR, weight_kg = ASSUMED_PATIENT_WEIGHT_KG,
                           cycle_weeks = 2, annual_rate = 0.03, apply_cap = TRUE, cap_cycle = 52,
                           raw_dir = "data/raw", proc_dir = "data/processed") {
  matrices <- build_all_transition_matrices(raw_dir)

  comparator_summaries <- stats::setNames(
    lapply(COMPARATOR_THERAPIES, function(tx) {
      run_comparator_arm_lifetime(tx, n_cycles, matrices, weight_kg, cycle_weeks, annual_rate,
                                   apply_cap, cap_cycle, raw_dir, proc_dir)
    }),
    COMPARATOR_THERAPIES
  )

  treg_price <- load_treg_dose_acquisition_cost(proc_dir)
  treg_summary <- run_treg_arm_lifetime(
    n_cycles, pi_sdr = 0, relapse_hazard_annual = 0, price_usd = treg_price, matrices = matrices,
    weight_kg = weight_kg, cycle_weeks = cycle_weeks, annual_rate = annual_rate,
    cap_cycle = cap_cycle, apply_cap = apply_cap, raw_dir = raw_dir, proc_dir = proc_dir
  )

  all_summaries <- c(comparator_summaries, list(TREG = treg_summary))

  rows <- lapply(names(all_summaries), function(nm) {
    s <- all_summaries[[nm]]
    nmb_cols <- stats::setNames(
      as.list(vapply(WTP_THRESHOLDS_USD, function(wtp) net_monetary_benefit(s, wtp), numeric(1))),
      # format(..., scientific = FALSE): paste0() left to its own devices renders 100000 as
      # "1e+05", giving a column named "nmb_at_1e.05" once rbind.data.frame() make.names()-sanitises
      # it -- not obviously wrong until you go looking for "nmb_at_100000" and it isn't there.
      paste0("nmb_at_", format(WTP_THRESHOLDS_USD, scientific = FALSE, trim = TRUE))
    )
    c(list(intervention = nm, qalys = s$qalys, total_cost = s$total_cost), nmb_cols)
  })
  do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE))
}

# ---- Headroom frontier (Aim 4) ---------------------------------------------------------------------

#' Minimum durable cure fraction pi* required for Treg to be cost-effective at a given price and
#' WTP threshold (analysis_plan.md §7.2 point 2 / Aim 4), holding the relapse hazard fixed. Treg's
#' NMB is monotone non-decreasing in pi (higher pi shifts more mass into the SDR state, which has
#' the same or better utility and the same or lower cost than the ongoing non-cured Markov trace it
#' replaces -- SDR uses Remission's utility/cost, halved after the cap boundary; nothing about
#' increasing pi ever makes the trace more expensive or less healthy), so a root exists and is
#' unique whenever NMB(pi=0) < target <= NMB(pi=1); uniroot() finds it directly rather than a grid
#' search. Returns:
#'   - `pi_star = 0` (already cost-effective with no cure at all -- report this as "no cure needed
#'     at this price", not as an actual root)
#'   - a value in (0, 1] (the genuine headroom result)
#'   - `NA`, `feasible = FALSE` (not cost-effective at this price even at pi = 1 -- the price is
#'     simply too high for any cure fraction to rescue, and that is itself a reportable finding,
#'     not an error)
headroom_pi_star <- function(price_usd, wtp_usd, relapse_hazard_annual = 0,
                              n_cycles = HORIZON_CYCLES_6YR, matrices = NULL,
                              comparator_summaries = NULL, weight_kg = ASSUMED_PATIENT_WEIGHT_KG,
                              cycle_weeks = 2, annual_rate = 0.03, apply_cap = TRUE, cap_cycle = 52,
                              raw_dir = "data/raw", proc_dir = "data/processed") {
  if (is.null(matrices)) matrices <- build_all_transition_matrices(raw_dir)
  if (is.null(comparator_summaries)) {
    comparator_summaries <- stats::setNames(
      lapply(COMPARATOR_THERAPIES, function(tx) {
        run_comparator_arm_lifetime(tx, n_cycles, matrices, weight_kg, cycle_weeks, annual_rate,
                                     apply_cap, cap_cycle, raw_dir, proc_dir)
      }),
      COMPARATOR_THERAPIES
    )
  }
  target_nmb <- best_comparator_nmb(comparator_summaries, wtp_usd)$nmb

  treg_nmb_at_pi <- function(pi_sdr) {
    s <- run_treg_arm_lifetime(
      n_cycles, pi_sdr, relapse_hazard_annual, price_usd, matrices, weight_kg, cycle_weeks,
      annual_rate, cap_cycle = cap_cycle, apply_cap = apply_cap, raw_dir = raw_dir, proc_dir = proc_dir
    )
    net_monetary_benefit(s, wtp_usd)
  }

  # Boundary comparisons use a tiny absolute tolerance (dollars, not proportion of NMB) -- without
  # it, a price solved from the closed-form algebra in R/08_ejp.R (a different code path computing
  # the same INMB=0 root) can land NMB(pi=1) a few 1e-10 below target_nmb from floating-point
  # rounding alone, which would wrongly report the true pi=1 boundary case as infeasible
  # (test-ejp.R's round-trip check). 1e-6 is far above any such noise and far below any NMB
  # difference anyone would treat as a real result.
  BOUNDARY_TOL_USD <- 1e-6

  nmb_at_0 <- treg_nmb_at_pi(0)
  if (nmb_at_0 >= target_nmb - BOUNDARY_TOL_USD) return(list(pi_star = 0, feasible = TRUE))

  nmb_at_1 <- treg_nmb_at_pi(1)
  if (nmb_at_1 < target_nmb - BOUNDARY_TOL_USD) return(list(pi_star = NA_real_, feasible = FALSE))
  # Within tolerance of target but not quite at/above it (the floating-point-noise case above):
  # clamp so uniroot() sees f(1) == 0 exactly rather than a same-signed tiny negative at both ends,
  # which would otherwise throw "values at end points not of opposite sign" for a price whose true
  # root is pi* = 1.
  if (nmb_at_1 < target_nmb) nmb_at_1 <- target_nmb

  root <- stats::uniroot(function(pi_sdr) {
    if (pi_sdr == 1) return(nmb_at_1 - target_nmb)
    treg_nmb_at_pi(pi_sdr) - target_nmb
  }, interval = c(0, 1), tol = 1e-4)
  list(pi_star = root$root, feasible = TRUE)
}

#' The (pi, price) frontier itself (Aim 4 / analysis_plan.md §10.1's "two-way analysis on (pi,
#' price)... the more informative display for this paper than the tornado"): headroom_pi_star()
#' evaluated across a grid of prices, at one WTP threshold and relapse hazard. Rebuilds the
#' transition matrices and comparator summaries ONCE for the whole grid (both are price/pi-
#' independent), rather than once per price point.
headroom_frontier <- function(price_grid_usd, wtp_usd, relapse_hazard_annual = 0,
                               n_cycles = HORIZON_CYCLES_6YR, weight_kg = ASSUMED_PATIENT_WEIGHT_KG,
                               cycle_weeks = 2, annual_rate = 0.03, apply_cap = TRUE, cap_cycle = 52,
                               raw_dir = "data/raw", proc_dir = "data/processed") {
  matrices <- build_all_transition_matrices(raw_dir)
  comparator_summaries <- stats::setNames(
    lapply(COMPARATOR_THERAPIES, function(tx) {
      run_comparator_arm_lifetime(tx, n_cycles, matrices, weight_kg, cycle_weeks, annual_rate,
                                   apply_cap, cap_cycle, raw_dir, proc_dir)
    }),
    COMPARATOR_THERAPIES
  )

  rows <- lapply(price_grid_usd, function(price) {
    res <- headroom_pi_star(price, wtp_usd, relapse_hazard_annual, n_cycles, matrices,
                             comparator_summaries, weight_kg, cycle_weeks, annual_rate, apply_cap,
                             cap_cycle, raw_dir, proc_dir)
    data.frame(price_usd = price, wtp_usd = wtp_usd, relapse_hazard_annual = relapse_hazard_annual,
               pi_star = res$pi_star, feasible = res$feasible)
  })
  do.call(rbind, rows)
}
