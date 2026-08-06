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
# **h is the exception to that, and was mishandled until 2026-08-06 (peer review's B3,
# docs/model_audit_v6.md A19).** Decision 4's "no numeric value" holds for pi, which genuinely has
# none. It stopped holding for h on 2026-08-05, when the h sweep landed
# DEFAULT_RELAPSE_DURATION_YEARS <- 10 (below) as this project's recorded, PolTREG-anchored
# base-case median SDR duration. But headroom_pi_star()/headroom_frontier() here -- and
# ejp_deterministic()/ejp_frontier() in R/08 -- kept defaulting `relapse_hazard_annual = 0`, so
# every headline number this file produced was computed as though Treg's cure were permanent for
# every patient forever. That is not a neutral placeholder and not "conservative": zero relapse
# hazard is the single MOST favourable assumption available to the intervention
# (R/03_cure_fraction_module.R's duration_to_hazard() docstring says so at length, and said so
# already while these defaults sat at 0 -- the sweep tooling was built as an add-on and the
# functions that actually produce results were never switched over). All four defaults are now
# duration_to_hazard(DEFAULT_RELAPSE_DURATION_YEARS), and R/06_psa.R samples h per draw around the
# same anchor. Every Treg-involving number produced before 2026-08-06 is superseded; the pi=0 rows
# are the sole exception, for the structural reason recorded at run_base_case()'s own Treg call.
#
# ---- Horizon: lifetime, now implemented (2026-08-05) -----------------------------------------
#
# analysis_plan.md §4.1 recommends a LIFETIME horizon (6.15-year and 10-year as comparability
# scenarios). The blocker this module's header used to record here -- "no life-table data has
# been sourced anywhere in this repository" -- is closed: data/raw/nchs_us_life_tables_2021.csv
# (NCHS United States Life Tables, 2021, transcribed from the published PDF; see
# R/utils/life_table.R's own module header for the full sourcing chain and the "why REPLACE, not
# ADD, Aliyev's own trial-cohort mortality" reasoning) plus R/utils/transition_matrix.R's
# age_adjust_matrix() and the death_prob_schedule wiring through R/02/R/03 (both files' own
# headers) give every arm-runner in this file an opt-in `baseline_age` parameter: NULL (default)
# reproduces the exact pre-lifetime-horizon behaviour (fixed matrices, Aliyev's own small trial
# mortality figure, every existing test and result unchanged); ASSUMED_PATIENT_AGE_YEARS
# (R/04_costs_utilities.R: 35, this study's own cohort mean) with `n_cycles =
# HORIZON_CYCLES_LIFETIME` runs the lifetime base case and headroom frontier with age- and
# sex-specific background mortality.
#
# HORIZON_CYCLES_6YR and HORIZON_CYCLES_10YR remain available as the comparability scenarios
# analysis_plan.md §4.1 itself asks for (both still default to `baseline_age = NULL`, i.e.
# Aliyev's own embedded mortality, unchanged from before this pass -- a 6-10 year run at a
# baseline age of 35 is close enough to Aliyev's own trial-supported range that swapping in
# life-table mortality there wouldn't materially change anything and isn't done by default; pass
# `baseline_age` explicitly to any horizon if a mortality-basis-matched comparison is wanted).
#
# What this pass does NOT do: wire baseline_age/lifetime horizon through R/06_psa.R (the PSA),
# R/07_evpi_evppi.R (EVPI/EVPPI), or R/08_ejp.R (EJP) -- those still run at HORIZON_CYCLES_6YR by
# default, deterministic-only for the lifetime horizon in this pass. Extending age_adjust_matrix()
# to every PSA draw (~40,000 arm-runner calls, each now potentially re-deriving up to 1,691
# age-adjusted matrix pairs) is a real performance question this pass hasn't benchmarked, not
# just a wiring exercise -- flagged as the next open item, same "not silently done" standard as
# every other scope boundary in this file's history.
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
#   (refractory multipliers, AE costs) don't have sourced ranges yet either, so a tornado built now
#   would mostly be placeholder bars. A future pass should add tornado_one_way() once more of
#   Gate 3's parameterisation work lands. (Half-cycle correction -- implemented 2026-08-05,
#   R/04_costs_utilities.R's half_cycle_weights() -- isn't a tornado-range input at all, it's a
#   structural on/off toggle now defaulted on everywhere; dropped from this list accordingly.)
# - Discount-rate scenarios (0%/1.5%/5%, §4.1) and the societal-perspective scenario -- the plumbing
#   (annual_rate is already a parameter throughout) supports them; no orchestration wrapper exists
#   yet to run and report them as named scenarios.
# - S1, S2, S4-S12 structural scenarios (§10.3, analysis/run_scenario_analyses.R) beyond the no-cap
#   variant that run_maintenance_arm(apply_cap = FALSE) already gives for free, and beyond S3
#   (refractory-population co-primary scenario, Decision 3) -- implemented 2026-08-05,
#   run_refractory_scenario() below, using UNITI-1-vs-UNITI-2-sourced multipliers
#   (R/utils/refractory_multipliers.R). Surgery-hazard elevation for the refractory population is
#   explicitly NOT part of that implementation -- see that file's own module header.

if (!exists("MAINTENANCE_STATES")) source("R/utils/transition_matrix.R")
if (!exists("load_published_induction")) source("R/00_derive_transition_probs.R")
if (!exists("run_decision_tree")) source("R/01_decision_tree.R")
if (!exists("run_maintenance_arm")) source("R/02_markov_engine.R")
if (!exists("run_treg_arm")) source("R/03_cure_fraction_module.R")
if (!exists("summarise_arm")) source("R/04_costs_utilities.R")
if (!exists("load_refractory_multipliers")) source("R/utils/refractory_multipliers.R")

# ---- Horizons and thresholds ---------------------------------------------------------------------

#' Aliyev's own 40 x 8-week-cycle span (analysis_plan.md §4.1's "current draft" horizon, kept here
#' as a scenario rather than the current draft's mistaken use of it as the ONLY horizon), converted
#' to this project's native 2-week cycle: 40 * 8 / 2 = 160. 320 weeks / 52 = 6.1538 years.
HORIZON_CYCLES_6YR <- 160

#' 10-year horizon scenario (analysis_plan.md §4.1): 10 * 52 / 2 = 260 native 2-week cycles.
HORIZON_CYCLES_10YR <- 260

#' Lifetime horizon (analysis_plan.md §4.1's recommended base case; see module header,
#' "Horizon: lifetime, now implemented"). Sized so the LAST cycle's transition lands the cohort
#' exactly at the life table's terminal age (data/raw/nchs_us_life_tables_2021.csv: age 100,
#' qx = 1) starting from ASSUMED_PATIENT_AGE_YEARS (R/04_costs_utilities.R: 35) -- this fully
#' exhausts the cohort by construction rather than leaving residual undiscounted mass alive
#' forever, so there is no need to run further cycles "just in case."  Derivation:
#' death_prob_schedule()'s attained-age formula (R/utils/life_table.R) is
#' floor(35 + (t-1)*2/52); solving floor(35 + (t-1)*2/52) = 100 for the smallest such t gives
#' t = 1691, i.e. (100 - 35) * 52 / 2 + 1 = 1691 native 2-week cycles = 65 years.
HORIZON_CYCLES_LIFETIME <- 1691

#' Willingness-to-pay thresholds to report (analysis_plan.md §4.1: "Report $50k, $100k, $150k").
WTP_THRESHOLDS_USD <- c(50000, 100000, 150000)

#' Median post-cure Sustained Deep Remission duration T (years), swept for the (pi, T, price)
#' headroom surface (`headroom_frontier_by_duration()` below; `docs/treg-cd_decision_resolutions_
#' 2026-08-05.md` §3.3). 10 years is the stated base-case anchor -- the PolTREG T1D cohort
#' supports a multi-year plateau in a subset out to 7-12 years, but nothing published supports
#' permanence -- swept from a pessimistic 2 years up to Inf (permanent remission, the
#' Ovasave/CATS1-style upper bound; R/03_cure_fraction_module.R's duration_to_hazard()). Converted
#' to an annual hazard via duration_to_hazard() at the point of use, not stored as hazards here,
#' so this grid stays in the units a reader reasons about.
RELAPSE_DURATION_GRID_YEARS <- c(2, 5, 10, 20, Inf)

#' The base-case anchor itself, in the units this project reports (years of median SDR duration).
#' Since 2026-08-06 (B3, docs/model_audit_v6.md A19) this is not merely the centre of the sweep
#' grid but the actual default relapse assumption behind every deterministic headline number:
#' `duration_to_hazard(DEFAULT_RELAPSE_DURATION_YEARS)` is the default `relapse_hazard_annual` of
#' headroom_pi_star()/headroom_frontier() below and of ejp_deterministic()/ejp_frontier()
#' (R/08_ejp.R), and the mean of R/06_psa.R's Gamma prior on h. Before that date all four
#' defaulted to 0 (permanent cure) while this constant sat here unused by any of them -- the
#' defect B3 named. Changing this constant therefore now moves published results, not just one
#' row of a sensitivity table; it is a base-case parameter, and the sweep across
#' RELAPSE_DURATION_GRID_YEARS above is what reports sensitivity to it.
DEFAULT_RELAPSE_DURATION_YEARS <- 10

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

#' The "Markov" half of run_comparator_arm_lifetime(), split out (2026-08-06, alongside the
#' horizon-focus change, README.md's Status section) so a caller that needs the same occupancy
#' trace many times over unchanged inputs -- R/06_psa.R's PSA loop, specifically -- can compute it
#' ONCE rather than once per draw. This is safe because nothing PSA currently samples
#' (utilities, Treg's price, pi) ever reaches this function: transition probabilities aren't
#' Dirichlet-sampled (module header's "Not yet implemented" list, R/06_psa.R), so a comparator
#' therapy's induction split and maintenance matrix -- and therefore its whole occupancy trace --
#' is identical across every PSA draw for a fixed (n_cycles, baseline_age, refractory) call.
#' Returns exactly the `arm` object (on_biologic, on_ct, total)
#' run_comparator_arm_lifetime() feeds to attach_maintenance_costs_utilities() next; this function
#' is a pure extraction, not a behavioural change -- run_comparator_arm_lifetime() below now calls
#' it internally and is otherwise byte-identical to before this split.
simulate_comparator_arm_lifetime <- function(therapy, n_cycles, matrices, cap_cycle = 52, apply_cap = TRUE,
                                              raw_dir = "data/raw", induction_data = NULL,
                                              baseline_age = NULL, life_table = NULL, cycle_weeks = 2,
                                              refractory = FALSE, refractory_multipliers = NULL) {
  therapy <- match.arg(therapy, COMPARATOR_THERAPIES)
  stopifnot(all(c(therapy, "CT") %in% names(matrices)))

  if (is.null(induction_data)) induction_data <- load_published_induction(raw_dir)
  induction_row <- induction_data[induction_data$therapy == therapy, ]

  therapy_matrix <- matrices[[therapy]]
  if (refractory) {
    if (is.null(refractory_multipliers)) refractory_multipliers <- load_refractory_multipliers(raw_dir)
    induction_row <- apply_refractory_multiplier_induction(
      induction_row, refractory_multipliers$induction_response_multiplier,
      refractory_multipliers$induction_remission_multiplier
    )
    therapy_matrix <- apply_refractory_multiplier_maintenance(
      therapy_matrix, refractory_multipliers$maintenance_remission_multiplier_cumulative,
      refractory_multipliers$maintenance_cumulative_n_cycles,
      surgery_hazard_multiplier_cumulative = refractory_multipliers$surgery_hazard_multiplier_cumulative,
      surgery_cumulative_n_cycles = refractory_multipliers$surgery_cumulative_n_cycles
    )
  }
  split <- run_decision_tree(induction_row)

  if (is.null(baseline_age)) {
    run_maintenance_arm(
      therapy_matrix, matrices[["CT"]], split$initial_on_biologic, n_cycles,
      cap_cycle = cap_cycle, apply_cap = apply_cap, initial_on_ct = split$initial_on_ct
    )
  } else {
    run_maintenance_arm_with_mortality(
      therapy_matrix, matrices[["CT"]], split$initial_on_biologic, n_cycles,
      baseline_age = baseline_age, cap_cycle = cap_cycle, apply_cap = apply_cap,
      initial_on_ct = split$initial_on_ct, life_table = life_table, raw_dir = raw_dir,
      cycle_weeks = cycle_weeks
    )
  }
}

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
#'
#' `baseline_age = NULL` (default): the induction+maintenance trace runs on a fixed matrix every
#' cycle, exactly as before this parameter existed (module header; every 6.15-year/10-year caller
#' and test is unaffected -- those horizons deliberately still run on Aliyev's own embedded
#' trial-cohort mortality). Pass a starting age (HORIZON_CYCLES_LIFETIME callers use
#' ASSUMED_PATIENT_AGE_YEARS, R/04_costs_utilities.R) to instead route through
#' run_maintenance_arm_with_mortality() (R/02_markov_engine.R) -- age- and sex-specific
#' background mortality, sourced from `life_table` (R/utils/life_table.R's load_life_table(),
#' loaded from `raw_dir` if not supplied).
#'
#' `half_cycle_correction = TRUE` (default, A12/analysis_plan.md §4.1) -- forwarded to
#' attach_maintenance_costs_utilities(); `FALSE` reproduces the exact pre-2026-08-05 uncorrected
#' behaviour (R/04_costs_utilities.R's trace_qalys()/trace_costs() docstrings), e.g. for a
#' comparability scenario against results computed before this correction existed.
#'
#' `refractory = FALSE` (default): reproduces the exact biologic-naive base case, byte-identical,
#' unaffected by anything in this paragraph. `refractory = TRUE` runs Decision 3 / scenario S3
#' instead (analysis_plan.md §10.3) -- the induction row and the therapy's own maintenance matrix
#' (CT's matrix is left untouched; R/utils/refractory_multipliers.R's module header explains why)
#' are both adjusted via the UNITI-1-vs-UNITI-2-sourced multipliers before the same downstream
#' pipeline runs, so cost/utility attachment, the M-S cap-switch to CT, and mortality all still
#' apply exactly as in the naive case -- only the two therapy-specific matrices differ.
#' `refractory_multipliers = NULL` (default): load once via load_refractory_multipliers(raw_dir);
#' pass a pre-loaded list to avoid re-reading the CSV across many calls (same performance pattern
#' as `utilities`/`schedule`/`prices` above).
#'
#' `perspective = "healthcare_sector"` (default, byte-identical to before this parameter existed):
#' scenario S9 (analysis_plan.md §4.1/§10.3) is `perspective = "societal"`, which swaps in
#' societal_monitoring_costs() (R/04_costs_utilities.R, Manceur et al. 2020-sourced productivity
#' cost added to every living state) in place of health_state_monitoring_costs() -- the only
#' change; no other part of this function's pipeline is affected, since the productivity cost is
#' folded into the same `monitoring_costs` named vector every downstream call already consumes.
run_comparator_arm_lifetime <- function(therapy, n_cycles, matrices, weight_kg = ASSUMED_PATIENT_WEIGHT_KG,
                                         cycle_weeks = 2, annual_rate = 0.03, apply_cap = TRUE,
                                         cap_cycle = 52, raw_dir = "data/raw", proc_dir = "data/processed",
                                         utilities = NULL, induction_data = NULL, schedule = NULL,
                                         prices = NULL, baseline_age = NULL, life_table = NULL,
                                         half_cycle_correction = TRUE, refractory = FALSE,
                                         refractory_multipliers = NULL,
                                         perspective = "healthcare_sector") {
  therapy <- match.arg(therapy, COMPARATOR_THERAPIES)

  arm <- simulate_comparator_arm_lifetime(
    therapy, n_cycles, matrices, cap_cycle = cap_cycle, apply_cap = apply_cap, raw_dir = raw_dir,
    induction_data = induction_data, baseline_age = baseline_age, life_table = life_table,
    cycle_weeks = cycle_weeks, refractory = refractory, refractory_multipliers = refractory_multipliers
  )

  if (is.null(utilities)) utilities <- load_health_state_utilities(proc_dir)
  perspective <- match.arg(perspective, c("healthcare_sector", "societal"))
  monitoring_costs <- if (perspective == "societal") {
    societal_monitoring_costs(cycle_weeks)
  } else {
    health_state_monitoring_costs()
  }
  if (is.null(schedule)) schedule <- load_dosing_schedule(raw_dir)
  if (is.null(prices)) prices <- load_drug_prices(raw_dir)

  attached <- attach_maintenance_costs_utilities(
    arm, utilities, monitoring_costs, cycle_weeks, annual_rate,
    therapy = therapy, weight_kg = weight_kg, schedule = schedule, prices = prices,
    half_cycle_correction = half_cycle_correction
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
#'
#' `baseline_age = NULL` (default): same meaning and same backward-compatibility guarantee as
#' run_comparator_arm_lifetime()'s equivalent parameter -- see its docstring. Routes through
#' run_treg_arm_with_mortality() (R/03_cure_fraction_module.R) when supplied, which applies
#' background mortality to the pre-landmark UST-equivalent track, the post-landmark non-cured
#' track, AND the SDR track itself (R/03's own module header on why SDR needed a death exit at
#' all for a lifetime horizon to be sound).
#'
#' `half_cycle_correction = TRUE` (default) -- same meaning as
#' run_comparator_arm_lifetime()'s equivalent parameter, forwarded to attach_treg_costs_utilities()
#' (which applies it identically to both the pre-landmark/non-cured Markov portion and the SDR
#' portion -- R/04's own docstring on that function).
#'
#' `relapse_destination = "Mild"` (default, matching R/03_cure_fraction_module.R's own default --
#' fully backward compatible): the health state a relapsed SDR patient re-enters. Scenario S11
#' (analysis_plan.md §10.3) is `relapse_destination = "Moderate-Severe Responder"` -- the only
#' other value R/03's `run_treg_arm()`/`run_treg_arm_with_mortality()` already accept; this
#' parameter just exposes their existing argument through this wrapper rather than adding new
#' mechanics.
#'
#' `non_cured_hazard_ratio = 1` (default, the identity transform -- fully backward compatible):
#' scenario S12 (analysis_plan.md §6.2/§10.3) is `non_cured_hazard_ratio < 1`, applying
#' R/03_cure_fraction_module.R's `apply_non_cured_hazard_ratio()` to a LOCAL COPY of UST's matrix
#' before this call's own induction/maintenance run -- `matrices[["UST"]]` itself is never
#' mutated, so the real UST comparator arm (run_comparator_arm_lifetime()) is unaffected even when
#' this function and that one share the same `matrices` list in the same caller.
#'
#' `perspective = "healthcare_sector"` (default) -- same meaning as
#' run_comparator_arm_lifetime()'s equivalent parameter (scenario S9); applies identically here
#' since societal_monitoring_costs() is a property of the health state, not the arm.
#'
#' `sdr_utility_source = "remission"` (default, byte-identical to before this parameter existed):
#' scenario S7 (analysis_plan.md §6.2/§7.1 item 21/§10.3) is `sdr_utility_source =
#' "general_population"`, which computes a general-population utility value/schedule
#' (R/utils/population_utility.R) and passes it to attach_treg_costs_utilities()'s `sdr_utility`
#' argument instead of leaving it NULL (which defaults to Remission's own utility there). Uses
#' `baseline_age` when supplied (a genuinely age-varying, cycle-by-cycle schedule for the lifetime
#' horizon) or a single reference-age value otherwise (the 6.15-year/10-year horizons, which don't
#' track attained age at all) -- population_utility.R's own function handles both cases via one
#' call, no branching needed here.
run_treg_arm_lifetime <- function(n_cycles, pi_sdr, relapse_hazard_annual, price_usd, matrices,
                                   weight_kg = ASSUMED_PATIENT_WEIGHT_KG, cycle_weeks = 2,
                                   annual_rate = 0.03, landmark_cycle = 28, cap_cycle = 52,
                                   apply_cap = TRUE, cyclophosphamide_dose_mg = 0,
                                   observation_stay_cost_usd = 0, raw_dir = "data/raw",
                                   proc_dir = "data/processed", utilities = NULL,
                                   induction_data = NULL, prices = NULL, baseline_age = NULL,
                                   life_table = NULL, half_cycle_correction = TRUE,
                                   relapse_destination = "Mild", non_cured_hazard_ratio = 1,
                                   perspective = "healthcare_sector",
                                   sdr_utility_source = "remission") {
  stopifnot(all(c("UST", "CT") %in% names(matrices)))

  if (is.null(induction_data)) induction_data <- load_published_induction(raw_dir)
  induction_row <- induction_data[induction_data$therapy == "UST", ]
  split <- run_decision_tree(induction_row)

  non_cured_matrix <- if (non_cured_hazard_ratio == 1) {
    matrices[["UST"]]
  } else {
    apply_non_cured_hazard_ratio(matrices[["UST"]], non_cured_hazard_ratio)
  }

  if (is.null(baseline_age)) {
    arm <- run_treg_arm(
      non_cured_matrix, matrices[["CT"]], split$initial_on_biologic, split$initial_on_ct, n_cycles,
      pi_sdr = pi_sdr, relapse_hazard_annual = relapse_hazard_annual, landmark_cycle = landmark_cycle,
      cap_cycle = cap_cycle, apply_cap = apply_cap, relapse_destination = relapse_destination
    )
  } else {
    arm <- run_treg_arm_with_mortality(
      non_cured_matrix, matrices[["CT"]], split$initial_on_biologic, split$initial_on_ct, n_cycles,
      pi_sdr = pi_sdr, relapse_hazard_annual = relapse_hazard_annual, baseline_age = baseline_age,
      landmark_cycle = landmark_cycle, cap_cycle = cap_cycle, apply_cap = apply_cap,
      life_table = life_table, raw_dir = raw_dir, cycle_weeks = cycle_weeks,
      relapse_destination = relapse_destination
    )
  }

  if (is.null(utilities)) utilities <- load_health_state_utilities(proc_dir)
  perspective <- match.arg(perspective, c("healthcare_sector", "societal"))
  monitoring_costs <- if (perspective == "societal") {
    societal_monitoring_costs(cycle_weeks)
  } else {
    health_state_monitoring_costs()
  }

  sdr_utility_source <- match.arg(sdr_utility_source, c("remission", "general_population"))
  sdr_utility <- if (sdr_utility_source == "general_population") {
    general_population_utility_schedule(n_cycles, cycle_weeks, baseline_age, raw_dir = raw_dir)
  } else {
    NULL
  }

  attached <- attach_treg_costs_utilities(
    arm, utilities, monitoring_costs, cycle_weeks, annual_rate, halve_after_cycle = cap_cycle,
    half_cycle_correction = half_cycle_correction, sdr_utility = sdr_utility
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
#' `baseline_age = NULL` (default): 6.15-year/10-year horizons, exactly as before this parameter
#' existed (module header). Pass ASSUMED_PATIENT_AGE_YEARS with `n_cycles = HORIZON_CYCLES_LIFETIME`
#' for the lifetime-horizon base case (analysis/run_full_analysis.R does this explicitly).
#'
#' `comparator_therapies = COMPARATOR_THERAPIES` (default, fully backward compatible): which
#' biologics to run as comparator arms. Scenario S2 (analysis_plan.md §10.3, Decision 2) is
#' `comparator_therapies = c("UST", "IFX")` -- ADA excluded -- without needing a separate
#' function; CT and TREG are unaffected either way (CT is not itself a comparator arm, and TREG's
#' own non-cured track is UST-equivalent regardless of which arms this call reports).
#'
#' `perspective = "healthcare_sector"` and `sdr_utility_source = "remission"` (both default,
#' fully backward compatible): forwarded unchanged to every arm's own
#' run_comparator_arm_lifetime()/run_treg_arm_lifetime() call (see their own docstrings) --
#' scenarios S9 and S7 respectively. Both apply to every arm in the same call (all arms report
#' under the same perspective; `sdr_utility_source` only actually changes TREG's own row, since
#' comparator arms have no SDR state at all).
run_base_case <- function(n_cycles = HORIZON_CYCLES_6YR, weight_kg = ASSUMED_PATIENT_WEIGHT_KG,
                           cycle_weeks = 2, annual_rate = 0.03, apply_cap = TRUE, cap_cycle = 52,
                           raw_dir = "data/raw", proc_dir = "data/processed", baseline_age = NULL,
                           life_table = NULL, comparator_therapies = COMPARATOR_THERAPIES,
                           perspective = "healthcare_sector", sdr_utility_source = "remission") {
  matrices <- build_all_transition_matrices(raw_dir)

  comparator_summaries <- stats::setNames(
    lapply(comparator_therapies, function(tx) {
      run_comparator_arm_lifetime(tx, n_cycles, matrices, weight_kg, cycle_weeks, annual_rate,
                                   apply_cap, cap_cycle, raw_dir, proc_dir,
                                   baseline_age = baseline_age, life_table = life_table,
                                   perspective = perspective)
    }),
    comparator_therapies
  )

  treg_price <- load_treg_dose_acquisition_cost(proc_dir)
  # `relapse_hazard_annual = 0` here is NOT the pre-2026-08-06 permanent-cure default this
  # function's sibling headroom/EJP functions were corrected away from (B3, module header) -- it is
  # STRUCTURALLY INERT at pi_sdr = 0 and deliberately left explicit. run_treg_arm()
  # (R/03_cure_fraction_module.R) seeds the SDR track with `remission_mass * pi_sdr`, which is
  # exactly 0 at pi_sdr = 0, and its per-cycle update only ever multiplies that mass (`dying` and
  # `relapsed` are both on_sdr[t] times something), so on_sdr stays identically 0 for the whole
  # horizon and the relapse probability multiplies nothing. Any h whatsoever gives a bit-identical
  # trace here; a test asserts this rather than leaving it as an argument
  # (tests/testthat/test-deterministic-results.R). Do not "fix" this to the new default -- passing
  # a non-zero h would falsely suggest the Null floor depends on the durability assumption, and
  # would silently make this line's behaviour depend on a constant it is provably immune to.
  treg_summary <- run_treg_arm_lifetime(
    n_cycles, pi_sdr = 0, relapse_hazard_annual = 0, price_usd = treg_price, matrices = matrices,
    weight_kg = weight_kg, cycle_weeks = cycle_weeks, annual_rate = annual_rate,
    cap_cycle = cap_cycle, apply_cap = apply_cap, raw_dir = raw_dir, proc_dir = proc_dir,
    baseline_age = baseline_age, life_table = life_table, perspective = perspective,
    sdr_utility_source = sdr_utility_source
  )

  all_summaries <- c(comparator_summaries, list(TREG = treg_summary))
  summaries_to_results_table(all_summaries)
}

#' Shared row-builder behind run_base_case() and run_refractory_scenario(): one row per named
#' arm summary (intervention, qalys, total_cost, NMB at each WTP threshold). Pulled out of
#' run_base_case() when the refractory scenario needed the identical table shape -- this function
#' itself changed nothing about run_base_case()'s own output, just where the code lives.
summaries_to_results_table <- function(all_summaries, extra_cols = list()) {
  rows <- lapply(names(all_summaries), function(nm) {
    s <- all_summaries[[nm]]
    nmb_cols <- stats::setNames(
      as.list(vapply(WTP_THRESHOLDS_USD, function(wtp) net_monetary_benefit(s, wtp), numeric(1))),
      # format(..., scientific = FALSE): paste0() left to its own devices renders 100000 as
      # "1e+05", giving a column named "nmb_at_1e.05" once rbind.data.frame() make.names()-sanitises
      # it -- not obviously wrong until you go looking for "nmb_at_100000" and it isn't there.
      paste0("nmb_at_", format(WTP_THRESHOLDS_USD, scientific = FALSE, trim = TRUE))
    )
    c(extra_cols, list(intervention = nm, qalys = s$qalys, total_cost = s$total_cost), nmb_cols)
  })
  do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE))
}

# ---- Refractory-population scenario (Decision 3 / S3) ----------------------------------------

#' Scenario S3 (analysis_plan.md §10.3, Decision 3): the same base-case pipeline run twice, once
#' on Aliyev's published biologic-naive matrices (population = "biologic_naive", identical to
#' run_base_case()'s own numbers) and once with UST/IFX/ADA's induction and maintenance matrices
#' adjusted by the UNITI-1-vs-UNITI-2-sourced refractory multipliers
#' (R/utils/refractory_multipliers.R; population = "refractory"). CT and TREG are unaffected in
#' both rows -- CT has no refractory-specific data sourced this pass (module header,
#' R/utils/refractory_multipliers.R), and Treg's own reference clinical programme is already a
#' refractory population by design (analysis_plan.md §3), so nothing about its parameterisation
#' changes between the two rows; its TREG row is included at the same pi=0 floor run_base_case()
#' uses, purely so both scenario's tables are directly NMB-comparable arm-for-arm.
#'
#' Returns one combined table (rbind of the two population's own summaries_to_results_table()
#' output) with a leading `population` column -- deliberately not two separate return values, so
#' a caller (or a written-out CSV) can group_by(population) directly rather than needing two
#' objects kept in sync.
run_refractory_scenario <- function(n_cycles = HORIZON_CYCLES_6YR, weight_kg = ASSUMED_PATIENT_WEIGHT_KG,
                                     cycle_weeks = 2, annual_rate = 0.03, apply_cap = TRUE, cap_cycle = 52,
                                     raw_dir = "data/raw", proc_dir = "data/processed",
                                     baseline_age = NULL, life_table = NULL) {
  matrices <- build_all_transition_matrices(raw_dir)
  refractory_multipliers <- load_refractory_multipliers(raw_dir)
  treg_price <- load_treg_dose_acquisition_cost(proc_dir)

  run_one_population <- function(refractory_flag) {
    comparator_summaries <- stats::setNames(
      lapply(COMPARATOR_THERAPIES, function(tx) {
        run_comparator_arm_lifetime(tx, n_cycles, matrices, weight_kg, cycle_weeks, annual_rate,
                                     apply_cap, cap_cycle, raw_dir, proc_dir,
                                     baseline_age = baseline_age, life_table = life_table,
                                     refractory = refractory_flag,
                                     refractory_multipliers = refractory_multipliers)
      }),
      COMPARATOR_THERAPIES
    )
    # `relapse_hazard_annual = 0`: structurally inert at pi_sdr = 0, same reasoning (and same
    # do-not-"fix"-this warning) as run_base_case()'s own Treg call above.
    treg_summary <- run_treg_arm_lifetime(
      n_cycles, pi_sdr = 0, relapse_hazard_annual = 0, price_usd = treg_price, matrices = matrices,
      weight_kg = weight_kg, cycle_weeks = cycle_weeks, annual_rate = annual_rate,
      cap_cycle = cap_cycle, apply_cap = apply_cap, raw_dir = raw_dir, proc_dir = proc_dir,
      baseline_age = baseline_age, life_table = life_table
    )
    all_summaries <- c(comparator_summaries, list(TREG = treg_summary))
    population_label <- if (refractory_flag) "refractory" else "biologic_naive"
    summaries_to_results_table(all_summaries, extra_cols = list(population = population_label))
  }

  rbind(run_one_population(FALSE), run_one_population(TRUE))
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
#' `baseline_age`/`life_table = NULL` (default): same meaning as
#' run_comparator_arm_lifetime()'s equivalent parameters, threaded through to both the
#' comparator summaries (if `comparator_summaries` isn't already supplied) and every
#' treg_nmb_at_pi() evaluation below.
#'
#' `relapse_hazard_annual` defaults to `duration_to_hazard(DEFAULT_RELAPSE_DURATION_YEARS)` --
#' a median 10 years of drug-free remission, this project's PolTREG-anchored base case -- **not
#' to 0, which is what it defaulted to before 2026-08-06 (B3, docs/model_audit_v6.md A19; see this
#' module's own header)**. h = 0 means the cure never wanes for anyone, ever; it is the most
#' favourable assumption available to Treg, not a cautious one, and pi* computed under it is
#' correspondingly the most optimistic figure this function can return. It remains available by
#' passing `relapse_hazard_annual = 0` explicitly (equivalently `duration_to_hazard(Inf)`), which
#' is exactly what headroom_frontier_by_duration()'s T = Inf row does -- an explicitly-labelled
#' upper bound in a sweep, which is the honest place for it, rather than a silent default.
headroom_pi_star <- function(price_usd, wtp_usd,
                              relapse_hazard_annual = duration_to_hazard(DEFAULT_RELAPSE_DURATION_YEARS),
                              n_cycles = HORIZON_CYCLES_6YR, matrices = NULL,
                              comparator_summaries = NULL, weight_kg = ASSUMED_PATIENT_WEIGHT_KG,
                              cycle_weeks = 2, annual_rate = 0.03, apply_cap = TRUE, cap_cycle = 52,
                              raw_dir = "data/raw", proc_dir = "data/processed", baseline_age = NULL,
                              life_table = NULL) {
  if (is.null(matrices)) matrices <- build_all_transition_matrices(raw_dir)
  if (is.null(comparator_summaries)) {
    comparator_summaries <- stats::setNames(
      lapply(COMPARATOR_THERAPIES, function(tx) {
        run_comparator_arm_lifetime(tx, n_cycles, matrices, weight_kg, cycle_weeks, annual_rate,
                                     apply_cap, cap_cycle, raw_dir, proc_dir,
                                     baseline_age = baseline_age, life_table = life_table)
      }),
      COMPARATOR_THERAPIES
    )
  }
  target_nmb <- best_comparator_nmb(comparator_summaries, wtp_usd)$nmb

  # Cached alongside its NMB so the caller (below, and headroom_frontier()) can report the
  # composite "expected discounted QALY gain per treated patient" quantity (pi* * incremental
  # QALYs over the pi=0 track) `docs/treg-cd_decision_resolutions_2026-08-05.md` §3.3 point 2
  # recommends as the primary headroom axis, without a second simulation pass at the resolved
  # root -- treg_nmb_at_pi() already runs the one Markov trace this needs.
  treg_summary_at_pi <- function(pi_sdr) {
    run_treg_arm_lifetime(
      n_cycles, pi_sdr, relapse_hazard_annual, price_usd, matrices, weight_kg, cycle_weeks,
      annual_rate, cap_cycle = cap_cycle, apply_cap = apply_cap, raw_dir = raw_dir, proc_dir = proc_dir,
      baseline_age = baseline_age, life_table = life_table
    )
  }
  treg_nmb_at_pi <- function(pi_sdr) net_monetary_benefit(treg_summary_at_pi(pi_sdr), wtp_usd)

  # Boundary comparisons use a tiny absolute tolerance (dollars, not proportion of NMB) -- without
  # it, a price solved from the closed-form algebra in R/08_ejp.R (a different code path computing
  # the same INMB=0 root) can land NMB(pi=1) a few 1e-10 below target_nmb from floating-point
  # rounding alone, which would wrongly report the true pi=1 boundary case as infeasible
  # (test-ejp.R's round-trip check). 1e-6 is far above any such noise and far below any NMB
  # difference anyone would treat as a real result.
  BOUNDARY_TOL_USD <- 1e-6

  s0 <- treg_summary_at_pi(0)
  qalys_at_0 <- s0$qalys
  nmb_at_0 <- net_monetary_benefit(s0, wtp_usd)
  if (nmb_at_0 >= target_nmb - BOUNDARY_TOL_USD) {
    return(list(pi_star = 0, feasible = TRUE, qaly_gain = 0))
  }

  s1 <- treg_summary_at_pi(1)
  nmb_at_1 <- net_monetary_benefit(s1, wtp_usd)
  if (nmb_at_1 < target_nmb - BOUNDARY_TOL_USD) {
    return(list(pi_star = NA_real_, feasible = FALSE, qaly_gain = NA_real_))
  }
  # Within tolerance of target but not quite at/above it (the floating-point-noise case above):
  # clamp so uniroot() sees f(1) == 0 exactly rather than a same-signed tiny negative at both ends,
  # which would otherwise throw "values at end points not of opposite sign" for a price whose true
  # root is pi* = 1.
  if (nmb_at_1 < target_nmb) nmb_at_1 <- target_nmb

  root <- stats::uniroot(function(pi_sdr) {
    if (pi_sdr == 1) return(nmb_at_1 - target_nmb)
    treg_nmb_at_pi(pi_sdr) - target_nmb
  }, interval = c(0, 1), tol = 1e-4)
  # One more evaluation at the resolved root -- uniroot() itself only ever sees the NMB
  # difference, not qalys, so this is the first (and only extra) call that actually needs it.
  qaly_gain_at_root <- treg_summary_at_pi(root$root)$qalys - qalys_at_0
  list(pi_star = root$root, feasible = TRUE, qaly_gain = qaly_gain_at_root)
}

#' The (pi, price) frontier itself (Aim 4 / analysis_plan.md §10.1's "two-way analysis on (pi,
#' price)... the more informative display for this paper than the tornado"): headroom_pi_star()
#' evaluated across a grid of prices, at one WTP threshold and relapse hazard. Rebuilds the
#' transition matrices and comparator summaries ONCE for the whole grid (both are price/pi-
#' independent), rather than once per price point.
#' `baseline_age`/`life_table = NULL` (default): same meaning as headroom_pi_star()'s equivalent
#' parameters, threaded through unchanged.
#' `comparator_therapies = COMPARATOR_THERAPIES` (default, backward compatible): same meaning as
#' run_base_case()'s equivalent parameter -- scenario S2 (ADA excluded) passes
#' `c("UST", "IFX")` here to see how the "next-best comparator" (best_comparator_nmb(), used
#' inside headroom_pi_star()) and therefore pi* shift when ADA isn't in the comparator set.
#' `relapse_hazard_annual`: same default and same 2026-08-06 change of it as headroom_pi_star()'s
#' own -- see that function's docstring; this one just forwards whatever it's given, once per
#' price point.
headroom_frontier <- function(price_grid_usd, wtp_usd,
                               relapse_hazard_annual = duration_to_hazard(DEFAULT_RELAPSE_DURATION_YEARS),
                               n_cycles = HORIZON_CYCLES_6YR, weight_kg = ASSUMED_PATIENT_WEIGHT_KG,
                               cycle_weeks = 2, annual_rate = 0.03, apply_cap = TRUE, cap_cycle = 52,
                               raw_dir = "data/raw", proc_dir = "data/processed", baseline_age = NULL,
                               life_table = NULL, comparator_therapies = COMPARATOR_THERAPIES) {
  matrices <- build_all_transition_matrices(raw_dir)
  comparator_summaries <- stats::setNames(
    lapply(comparator_therapies, function(tx) {
      run_comparator_arm_lifetime(tx, n_cycles, matrices, weight_kg, cycle_weeks, annual_rate,
                                   apply_cap, cap_cycle, raw_dir, proc_dir,
                                   baseline_age = baseline_age, life_table = life_table)
    }),
    comparator_therapies
  )

  rows <- lapply(price_grid_usd, function(price) {
    res <- headroom_pi_star(price, wtp_usd, relapse_hazard_annual, n_cycles, matrices,
                             comparator_summaries, weight_kg, cycle_weeks, annual_rate, apply_cap,
                             cap_cycle, raw_dir, proc_dir, baseline_age, life_table)
    data.frame(price_usd = price, wtp_usd = wtp_usd, relapse_hazard_annual = relapse_hazard_annual,
               pi_star = res$pi_star, feasible = res$feasible, qaly_gain = res$qaly_gain)
  })
  do.call(rbind, rows)
}

#' The (pi, T, price) headroom SURFACE (`docs/treg-cd_decision_resolutions_2026-08-05.md` §3.3
#' point 3: "h interacts with the horizon, so [the lifetime horizon] and the h sweep must be run
#' jointly -- at 6.15 years h barely matters; over a lifetime it dominates. A one-way sweep of
#' each in isolation will understate both"). Thin wrapper, not a new simulation: converts each
#' entry of `duration_grid_years` to its annual hazard via
#' R/03_cure_fraction_module.R's duration_to_hazard() and calls headroom_frontier() once per
#' duration, stacking the results with a `duration_years` column so a single data.frame carries
#' the full (T, price) -> (pi*, QALY gain) surface. Intended to be called with
#' `baseline_age`/`life_table` supplied (the lifetime horizon, per the memo's own framing above --
#' the joint-sweep point of this function is close to moot at HORIZON_CYCLES_6YR, where h barely
#' matters regardless of what it's set to).
#'
#' `qaly_gain` (from headroom_pi_star(), each row's exact simulated Treg-incremental-QALY at that
#' row's own pi*) is exactly the "expected discounted QALY gain per treated patient" quantity the
#' memo recommends as the primary reporting axis: because it already collapses whichever (T,
#' price) combination produced it into a single number, the SAME qaly_gain value at different
#' (T, price) points is the honest sense in which "the (pi, T) surface is behind it as a
#' supplementary panel" -- this function returns both, the caller (analysis/run_full_analysis.R)
#' decides which to lead with.
headroom_frontier_by_duration <- function(price_grid_usd, wtp_usd,
                                           duration_grid_years = RELAPSE_DURATION_GRID_YEARS,
                                           n_cycles = HORIZON_CYCLES_6YR,
                                           weight_kg = ASSUMED_PATIENT_WEIGHT_KG, cycle_weeks = 2,
                                           annual_rate = 0.03, apply_cap = TRUE, cap_cycle = 52,
                                           raw_dir = "data/raw", proc_dir = "data/processed",
                                           baseline_age = NULL, life_table = NULL) {
  rows <- lapply(duration_grid_years, function(duration_years) {
    h <- duration_to_hazard(duration_years)
    res <- headroom_frontier(price_grid_usd, wtp_usd, relapse_hazard_annual = h, n_cycles = n_cycles,
                              weight_kg = weight_kg, cycle_weeks = cycle_weeks, annual_rate = annual_rate,
                              apply_cap = apply_cap, cap_cycle = cap_cycle, raw_dir = raw_dir,
                              proc_dir = proc_dir, baseline_age = baseline_age, life_table = life_table)
    res$duration_years <- duration_years
    res
  })
  do.call(rbind, rows)
}

# ---- pi x g(h) factorisation check (resolutions memo §3.3 point 2) ---------------------------

#' Verify, numerically, the resolutions memo's (§3.3 point 2) proposed factorisation: Treg's
#' incremental QALY over its own pi=0 track equals pi * g(h), where g(h) is the discounted QALY
#' value of curing exactly one patient at this relapse hazard (the pi=1-vs-pi=0 gap).
#'
#' **Finding (checked here, not assumed): this is EXACT to floating-point precision, not merely
#' close.** The memo's own framing was cautious -- "very good but not exact," reasoning that
#' relapsed SDR patients re-entering the ordinary Markov trace (R/03_cure_fraction_module.R's own
#' module header) that non-cured Treg patients were already running in might break additivity.
#' It doesn't: run_treg_arm()'s landmark split (`remission_mass * pi_sdr` into SDR,
#' `remission_mass * (1 - pi_sdr)` staying in the ordinary track) is the ONLY place pi enters the
#' whole simulation, and every operation downstream of it -- matrix-vector transitions, the
#' relapse redistribution, the cap-boundary sweep, cost/utility attachment and discounting -- is
#' linear (matrix multiplication and addition throughout, no min/max clipping or other
#' pi-dependent nonlinearity anywhere in the recursion). A linear map applied to a quantity that
#' is itself exactly linear in pi is exactly linear in pi; observed residuals across every (price,
#' horizon, duration) combination tried are ~1e-13 to 1e-16 relative -- ordinary double-precision
#' floating-point noise from accumulating ~1,700 cycles of arithmetic in a different order for
#' each pi, not a sign of genuine model nonlinearity. This function exists to make that check
#' reproducible and visible, not to assert it from algebra alone -- the manuscript should cite
#' this verification, not just the structural argument, since the argument alone was already
#' available before this function existed and the memo still asked for a numeric check.
#'
#' Returns one row per `pi_grid` value: the exact simulated incremental QALY gain over pi=0
#' (`observed_qaly_gain`), the pi * g(h) prediction (`predicted_qaly_gain`), and both the absolute
#' and relative error between them.
verify_pi_factorization <- function(price_usd, relapse_hazard_annual, wtp_usd = WTP_THRESHOLDS_USD[1],
                                     pi_grid = seq(0, 1, by = 0.1), n_cycles = HORIZON_CYCLES_6YR,
                                     weight_kg = ASSUMED_PATIENT_WEIGHT_KG, cycle_weeks = 2,
                                     annual_rate = 0.03, apply_cap = TRUE, cap_cycle = 52,
                                     raw_dir = "data/raw", proc_dir = "data/processed",
                                     baseline_age = NULL, life_table = NULL) {
  matrices <- build_all_transition_matrices(raw_dir)

  qalys_at <- function(pi_sdr) {
    run_treg_arm_lifetime(
      n_cycles, pi_sdr, relapse_hazard_annual, price_usd, matrices, weight_kg, cycle_weeks,
      annual_rate, cap_cycle = cap_cycle, apply_cap = apply_cap, raw_dir = raw_dir, proc_dir = proc_dir,
      baseline_age = baseline_age, life_table = life_table
    )$qalys
  }

  q0 <- qalys_at(0)
  g_h <- qalys_at(1) - q0   # the pi=1-vs-pi=0 gap: the discounted QALY value of one cure at this h

  observed <- vapply(pi_grid, function(pi) qalys_at(pi) - q0, numeric(1))
  predicted <- pi_grid * g_h
  abs_error <- abs(observed - predicted)

  data.frame(
    pi = pi_grid, relapse_hazard_annual = relapse_hazard_annual, g_h = g_h,
    observed_qaly_gain = observed, predicted_qaly_gain = predicted, abs_error = abs_error,
    # Relative error is undefined (0/0) at pi=0, where both gains are exactly 0 by construction --
    # reported as 0, not NaN, since "no error" is the correct characterisation there, not "no
    # data."
    rel_error = ifelse(abs(observed) < 1e-9, 0, abs_error / abs(observed))
  )
}
