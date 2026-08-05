# Structural scenarios S1-S12 (analysis_plan.md §10.3), first pass, 2026-08-05.
#
# This file covers the scenarios that needed NO new sourcing and little-to-no new mechanics --
# every underlying toggle (apply_cap, comparator_therapies, annual_rate, relapse_destination,
# HORIZON_CYCLES_10YR) already existed in R/05_deterministic_results.R or was a one-line,
# backward-compatible parameter addition there. What follows is orchestration: run the existing
# base case under each named scenario setting and stack the results with a labelled column.
#
# ---- Scenarios covered here -------------------------------------------------------------------
#
# - S1  (2-year maintenance cap on vs. off, Decision 1) -- run_scenario_s1_cap()
# - S2  (ADA included vs. excluded, Decision 2) -- run_scenario_s2_comparator_set()
# - S4  (cure-fraction anchor: named pi points x relapse-duration grid, Decision 4) --
#   run_cure_fraction_anchor_table(). Not new mechanics: the continuous pi sweep
#   (headroom_pi_star()) and the T-duration sweep (headroom_frontier_by_duration()) already exist
#   from the h-sweep work; this reports Treg's OWN outcomes at named (pi, T) points against its
#   SOURCED price, which those two functions don't directly give (they solve for pi* or price*
#   given a target, not "what happens at pi=X").
# - S5  (time horizon: lifetime / 10-year / 6.15-year, analysis_plan.md §4.1/§4.3) --
#   run_scenario_s5_horizon(). Lifetime itself was implemented 2026-08-05 (R/05's own module
#   header); this just stacks all three horizons the base case already supports into one table.
# - S10 (discount rate: 0% / 1.5% / 3% / 5%) -- run_scenario_s10_discount_rate(). annual_rate was
#   already a parameter throughout R/04/R/05; this was purely a missing orchestration wrapper
#   (R/05's own "Not yet implemented" list named exactly this).
# - S11 (SDR relapse re-enters Mild vs. Moderate-Severe Responder) --
#   run_scenario_s11_relapse_destination(). relapse_destination already existed in
#   R/03_cure_fraction_module.R's run_treg_arm()/run_treg_arm_with_mortality(); R/05's
#   run_treg_arm_lifetime() gained a one-line passthrough for it (same commit as this file).
#   Evaluated at the SAME named pi grid as S4 (PI_ANCHOR_GRID) rather than at pi=0: run_base_case()'s
#   own pi=0 Null floor has ZERO patients ever entering SDR, so relapse_destination has
#   (correctly) no effect there at all -- this scenario is only visible once some patients are
#   actually cured, i.e. pi > 0.
# - S12 (non-cured Treg = UST-equivalent vs. HR-advantaged, analysis_plan.md §6.2, added
#   2026-08-05) -- run_scenario_s12_non_cured_hr(). §6.2 itself names the mechanism ("model it
#   explicitly as a hazard ratio on the M-S and Surgery transitions"); R/03_cure_fraction_module.R's
#   new apply_non_cured_hazard_ratio() is exactly that, wired through R/05's run_treg_arm_lifetime()
#   as an opt-in non_cured_hazard_ratio = 1 parameter (identity transform, backward compatible).
#   Evaluated at pi = 0 (unlike S11, this scenario is MOST visible there -- at pi=0, 100% of Treg
#   patients are on the non-cured track this HR modifies). **Genuinely counterintuitive finding,
#   verified not a bug:** QALYs move very slightly in the "wrong" direction as the advantage
#   strengthens (a ~0.01% relative decrease from HR=1 to HR=0.05), because Surgery's own row in
#   Aliyev's published matrix has an 86.7% chance of landing back in Remission the very next cycle
#   -- a faster one-step return to Remission than continuing to cycle through Moderate-Severe
#   Responder or Mild offers. Reducing entry into Surgery therefore isn't unambiguously a QALY win
#   in THIS matrix, even though it is unambiguously a cost win (Surgery is expensive) -- NMB is
#   monotonically higher as the advantage strengthens at every HR tried, since the cost saving
#   dominates the tiny QALY wobble. See R/03's own module header for the full derivation.
#
# ---- Scenarios explicitly NOT covered here, a real flagged gap, not an oversight ---------------
#
# - S3 is done separately (R/utils/refractory_multipliers.R, R/05's run_refractory_scenario()) --
#   not repeated here.
# - S6 (Treg 1 dose vs. 2 doses): treg_dose_cost()'s own module header (R/04) already explains why
#   this needs new per-cycle-trace logic a one-time cost function can't provide -- not built here.
# - S7 (SDR utility = Remission vs. general-population): needs sourcing US general-population
#   utility norms -- not done here (R/04's own docstring already flags this as unimplemented).
# - S9 (healthcare-sector vs. societal perspective): needs the Manceur et al. 2020 productivity/
#   absenteeism cost data actually wired into R/04, not just referenced in analysis_plan.md §4.1 --
#   not done here.

if (!exists("run_base_case")) source("R/05_deterministic_results.R")

#' Named pi anchor points for S4/S11 (Decision 4's own recorded resolution, analysis_plan.md
#' §7.2/§15: "report named scenario points at pi = 0% (null/floor case), 50%, 75%, and 90%").
#' 0% reproduces run_base_case()'s own Null floor exactly (no patients ever cured); the other
#' three are genuinely new evaluations (run_base_case() never runs Treg at pi > 0).
PI_ANCHOR_GRID <- c(0, 0.5, 0.75, 0.9)

# ---- S1: 2-year maintenance cap on vs. off -----------------------------------------------------

#' Stacks run_base_case() with `apply_cap = TRUE` (the actual base case) and `apply_cap = FALSE`
#' (the no-cap variant R/02_markov_engine.R's own docstring already gives "for free" -- this
#' function is the orchestration that was still missing, not new simulation logic). `...` forwards
#' to run_base_case() (e.g. `n_cycles`, `baseline_age` for a non-default horizon).
run_scenario_s1_cap <- function(...) {
  with_cap <- run_base_case(apply_cap = TRUE, ...)
  without_cap <- run_base_case(apply_cap = FALSE, ...)
  rbind(
    cbind(cap = "on", with_cap),
    cbind(cap = "off", without_cap)
  )
}

# ---- S2: ADA included vs. excluded ------------------------------------------------------------

#' Stacks run_base_case() with the full comparator set (UST/IFX/ADA, the actual base case) and
#' with ADA excluded (UST/IFX only, Decision 2's scenario). Both rows use `comparator_therapies`,
#' the one-line backward-compatible parameter run_base_case() gained for exactly this scenario.
run_scenario_s2_comparator_set <- function(...) {
  with_ada <- run_base_case(comparator_therapies = COMPARATOR_THERAPIES, ...)
  without_ada <- run_base_case(comparator_therapies = setdiff(COMPARATOR_THERAPIES, "ADA"), ...)
  rbind(
    cbind(comparator_set = "UST_IFX_ADA", with_ada),
    cbind(comparator_set = "UST_IFX", without_ada)
  )
}

#' The more decision-relevant S2 number: required durable cure fraction pi* at Treg's SOURCED
#' price, under each comparator set, at each WTP threshold -- shows directly whether excluding ADA
#' (i.e. a different "next-best comparator" feeding best_comparator_nmb() inside
#' headroom_pi_star()) changes what cure fraction Treg would need to justify its actual price,
#' which run_scenario_s2_comparator_set()'s NMB table alone doesn't make as legible.
#' `headroom_pi_star()` takes a pre-built `comparator_summaries` list, not a set of therapy
#' names, so each comparator set's summaries are built here (once per set, reused across the WTP
#' grid) rather than passed as a therapy vector.
run_scenario_s2_headroom_at_sourced_price <- function(wtp_grid_usd = WTP_THRESHOLDS_USD,
                                                       n_cycles = HORIZON_CYCLES_6YR,
                                                       weight_kg = ASSUMED_PATIENT_WEIGHT_KG,
                                                       cycle_weeks = 2, annual_rate = 0.03,
                                                       apply_cap = TRUE, cap_cycle = 52,
                                                       raw_dir = "data/raw",
                                                       proc_dir = "data/processed",
                                                       baseline_age = NULL, life_table = NULL) {
  matrices <- build_all_transition_matrices(raw_dir)
  treg_price <- load_treg_dose_acquisition_cost(proc_dir)
  comparator_sets <- list(UST_IFX_ADA = COMPARATOR_THERAPIES,
                           UST_IFX = setdiff(COMPARATOR_THERAPIES, "ADA"))

  rows <- lapply(names(comparator_sets), function(set_name) {
    therapies <- comparator_sets[[set_name]]
    comparator_summaries <- stats::setNames(
      lapply(therapies, function(tx) {
        run_comparator_arm_lifetime(tx, n_cycles, matrices, weight_kg, cycle_weeks, annual_rate,
                                     apply_cap, cap_cycle, raw_dir, proc_dir,
                                     baseline_age = baseline_age, life_table = life_table)
      }),
      therapies
    )
    do.call(rbind, lapply(wtp_grid_usd, function(wtp) {
      res <- headroom_pi_star(treg_price, wtp, n_cycles = n_cycles, matrices = matrices,
                               comparator_summaries = comparator_summaries, weight_kg = weight_kg,
                               cycle_weeks = cycle_weeks, annual_rate = annual_rate,
                               apply_cap = apply_cap, cap_cycle = cap_cycle, raw_dir = raw_dir,
                               proc_dir = proc_dir, baseline_age = baseline_age,
                               life_table = life_table)
      data.frame(comparator_set = set_name, wtp_usd = wtp, price_usd = treg_price,
                 pi_star = res$pi_star, feasible = res$feasible)
    }))
  })
  do.call(rbind, rows)
}

# ---- S4: cure-fraction anchor (named pi points x relapse-duration grid) -----------------------

#' Treg's own outcomes (qalys, total_cost, NMB at each WTP threshold, and whether it beats the
#' best comparator) at every combination of PI_ANCHOR_GRID x RELAPSE_DURATION_GRID_YEARS, at
#' Treg's SOURCED acquisition price -- the "what would actually happen" complement to
#' headroom_pi_star()/headroom_frontier_by_duration()'s "what price/pi would be required" framing.
#' Reports the full grid rather than picking one hand-chosen (pi, T) pairing per named
#' optimistic/moderate/pessimistic label: analysis_plan.md §7.2's three published analogs inform
#' the PLAUSIBLE RANGE, not a single point estimate, so collapsing this to one row per label would
#' misrepresent exactly the uncertainty Decision 4 exists to make visible (same reasoning
#' verify_pi_factorization()/evppi_prior_sensitivity() already apply elsewhere in this codebase).
run_cure_fraction_anchor_table <- function(n_cycles = HORIZON_CYCLES_6YR,
                                            weight_kg = ASSUMED_PATIENT_WEIGHT_KG, cycle_weeks = 2,
                                            annual_rate = 0.03, apply_cap = TRUE, cap_cycle = 52,
                                            raw_dir = "data/raw", proc_dir = "data/processed",
                                            baseline_age = NULL, life_table = NULL,
                                            pi_grid = PI_ANCHOR_GRID,
                                            duration_grid_years = RELAPSE_DURATION_GRID_YEARS) {
  matrices <- build_all_transition_matrices(raw_dir)
  comparator_summaries <- stats::setNames(
    lapply(COMPARATOR_THERAPIES, function(tx) {
      run_comparator_arm_lifetime(tx, n_cycles, matrices, weight_kg, cycle_weeks, annual_rate,
                                   apply_cap, cap_cycle, raw_dir, proc_dir,
                                   baseline_age = baseline_age, life_table = life_table)
    }),
    COMPARATOR_THERAPIES
  )
  treg_price <- load_treg_dose_acquisition_cost(proc_dir)

  grid <- expand.grid(pi_sdr = pi_grid, duration_years = duration_grid_years)
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    pi_sdr <- grid$pi_sdr[i]
    duration_years <- grid$duration_years[i]
    h <- duration_to_hazard(duration_years)
    s <- run_treg_arm_lifetime(
      n_cycles, pi_sdr = pi_sdr, relapse_hazard_annual = h, price_usd = treg_price,
      matrices = matrices, weight_kg = weight_kg, cycle_weeks = cycle_weeks,
      annual_rate = annual_rate, cap_cycle = cap_cycle, apply_cap = apply_cap,
      raw_dir = raw_dir, proc_dir = proc_dir, baseline_age = baseline_age, life_table = life_table
    )
    nmb_cols <- stats::setNames(
      as.list(vapply(WTP_THRESHOLDS_USD, function(wtp) net_monetary_benefit(s, wtp), numeric(1))),
      paste0("nmb_at_", format(WTP_THRESHOLDS_USD, scientific = FALSE, trim = TRUE))
    )
    best_cols <- stats::setNames(
      as.list(vapply(WTP_THRESHOLDS_USD, function(wtp) best_comparator_nmb(comparator_summaries, wtp)$nmb, numeric(1))),
      paste0("best_comparator_nmb_at_", format(WTP_THRESHOLDS_USD, scientific = FALSE, trim = TRUE))
    )
    c(list(pi_sdr = pi_sdr, duration_years = duration_years, relapse_hazard_annual = h,
           price_usd = treg_price, qalys = s$qalys, total_cost = s$total_cost),
      nmb_cols, best_cols)
  })
  do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE))
}

# ---- S5: time horizon (6.15-year / 10-year / lifetime) -----------------------------------------

#' Stacks run_base_case() at all three named horizons analysis_plan.md §4.1/§4.3 discuss: the
#' 6.15-year current-draft horizon, the 10-year comparability horizon (HORIZON_CYCLES_10YR,
#' Aliyev's own embedded trial mortality, `baseline_age = NULL` for both of the first two), and
#' the lifetime horizon (HORIZON_CYCLES_LIFETIME, life-table mortality from ASSUMED_PATIENT_AGE_YEARS
#' -- implemented 2026-08-05, R/05's own module header). Not new simulation -- every one of these
#' three calls already works today; this is the missing side-by-side comparison table.
run_scenario_s5_horizon <- function(...) {
  horizon_6yr <- run_base_case(n_cycles = HORIZON_CYCLES_6YR, ...)
  horizon_10yr <- run_base_case(n_cycles = HORIZON_CYCLES_10YR, ...)
  horizon_lifetime <- run_base_case(n_cycles = HORIZON_CYCLES_LIFETIME,
                                     baseline_age = ASSUMED_PATIENT_AGE_YEARS, ...)
  rbind(
    cbind(horizon = "6.15_year", horizon_6yr),
    cbind(horizon = "10_year", horizon_10yr),
    cbind(horizon = "lifetime", horizon_lifetime)
  )
}

# ---- S10: discount rate (0% / 1.5% / 3% / 5%) --------------------------------------------------

#' Stacks run_base_case() across the four discount rates analysis_plan.md §4.1 names (3% is the
#' base case, already run everywhere else in this codebase; 0%/1.5%/5% are the scenario). Applies
#' to BOTH costs and QALYs identically, matching every other annual_rate call site in this project
#' (R/04_costs_utilities.R's trace_qalys()/trace_costs() docstrings) -- a single differential
#' discount rate for costs vs. QALYs was never part of this project's design and isn't introduced
#' here either.
DISCOUNT_RATE_GRID <- c(0, 0.015, 0.03, 0.05)

run_scenario_s10_discount_rate <- function(discount_rates = DISCOUNT_RATE_GRID, ...) {
  rows <- lapply(discount_rates, function(r) {
    cbind(annual_discount_rate = r, run_base_case(annual_rate = r, ...))
  })
  do.call(rbind, rows)
}

# ---- S11: SDR relapse re-enters Mild vs. Moderate-Severe Responder -----------------------------

#' Treg's own outcomes at PI_ANCHOR_GRID x relapse_destination ∈ {"Mild", "Moderate-Severe
#' Responder"}, at the base-case relapse duration (DEFAULT_RELAPSE_DURATION_YEARS) and Treg's
#' sourced price. Comparator arms are UNAFFECTED by relapse_destination (it's purely a property of
#' Treg's own SDR relapse mechanic, R/03_cure_fraction_module.R), so only Treg's row is reported --
#' unlike run_cure_fraction_anchor_table(), a best-comparator reference column isn't needed here
#' since comparator NMB is identical across both relapse_destination values and would just repeat.
run_scenario_s11_relapse_destination <- function(n_cycles = HORIZON_CYCLES_6YR,
                                                  weight_kg = ASSUMED_PATIENT_WEIGHT_KG,
                                                  cycle_weeks = 2, annual_rate = 0.03,
                                                  apply_cap = TRUE, cap_cycle = 52,
                                                  raw_dir = "data/raw", proc_dir = "data/processed",
                                                  baseline_age = NULL, life_table = NULL,
                                                  pi_grid = PI_ANCHOR_GRID,
                                                  duration_years = DEFAULT_RELAPSE_DURATION_YEARS) {
  matrices <- build_all_transition_matrices(raw_dir)
  treg_price <- load_treg_dose_acquisition_cost(proc_dir)
  h <- duration_to_hazard(duration_years)

  destinations <- c("Mild", "Moderate-Severe Responder")
  grid <- expand.grid(pi_sdr = pi_grid, relapse_destination = destinations,
                       stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    pi_sdr <- grid$pi_sdr[i]
    relapse_destination <- grid$relapse_destination[i]
    s <- run_treg_arm_lifetime(
      n_cycles, pi_sdr = pi_sdr, relapse_hazard_annual = h, price_usd = treg_price,
      matrices = matrices, weight_kg = weight_kg, cycle_weeks = cycle_weeks,
      annual_rate = annual_rate, cap_cycle = cap_cycle, apply_cap = apply_cap,
      raw_dir = raw_dir, proc_dir = proc_dir, baseline_age = baseline_age, life_table = life_table,
      relapse_destination = relapse_destination
    )
    nmb_cols <- stats::setNames(
      as.list(vapply(WTP_THRESHOLDS_USD, function(wtp) net_monetary_benefit(s, wtp), numeric(1))),
      paste0("nmb_at_", format(WTP_THRESHOLDS_USD, scientific = FALSE, trim = TRUE))
    )
    c(list(pi_sdr = pi_sdr, relapse_destination = relapse_destination,
           duration_years = duration_years, qalys = s$qalys, total_cost = s$total_cost),
      nmb_cols)
  })
  do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE))
}

# ---- S12: non-cured Treg = UST-equivalent vs. HR-advantaged --------------------------------------

#' Illustrative hazard-ratio grid for S12 -- 1.0 is the recorded base case (§6.2: "efficacy-
#' equivalent to ustekinumab"); the rest are NOT sourced from any trial (Treg has none) and are
#' explicitly illustrative "what if" values, the same "vary widely" spirit as this project's other
#' assumption-built scenarios (S3, S4).
NON_CURED_HR_GRID <- c(1, 0.9, 0.8, 0.7)

#' Treg's own outcomes at each hazard ratio in NON_CURED_HR_GRID, at pi = 0 (the Null floor --
#' unlike S11, this is where S12 is MOST visible: at pi=0 every Treg patient is on the non-cured
#' track apply_non_cured_hazard_ratio() modifies, not a diminishing (1-pi) fraction of one). Same
#' table shape as run_base_case()'s own Treg row, stacked across the HR grid; comparator rows
#' aren't repeated here since they're identical across every HR value (this parameter only touches
#' the copy of UST's matrix Treg's own non-cured track runs on, R/05's own docstring on
#' run_treg_arm_lifetime()).
run_scenario_s12_non_cured_hr <- function(n_cycles = HORIZON_CYCLES_6YR,
                                           weight_kg = ASSUMED_PATIENT_WEIGHT_KG, cycle_weeks = 2,
                                           annual_rate = 0.03, apply_cap = TRUE, cap_cycle = 52,
                                           raw_dir = "data/raw", proc_dir = "data/processed",
                                           baseline_age = NULL, life_table = NULL,
                                           hazard_ratio_grid = NON_CURED_HR_GRID) {
  matrices <- build_all_transition_matrices(raw_dir)
  treg_price <- load_treg_dose_acquisition_cost(proc_dir)

  rows <- lapply(hazard_ratio_grid, function(hr) {
    s <- run_treg_arm_lifetime(
      n_cycles, pi_sdr = 0, relapse_hazard_annual = 0, price_usd = treg_price, matrices = matrices,
      weight_kg = weight_kg, cycle_weeks = cycle_weeks, annual_rate = annual_rate,
      cap_cycle = cap_cycle, apply_cap = apply_cap, raw_dir = raw_dir, proc_dir = proc_dir,
      baseline_age = baseline_age, life_table = life_table, non_cured_hazard_ratio = hr
    )
    nmb_cols <- stats::setNames(
      as.list(vapply(WTP_THRESHOLDS_USD, function(wtp) net_monetary_benefit(s, wtp), numeric(1))),
      paste0("nmb_at_", format(WTP_THRESHOLDS_USD, scientific = FALSE, trim = TRUE))
    )
    c(list(non_cured_hazard_ratio = hr, qalys = s$qalys, total_cost = s$total_cost), nmb_cols)
  })
  do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE))
}
