# Cohort Markov maintenance engine (analysis_plan.md §6.1, §12.1): a custom matrix-based engine
# (state-occupancy vector multiplied through a list of transition matrices), not heemod, since
# the week-56 cure landmark (R/03) and arm-specific CT-switch rules are awkward in heemod's
# formula interface.
#
# Runs on Aliyev's native 2-week cycle (see R/00_derive_transition_probs.R's third-revision
# header comment for why the project dropped its original 8-week design) using his published
# Supplementary Table 4 matrices unmodified -- no cycle-length conversion happens anywhere in
# this pipeline any more.
#
# This module is deliberately undiscounted and dollar-free: it produces a per-cycle state-
# occupancy trace, nothing else. Discounting/half-cycle correction and cost/utility attachment
# are aggregation-time concerns for R/04_costs_utilities.R and R/05_deterministic_results.R,
# which consume this trace; keeping them out of here is what keeps this file auditable.

if (!exists("validate_row_sums")) source("R/utils/transition_matrix.R")
if (!exists("death_prob_schedule")) source("R/utils/life_table.R")

MAINTENANCE_STATES <- c("Moderate-Severe", "Moderate-Severe Responder", "Mild",
                         "Remission", "Surgery", "Death")

# ---- Core primitive ---------------------------------------------------------

#' Simulate a cohort through `n_cycles` of a single transition matrix, returning the full
#' per-cycle occupancy trace (not just the endpoint -- needed for half-cycle correction
#' downstream, and for the CT-switch bookkeeping in run_maintenance_arm() below).
#'
#' Row 1 of the returned matrix is `initial_state` (cycle 0, before any transition); row i+1 is
#' the occupancy after i cycles. Cohort conservation (each row summing to whatever
#' `initial_state` summed to) follows directly from `m`'s rows summing to 1 -- validated at load
#' time by validate_row_sums(), not re-checked here on every multiply.
simulate_cohort <- function(m, initial_state, n_cycles) {
  states <- rownames(m)
  stopifnot(!is.null(states), length(initial_state) == length(states), n_cycles >= 0)
  trace <- matrix(0, nrow = n_cycles + 1, ncol = length(states), dimnames = list(NULL, states))
  trace[1, ] <- initial_state
  for (t in seq_len(n_cycles)) {
    trace[t + 1, ] <- trace[t, ] %*% m
  }
  trace
}

# ---- Two-track maintenance arm (biologic + CT switch) -----------------------

#' One cycle of the biologic/CT two-track mechanic: advance both tracks under their own
#' matrices, switch whatever mass lands in Moderate-Severe under the biologic this cycle to the
#' CT track's Moderate-Severe bucket ("switched to the CT track at the end of that cycle",
#' §6.1 -- every cycle, cap or no cap), and, if `apply_cap_now`, move ALL remaining
#' biologic-track mass to CT wholesale (the 2-year cap, §6.4 -- Decision 1).
#'
#' Extracted out of run_maintenance_arm() so R/03_cure_fraction_module.R's post-landmark phase
#' can reuse the exact same switch/cap mechanic for Treg's non-cured track, rather than
#' duplicating it -- a second, drifted copy of this logic would be exactly the kind of bug that's
#' easy to introduce and hard to notice.
step_maintenance_cycle <- function(on_biologic_prev, on_ct_prev, biologic_matrix, ct_matrix,
                                    ms_idx, apply_cap_now) {
  bio_next <- as.numeric(on_biologic_prev %*% biologic_matrix)
  ct_next <- as.numeric(on_ct_prev %*% ct_matrix)

  ct_next[ms_idx] <- ct_next[ms_idx] + bio_next[ms_idx]
  bio_next[ms_idx] <- 0

  if (apply_cap_now) {
    ct_next <- ct_next + bio_next
    bio_next[] <- 0
  }

  list(on_biologic = bio_next, on_ct = ct_next)
}

#' Simulate a biologic maintenance arm (UST/IFX/ADA) alongside the CT track its non-responders
#' and, at the 2-year cap, its long-term survivors switch onto (analysis_plan.md §6.1, §6.4 --
#' Decision 1, recorded 2026-08-04: reinstate the cap for UST/IFX/ADA and non-cured Treg
#' responders; never for cured/SDR Treg).
#'
#' Tracks two parallel occupancy vectors per cycle (via step_maintenance_cycle(), above) rather
#' than one compound (arm x state) block matrix, because the Moderate-Severe switch is a
#' *state-conditional* reassignment applied after each matrix multiply, not a fixed transition
#' probability a block matrix could encode directly -- this gives the switching rule as one
#' auditable function rather than an implicit consequence of matrix structure.
#'
#' `apply_cap = FALSE` gives the no-cap structural scenario (§6.4) for free.
#'
#' `initial_on_ct` seeds `on_ct`'s cycle-0 row (defaults to zero: "everyone starts on the
#' biologic"). The default is only correct on its own for a placeholder/synthetic initial
#' state; a real induction split (R/01_decision_tree.R) has non-responders entering directly on
#' CT at cycle 0, not accumulating there only through the cycle-by-cycle Moderate-Severe switch
#' -- pass that split's `initial_on_ct` here to represent that correctly.
#'
#' `death_prob_schedule = NULL` (default): both matrices are used fixed, every cycle, exactly as
#' before this parameter existed -- fully backward compatible, every pre-lifetime-horizon caller
#' and test is unaffected. Pass a length-`n_cycles` vector (R/utils/life_table.R's
#' death_prob_schedule()) to instead re-derive both matrices EVERY cycle via age_adjust_matrix()
#' (R/utils/transition_matrix.R) at that cycle's background-mortality probability -- the
#' lifetime-horizon path (R/05_deterministic_results.R's HORIZON_CYCLES_LIFETIME). Not called
#' directly for that path in practice -- run_maintenance_arm_with_mortality() below is the
#' entry point that actually builds the schedule (per sex) and calls this.
run_maintenance_arm <- function(biologic_matrix, ct_matrix, initial_state, n_cycles,
                                 cap_cycle = 52, apply_cap = TRUE, initial_on_ct = NULL,
                                 death_prob_schedule = NULL) {
  states <- rownames(biologic_matrix)
  stopifnot(
    !is.null(states), identical(states, rownames(ct_matrix)),
    length(initial_state) == length(states), n_cycles >= 0, cap_cycle > 0,
    is.null(death_prob_schedule) || length(death_prob_schedule) == n_cycles
  )
  ms_idx <- match("Moderate-Severe", states)
  stopifnot(!is.na(ms_idx))
  if (is.null(initial_on_ct)) initial_on_ct <- rep(0, length(states))
  stopifnot(length(initial_on_ct) == length(states))

  on_biologic <- matrix(0, nrow = n_cycles + 1, ncol = length(states), dimnames = list(NULL, states))
  on_ct <- matrix(0, nrow = n_cycles + 1, ncol = length(states), dimnames = list(NULL, states))
  on_biologic[1, ] <- initial_state
  on_ct[1, ] <- initial_on_ct

  for (t in seq_len(n_cycles)) {
    if (is.null(death_prob_schedule)) {
      bio_m <- biologic_matrix
      ct_m <- ct_matrix
    } else {
      bio_m <- age_adjust_matrix(biologic_matrix, death_prob_schedule[t])
      ct_m <- age_adjust_matrix(ct_matrix, death_prob_schedule[t])
    }
    step <- step_maintenance_cycle(
      on_biologic[t, ], on_ct[t, ], bio_m, ct_m, ms_idx,
      apply_cap_now = apply_cap && t == cap_cycle
    )
    on_biologic[t + 1, ] <- step$on_biologic
    on_ct[t + 1, ] <- step$on_ct
  }

  list(on_biologic = on_biologic, on_ct = on_ct, total = on_biologic + on_ct)
}

# ---- Lifetime horizon: age- and sex-specific background mortality (R/utils/life_table.R) ----

#' Run run_maintenance_arm() twice, once per sex, at half the cohort's initial mass each (this
#' study's own 50% male / 50% female population, analysis_plan.md §5), each against that sex's
#' own death_prob_schedule(), and sum the two occupancy traces. This is the only place sex
#' enters the model -- see R/utils/life_table.R's module header for why summing two sub-cohorts
#' is exact here, not an approximation, and why background mortality REPLACES rather than adds
#' to Aliyev's own trial-cohort mortality figure.
run_maintenance_arm_with_mortality <- function(biologic_matrix, ct_matrix, initial_state, n_cycles,
                                                baseline_age, cap_cycle = 52, apply_cap = TRUE,
                                                initial_on_ct = NULL, life_table = NULL,
                                                raw_dir = "data/raw", cycle_weeks = 2) {
  states <- rownames(biologic_matrix)
  if (is.null(initial_on_ct)) initial_on_ct <- rep(0, length(states))
  if (is.null(life_table)) life_table <- load_life_table(raw_dir)

  runs <- lapply(c("male", "female"), function(sex) {
    schedule <- death_prob_schedule(baseline_age, sex, n_cycles, cycle_weeks, life_table)
    run_maintenance_arm(
      biologic_matrix, ct_matrix, initial_state * 0.5, n_cycles,
      cap_cycle = cap_cycle, apply_cap = apply_cap, initial_on_ct = initial_on_ct * 0.5,
      death_prob_schedule = schedule
    )
  })

  list(
    on_biologic = runs[[1]]$on_biologic + runs[[2]]$on_biologic,
    on_ct = runs[[1]]$on_ct + runs[[2]]$on_ct,
    total = runs[[1]]$total + runs[[2]]$total
  )
}

# ---- Discounting helper -------------------------------------------------------

#' Discount factor for a cycle `n_cycles` from baseline, given the cycle length in weeks and an
#' annual discount rate (analysis_plan.md §4.1: 3% base case, 0%/1.5%/5% scenarios). A minimal,
#' pure helper -- full QALY/cost discounting (including half-cycle correction) is
#' R/04_costs_utilities.R and R/05_deterministic_results.R's job; this is just the closed-form
#' piece those steps will multiply by.
discount_factor <- function(cycle, cycle_weeks = 2, annual_rate = 0.03) {
  stopifnot(cycle >= 0, cycle_weeks > 0, annual_rate >= 0)
  years <- cycle * cycle_weeks / 52
  1 / (1 + annual_rate)^years
}
