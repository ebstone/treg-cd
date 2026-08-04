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

#' Simulate a biologic maintenance arm (UST/IFX/ADA) alongside the CT track its non-responders
#' and, at the 2-year cap, its long-term survivors switch onto (analysis_plan.md §6.1, §6.4 --
#' Decision 1, recorded 2026-08-04: reinstate the cap for UST/IFX/ADA and non-cured Treg
#' responders; never for cured/SDR Treg).
#'
#' Tracks two parallel occupancy vectors per cycle rather than one compound (arm x state) block
#' matrix, because the Moderate-Severe switch is a *state-conditional* reassignment applied
#' after each matrix multiply, not a fixed transition probability a block matrix could encode
#' directly -- this gives the switching rule as one auditable line rather than an implicit
#' consequence of matrix structure.
#'
#' Every cycle: 1) each track advances under its own matrix; 2) whatever mass lands in
#' Moderate-Severe under the biologic this cycle exits to the CT track's Moderate-Severe bucket
#' ("switched to the CT track at the end of that cycle", §6.1) -- this happens every cycle,
#' cap or no cap; 3) if `apply_cap` and this is `cap_cycle`, ALL remaining biologic-track mass
#' (not just Moderate-Severe) switches to CT wholesale, one time, and the biologic track is zero
#' for the rest of the horizon. `apply_cap = FALSE` gives the no-cap structural scenario
#' (§6.4) for free.
#'
#' Everyone starts on the biologic: `on_ct`'s cycle-0 row is zero regardless of `initial_state`.
run_maintenance_arm <- function(biologic_matrix, ct_matrix, initial_state, n_cycles,
                                 cap_cycle = 52, apply_cap = TRUE) {
  states <- rownames(biologic_matrix)
  stopifnot(
    !is.null(states), identical(states, rownames(ct_matrix)),
    length(initial_state) == length(states), n_cycles >= 0, cap_cycle > 0
  )
  ms_idx <- match("Moderate-Severe", states)
  stopifnot(!is.na(ms_idx))

  on_biologic <- matrix(0, nrow = n_cycles + 1, ncol = length(states), dimnames = list(NULL, states))
  on_ct <- matrix(0, nrow = n_cycles + 1, ncol = length(states), dimnames = list(NULL, states))
  on_biologic[1, ] <- initial_state

  for (t in seq_len(n_cycles)) {
    bio_next <- as.numeric(on_biologic[t, ] %*% biologic_matrix)
    ct_next <- as.numeric(on_ct[t, ] %*% ct_matrix)

    ct_next[ms_idx] <- ct_next[ms_idx] + bio_next[ms_idx]
    bio_next[ms_idx] <- 0

    if (apply_cap && t == cap_cycle) {
      ct_next <- ct_next + bio_next
      bio_next[] <- 0
    }

    on_biologic[t + 1, ] <- bio_next
    on_ct[t + 1, ] <- ct_next
  }

  list(on_biologic = on_biologic, on_ct = on_ct, total = on_biologic + on_ct)
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
