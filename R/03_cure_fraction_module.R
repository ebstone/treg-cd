# Mixture-cure extension (analysis_plan.md §6.2, §6.3): week-56 landmark split into Sustained
# Deep Remission (SDR), conditional on being in Remission under the biologic track at that
# point; constant annual relapse hazard h from SDR back into the ordinary Markov. Cured
# patients never enter CT.
#
# "Non-cured Treg patients: efficacy-equivalent to ustekinumab... drop the 20% adverse-transition
# reduction entirely" (§6.2) means Treg's non-cured population uses UST's own induction split
# (R/01_decision_tree.R) and UST's own maintenance matrix (R/02_markov_engine.R) unmodified --
# there is no separate Treg transition matrix to source. This module is specifically the
# cure-branching layer grafted onto that existing UST pipeline: one landmark-conditional split,
# one new (simple) decaying state, one re-entry rule.
#
# Landmark: week 56. At this project's native 2-week cycle that's cycle 28 (56/2), not "cycle 7"
# as the plan text's original 8-week-cycle design had it (56/8) -- recomputed here for the
# 2-week-cycle switch (R/00_derive_transition_probs.R's third-revision header comment).
#
# 2-year cap interaction (confirmed 2026-08-04, not spelled out in the plan text): Decision 1
# (§6.4) applies the same cap to non-cured Treg responders as to UST/IFX/ADA specifically so
# Treg doesn't get an unearned advantage. A relapse-from-SDR event that happens *after* the cap
# has already fired routes directly to CT rather than back to the (no-longer-available) biologic
# track -- consistent with the cap's own rationale ("exhaustion of a maintenance drug"), and
# because letting only relapsers keep indefinite post-cap biologic-track access would reopen
# exactly the asymmetry Decision 1 exists to close.

if (!exists("MAINTENANCE_STATES")) source("R/utils/transition_matrix.R")
if (!exists("run_maintenance_arm")) source("R/02_markov_engine.R")

#' Convert an assumed annual relapse hazard into a per-cycle probability. `h` is a modelling
#' assumption (still TBD, analysis_plan.md §15 Decision 4 -- no elicitation run yet), not a
#' published cumulative probability, so this is a unit conversion, not a re-derivation -- same
#' constant-hazard formula as R/02_markov_engine.R's discount_factor() reasoning.
hazard_to_cycle_probability <- function(annual_rate, cycle_weeks = 2) {
  stopifnot(annual_rate >= 0, cycle_weeks > 0)
  1 - exp(-annual_rate * cycle_weeks / 52)
}

#' Simulate the Treg arm: UST-equivalent up to the week-56 landmark, then split into a cured
#' (SDR) track and a continuing non-cured track that behaves exactly like run_maintenance_arm()
#' would from there (same switch/cap mechanic, reused via step_maintenance_cycle() rather than
#' duplicated), plus SDR's own decay and cap-aware relapse re-entry.
#'
#' `on_sdr` is tracked as a plain numeric vector (one scalar per cycle), not a 6-column matrix:
#' SDR has no internal state structure -- patients only ever stay or relapse out entirely -- so
#' the matrix shape would just be five permanently-zero columns.
#'
#' If `n_cycles <= landmark_cycle`, the cure never has a chance to manifest within the horizon;
#' this returns the plain UST-equivalent trace (on_sdr all zero), not an error.
run_treg_arm <- function(ust_matrix, ct_matrix, initial_on_biologic, initial_on_ct, n_cycles,
                          pi_sdr, relapse_hazard_annual, landmark_cycle = 28, cap_cycle = 52,
                          apply_cap = TRUE, relapse_destination = "Mild") {
  states <- rownames(ust_matrix)
  stopifnot(
    !is.null(states), identical(states, rownames(ct_matrix)),
    length(initial_on_biologic) == length(states), length(initial_on_ct) == length(states),
    pi_sdr >= 0, pi_sdr <= 1, relapse_hazard_annual >= 0,
    n_cycles >= 0, landmark_cycle > 0, cap_cycle > 0
  )
  ms_idx <- match("Moderate-Severe", states)
  remission_idx <- match("Remission", states)
  relapse_idx <- match(relapse_destination, states)
  stopifnot(!is.na(ms_idx), !is.na(remission_idx), !is.na(relapse_idx))

  phase1_cycles <- min(landmark_cycle, n_cycles)
  phase1 <- run_maintenance_arm(
    ust_matrix, ct_matrix, initial_on_biologic, phase1_cycles,
    cap_cycle = cap_cycle, apply_cap = apply_cap, initial_on_ct = initial_on_ct
  )

  if (n_cycles <= landmark_cycle) {
    return(list(
      on_biologic = phase1$on_biologic, on_ct = phase1$on_ct,
      on_sdr = rep(0, n_cycles + 1), total = phase1$total
    ))
  }

  on_biologic <- matrix(0, nrow = n_cycles + 1, ncol = length(states), dimnames = list(NULL, states))
  on_ct <- matrix(0, nrow = n_cycles + 1, ncol = length(states), dimnames = list(NULL, states))
  on_sdr <- rep(0, n_cycles + 1)
  on_biologic[seq_len(landmark_cycle + 1), ] <- phase1$on_biologic
  on_ct[seq_len(landmark_cycle + 1), ] <- phase1$on_ct

  # Landmark split (§6.2): peel pi_sdr's share of whoever is in Remission under the biologic
  # track at the landmark into SDR; the rest continues unchanged in the ordinary Markov.
  remission_mass <- on_biologic[landmark_cycle + 1, remission_idx]
  on_sdr[landmark_cycle + 1] <- remission_mass * pi_sdr
  on_biologic[landmark_cycle + 1, remission_idx] <- remission_mass * (1 - pi_sdr)

  relapse_prob <- hazard_to_cycle_probability(relapse_hazard_annual)

  for (t in (landmark_cycle + 1):n_cycles) {
    step <- step_maintenance_cycle(
      on_biologic[t, ], on_ct[t, ], ust_matrix, ct_matrix, ms_idx,
      apply_cap_now = apply_cap && t == cap_cycle
    )

    relapsed <- on_sdr[t] * relapse_prob
    on_sdr[t + 1] <- on_sdr[t] - relapsed

    if (!apply_cap || t < cap_cycle) {
      step$on_biologic[relapse_idx] <- step$on_biologic[relapse_idx] + relapsed
    } else {
      # t >= cap_cycle: the biologic track has already been (or is this cycle being) swept to
      # CT for everyone else -- a relapse landing now has no biologic to rejoin either.
      step$on_ct[relapse_idx] <- step$on_ct[relapse_idx] + relapsed
    }

    on_biologic[t + 1, ] <- step$on_biologic
    on_ct[t + 1, ] <- step$on_ct
  }

  list(
    on_biologic = on_biologic, on_ct = on_ct, on_sdr = on_sdr,
    total = rowSums(on_biologic) + rowSums(on_ct) + on_sdr
  )
}
