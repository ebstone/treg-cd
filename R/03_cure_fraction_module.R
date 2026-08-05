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

#' Convert a median duration of Sustained Deep Remission T (years) into the constant annual
#' relapse hazard `hazard_to_cycle_probability()` (and everything downstream of it) actually
#' consumes, via h = ln(2)/T -- the standard median-to-rate conversion for a constant-hazard
#' (exponential) survival model, S(T) = exp(-h*T) = 0.5.
#'
#' Added per `docs/treg-cd_decision_resolutions_2026-08-05.md` §3.3: the briefing this project
#' worked from described the module's prior `h = 0` default as "conservative, no relapse" --
#' backwards. Zero relapse hazard is the single MOST favourable assumption available to the
#' intervention (it maximises time spent in SDR, which is never worse, cost- or utility-wise,
#' than the ongoing non-cured track it replaces -- headroom_pi_star()'s own docstring,
#' R/05_deterministic_results.R), so shipping it silently as "the conservative choice" would have
#' been a correctable error a reviewer would find. T is the quantity actually surfaced to readers
#' and swept across a grid (R/05_deterministic_results.R's RELAPSE_DURATION_GRID_YEARS) --
#' "a median 10 years of drug-free remission" is something a clinician or reviewer can reason
#' about directly; a bare hazard rate of 0.0693/year is not -- while h remains the only quantity
#' any Markov step in this file actually needs.
#'
#' `T = Inf` (permanent remission, h = 0) is an explicit edge case, not something the algebra
#' handles for free: `log(2) / Inf` evaluates to `0` in R's own IEEE arithmetic (not NaN, as the
#' analogous 0/Inf pattern might suggest), so this branch is defensive documentation of that fact,
#' not strictly required for correctness -- kept explicit so the T=Inf case reads as a deliberate
#' modelling point (the Ovasave/CATS1-style "cure never wanes" upper bound), not as a numeric
#' coincidence.
duration_to_hazard <- function(median_duration_years) {
  stopifnot(median_duration_years > 0)
  if (is.infinite(median_duration_years)) return(0)
  log(2) / median_duration_years
}

#' Inverse of duration_to_hazard() -- annual hazard back to its implied median SDR duration in
#' years, for reporting a fitted or swept hazard value back in the units a reader can reason
#' about. h = 0 (permanent remission) maps to T = Inf, not an error or an undefined division.
hazard_to_duration <- function(annual_hazard) {
  stopifnot(annual_hazard >= 0)
  if (annual_hazard == 0) return(Inf)
  log(2) / annual_hazard
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
#'
#' `death_prob_schedule = NULL` (default): fully backward compatible, no mortality applied
#' anywhere in this function -- on_sdr's only exit is relapse, exactly as before this parameter
#' existed. Pass a length-`n_cycles` vector (R/utils/life_table.R) to apply background mortality
#' throughout, phase 1 AND phase 2 alike, INCLUDING to on_sdr itself: SDR patients face the same
#' background mortality as anyone else (R/utils/life_table.R's module header -- "no CD excess
#' mortality" is not "no mortality"), competing each cycle against the relapse hazard for the
#' same pool of mass. This was a silent gap in the pre-lifetime-horizon model (on_sdr had no
#' death exit at all), harmless at a 6.15-year horizon where the compounded miss is negligible,
#' but load-bearing at a lifetime one -- not called directly for that path in practice;
#' run_treg_arm_with_mortality() below is the entry point that actually builds the schedule (per
#' sex) and calls this.
run_treg_arm <- function(ust_matrix, ct_matrix, initial_on_biologic, initial_on_ct, n_cycles,
                          pi_sdr, relapse_hazard_annual, landmark_cycle = 28, cap_cycle = 52,
                          apply_cap = TRUE, relapse_destination = "Mild",
                          death_prob_schedule = NULL) {
  states <- rownames(ust_matrix)
  stopifnot(
    !is.null(states), identical(states, rownames(ct_matrix)),
    length(initial_on_biologic) == length(states), length(initial_on_ct) == length(states),
    pi_sdr >= 0, pi_sdr <= 1, relapse_hazard_annual >= 0,
    n_cycles >= 0, landmark_cycle > 0, cap_cycle > 0,
    is.null(death_prob_schedule) || length(death_prob_schedule) == n_cycles
  )
  ms_idx <- match("Moderate-Severe", states)
  remission_idx <- match("Remission", states)
  relapse_idx <- match(relapse_destination, states)
  stopifnot(!is.na(ms_idx), !is.na(remission_idx), !is.na(relapse_idx))

  phase1_cycles <- min(landmark_cycle, n_cycles)
  phase1_schedule <- if (is.null(death_prob_schedule)) NULL else death_prob_schedule[seq_len(phase1_cycles)]
  phase1 <- run_maintenance_arm(
    ust_matrix, ct_matrix, initial_on_biologic, phase1_cycles,
    cap_cycle = cap_cycle, apply_cap = apply_cap, initial_on_ct = initial_on_ct,
    death_prob_schedule = phase1_schedule
  )

  if (n_cycles <= landmark_cycle) {
    return(list(
      on_biologic = phase1$on_biologic, on_ct = phase1$on_ct,
      on_sdr = rep(0, n_cycles + 1), on_sdr_dead = rep(0, n_cycles + 1), total = phase1$total
    ))
  }

  on_biologic <- matrix(0, nrow = n_cycles + 1, ncol = length(states), dimnames = list(NULL, states))
  on_ct <- matrix(0, nrow = n_cycles + 1, ncol = length(states), dimnames = list(NULL, states))
  on_sdr <- rep(0, n_cycles + 1)
  # Cumulative mass that has died OUT OF the SDR pool specifically (competing-risk "dying" below).
  # SDR has no internal state structure (module header) and so, unlike on_biologic/on_ct, no
  # Death COLUMN of its own for that mass to land in -- tracked here purely so `total` below
  # still conserves cohort mass exactly (this codebase treats that as a tested invariant
  # everywhere else, e.g. test-markov-engine.R's "cohort size is conserved"), not because
  # anything downstream needs the split. attach_treg_costs_utilities() (R/04_costs_utilities.R)
  # never reads `total` -- Death contributes zero cost/utility whether tracked here or not, so
  # this has no effect on any QALY/cost number, only on the accounting check.
  on_sdr_dead <- rep(0, n_cycles + 1)
  on_biologic[seq_len(landmark_cycle + 1), ] <- phase1$on_biologic
  on_ct[seq_len(landmark_cycle + 1), ] <- phase1$on_ct

  # Landmark split (§6.2): peel pi_sdr's share of whoever is in Remission under the biologic
  # track at the landmark into SDR; the rest continues unchanged in the ordinary Markov.
  remission_mass <- on_biologic[landmark_cycle + 1, remission_idx]
  on_sdr[landmark_cycle + 1] <- remission_mass * pi_sdr
  on_biologic[landmark_cycle + 1, remission_idx] <- remission_mass * (1 - pi_sdr)

  relapse_prob <- hazard_to_cycle_probability(relapse_hazard_annual)

  for (t in (landmark_cycle + 1):n_cycles) {
    if (is.null(death_prob_schedule)) {
      bio_m <- ust_matrix
      ct_m <- ct_matrix
      sdr_death_prob <- 0
    } else {
      bio_m <- age_adjust_matrix(ust_matrix, death_prob_schedule[t])
      ct_m <- age_adjust_matrix(ct_matrix, death_prob_schedule[t])
      sdr_death_prob <- death_prob_schedule[t]
    }
    step <- step_maintenance_cycle(
      on_biologic[t, ], on_ct[t, ], bio_m, ct_m, ms_idx,
      apply_cap_now = apply_cap && t == cap_cycle
    )

    # Competing risks on the SDR pool this cycle: death (sdr_death_prob, 0 unless a schedule was
    # supplied -- see this function's own docstring) and relapse (relapse_prob) draw from the
    # same on_sdr[t] mass. Standard competing-risk decomposition: death first, relapse among the
    # survivors -- order doesn't matter to first order at a 2-week cycle, but this is the exact
    # form, not an approximation.
    dying <- on_sdr[t] * sdr_death_prob
    relapsed <- on_sdr[t] * (1 - sdr_death_prob) * relapse_prob
    on_sdr[t + 1] <- on_sdr[t] - dying - relapsed
    on_sdr_dead[t + 1] <- on_sdr_dead[t] + dying

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
    on_biologic = on_biologic, on_ct = on_ct, on_sdr = on_sdr, on_sdr_dead = on_sdr_dead,
    total = rowSums(on_biologic) + rowSums(on_ct) + on_sdr + on_sdr_dead
  )
}

# ---- Lifetime horizon: age- and sex-specific background mortality (R/utils/life_table.R) ------

#' Run run_treg_arm() twice, once per sex, at half the cohort's initial mass each, each against
#' that sex's own death_prob_schedule(), and sum the three occupancy tracks (on_biologic, on_ct,
#' on_sdr). Same rationale and exactness argument as
#' R/02_markov_engine.R's run_maintenance_arm_with_mortality() -- see its header, and
#' R/utils/life_table.R's module header for why summing sub-cohorts is exact here.
run_treg_arm_with_mortality <- function(ust_matrix, ct_matrix, initial_on_biologic, initial_on_ct,
                                         n_cycles, pi_sdr, relapse_hazard_annual, baseline_age,
                                         landmark_cycle = 28, cap_cycle = 52, apply_cap = TRUE,
                                         relapse_destination = "Mild", life_table = NULL,
                                         raw_dir = "data/raw", cycle_weeks = 2) {
  if (is.null(life_table)) life_table <- load_life_table(raw_dir)

  runs <- lapply(c("male", "female"), function(sex) {
    schedule <- death_prob_schedule(baseline_age, sex, n_cycles, cycle_weeks, life_table)
    run_treg_arm(
      ust_matrix, ct_matrix, initial_on_biologic * 0.5, initial_on_ct * 0.5, n_cycles,
      pi_sdr = pi_sdr, relapse_hazard_annual = relapse_hazard_annual,
      landmark_cycle = landmark_cycle, cap_cycle = cap_cycle, apply_cap = apply_cap,
      relapse_destination = relapse_destination, death_prob_schedule = schedule
    )
  })

  list(
    on_biologic = runs[[1]]$on_biologic + runs[[2]]$on_biologic,
    on_ct = runs[[1]]$on_ct + runs[[2]]$on_ct,
    on_sdr = runs[[1]]$on_sdr + runs[[2]]$on_sdr,
    on_sdr_dead = runs[[1]]$on_sdr_dead + runs[[2]]$on_sdr_dead,
    total = runs[[1]]$total + runs[[2]]$total
  )
}
