# Refractory-population multipliers for scenario S3 (analysis_plan.md §10.3, Decision 3) --
# the resolutions memo's own priority-order item 6 (docs/treg-cd_decision_resolutions_2026-08-05.md
# §11), closed 2026-08-05.
#
# ---- What this is, and what it deliberately is NOT --------------------------------------------
#
# analysis_plan.md §5 (the paragraph introducing Decision 3) is explicit about method: "do not
# attempt to re-derive refractory transition probabilities from first principles. Instead apply
# an explicit, sourced set of multipliers to the biologic-naive probabilities (reduced response
# and remission probabilities, elevated surgery hazard) drawn from published second-line and
# post-anti-TNF-failure CD cohorts, and vary those multipliers widely. State the approach as a
# scenario built on assumption, not as a separately estimated model." This file is exactly that:
# a multiplicative adjustment to Aliyev's own published (biologic-naive) UST/IFX/ADA matrices,
# not a from-scratch refractory model. CT's own matrix is left untouched: CT already functions as
# this model's "failed biologic" fallback track, and no refractory-specific CT data was sought.
#
# ---- Surgery hazard: sourced, but a genuinely imperfect proxy -- read before using -------------
#
# Response/remission (UNITI-1 vs UNITI-2, below) is an unusually clean pairing: same drug, same
# dose, same endpoint, same trial programme, differing only in prior-anti-TNF-failure eligibility.
# No equivalently clean pairing exists for surgery. Trials report response/remission because
# that's what they're powered for; real-world bio-naive-vs-bio-experienced cohorts (checked:
# Parra et al. 2022 BMC Gastroenterol, the ICC registry, several others) report clinical/
# endoscopic remission by that split but not surgery. What IS available
# (data/raw/refractory_surgery_hazard_multiplier.csv) is two SEPARATE studies at two DIFFERENT
# severity tiers: SOJOURN's 11.6% one-year cumulative surgery incidence in biologic-naive
# ustekinumab-treated patients (the right severity tier), and Kassouri et al. 2020's 23.5%
# surgery rate at week 48 in patients who had already failed anti-TNF AND a second-line biologic
# (vedolizumab or ustekinumab) -- i.e. starting a THIRD line, a MORE severe population than this
# project's refractory scenario (which is UNITI-1-equivalent: single anti-TNF failure, starting a
# SECOND line). Ratio-ing these two gives a multiplier (2.03x) that is a deliberate, documented
# upper-bound/conservative proxy for surgery-hazard elevation, not a matched point estimate at the
# same severity tier as the response/remission multipliers -- flagged explicitly at every point
# this multiplier surfaces (loaded value, applied matrices, downstream tables/docs), per Eric's
# own instruction (2026-08-05) to use it with that caveat stated plainly rather than search
# indefinitely for a better-matched source that may not exist in the published literature.
#
# ---- Source: UNITI-1 vs UNITI-2 (Feagan et al. 2016, NEJM) -------------------------------------
#
# UNITI-1 and UNITI-2 are ustekinumab's own two induction trials, run identically (same drug,
# same doses, same endpoints, same 8-week window, same investigators) and differing ONLY in
# which patients were eligible: UNITI-1 enrolled patients with primary/secondary nonresponse to
# >=1 prior anti-TNF agent or unacceptable side effects (this project's "refractory" population);
# UNITI-2 enrolled patients in whom conventional (immunosuppressant/glucocorticoid) therapy
# failed, 68.6% of whom were anti-TNF-naive -- much closer to (not identical to) this project's
# biologic-naive base-case population. Comparing UNITI-1 to UNITI-2, same drug and dose, isolates
# the refractory-population effect about as cleanly as published data allows, which is exactly
# why this project uses this pair rather than trying to synthesise a multiplier from unrelated
# trials/drugs. Full citation and exact figures: data/raw/uniti_refractory_population_multipliers.csv.
#
# Two multipliers, both ustekinumab-specific but applied to IFX/ADA's matrices too (a stated
# cross-drug proxy assumption -- IFX/ADA don't have an equivalently clean paired refractory/naive
# trial design available; this is the "vary...widely" scenario the plan calls for, not a claim
# that IFX/ADA's own refractory effect is identical to UST's):
#
# - **Induction** (week 8, 6 mg/kg IV dose, the FDA-labelled induction dose Aliyev's induction
#   matrix is keyed to): response multiplier = 37.8% / 57.9% (UNITI-1/UNITI-2); remission
#   multiplier = 20.9% / 40.2%.
# - **Maintenance** (IM-UNITI week 44, 90 mg SC every-8-weeks dosing -- matches this project's own
#   UST maintenance cadence, R/04_costs_utilities.R): remission multiplier = 41.1% / 62.5%,
#   **cumulative over the week-8-to-week-44 (36-week / 18 native-2-week-cycle) follow-up window**,
#   not a per-cycle figure -- converted to a per-cycle multiplier below before use, since applying
#   a 36-week cumulative ratio as if it were a single-cycle multiplier would compound far too
#   aggressively over a lifetime horizon's hundreds of cycles (the same reasoning
#   R/utils/life_table.R's age_adjust_matrix() already applies per-cycle, not once).
#
# ---- Why response and remission get separate multipliers, and how "worse" is defined ----------
#
# Response and remission are nested (every patient in remission is also, definitionally, a
# responder) but the two multipliers differ (0.5 vs 0.65-ish) -- refractory disease depresses the
# CHANCE OF REMISSION specifically more than it depresses the chance of any response at all, which
# is itself a real, sourced finding, not noise, and is preserved rather than averaged away.
# apply_refractory_multiplier_induction() below shrinks total-response and remission separately
# then rebuilds the response-but-not-remission share as the difference, floor-clamped so remission
# never exceeds total response (a probability-mass constraint, not a modelling choice).
# apply_refractory_multiplier_maintenance() shrinks each cycle's flow into Remission and
# redistributes the removed mass into Mild specifically (the state immediately below Remission on
# the severity axis in every row this project's matrices have a nonzero Remission entry for,
# including the Surgery row) -- a deliberately narrow, sourced-quantity-only adjustment; this
# project does NOT also reshuffle Moderate-Severe-Responder<->Mild flows, since no maintenance-
# phase multiplier for those was sourced this pass (only Remission's rate was reported by
# refractory-subgroup in Figure 3B). The same function's optional surgery-hazard step (below) runs
# AFTER the Remission/Mild adjustment, in the same row, using age_adjust_matrix()'s own
# proportional-rescale mechanic (R/utils/transition_matrix.R) to install a higher Surgery
# probability: every OTHER column in the row (including the just-adjusted Remission/Mild, and
# Death) is shrunk by the same factor, `(1 - new_surgery) / (1 - old_surgery)`, so the row still
# sums to 1 with no separately-chosen "which column pays for it" rule -- exactly the same
# mechanic this codebase already uses and tests for installing a new Death probability.

if (!exists("build_transition_matrix")) source("R/utils/transition_matrix.R")

#' Load the sourced UNITI-1/UNITI-2 (+ IM-UNITI) rates and derive both multipliers. Returns a
#' list, not just a data.frame, so callers get ready-to-use scalars without re-deriving the
#' induction/response/remission arithmetic at every call site.
load_refractory_multipliers <- function(raw_dir = "data/raw") {
  df <- utils::read.csv(
    file.path(raw_dir, "uniti_refractory_population_multipliers.csv"),
    stringsAsFactors = FALSE
  )

  rate <- function(phase, population, endpoint, dose = NULL) {
    sel <- df$trial_phase == phase & df$population == population & df$endpoint == endpoint
    if (!is.null(dose)) sel <- sel & df$dose == dose
    stopifnot(sum(sel) == 1)
    df$rate_pct[sel] / 100
  }

  induction_response_uniti1 <- rate("induction", "UNITI-1", "clinical response (CDAI decrease >=100 or CDAI <150)")
  induction_response_uniti2 <- rate("induction", "UNITI-2", "clinical response (CDAI decrease >=100 or CDAI <150)")
  induction_remission_uniti1 <- rate("induction", "UNITI-1", "clinical remission (CDAI <150)")
  induction_remission_uniti2 <- rate("induction", "UNITI-2", "clinical remission (CDAI <150)")

  maintenance_remission_uniti1 <- rate("maintenance", "UNITI-1 population (IM-UNITI subgroup)",
                                        "clinical remission (CDAI <150)", dose = "90 mg SC q8w")
  maintenance_remission_uniti2 <- rate("maintenance", "UNITI-2 population (IM-UNITI subgroup)",
                                        "clinical remission (CDAI <150)", dose = "90 mg SC q8w")

  surgery_df <- utils::read.csv(
    file.path(raw_dir, "refractory_surgery_hazard_multiplier.csv"),
    stringsAsFactors = FALSE
  )
  surgery_naive <- surgery_df$rate_pct[surgery_df$population == "biologic_naive"] / 100
  surgery_refractory <- surgery_df$rate_pct[surgery_df$population == "refractory_double_biologic_failure"] / 100
  stopifnot(length(surgery_naive) == 1, length(surgery_refractory) == 1)

  list(
    raw = df,
    surgery_raw = surgery_df,
    induction_response_multiplier = induction_response_uniti1 / induction_response_uniti2,
    induction_remission_multiplier = induction_remission_uniti1 / induction_remission_uniti2,
    maintenance_remission_multiplier_cumulative = maintenance_remission_uniti1 / maintenance_remission_uniti2,
    maintenance_cumulative_n_cycles = 18,  # week 8 -> week 44 = 36 weeks = 18 native 2-week cycles
    # UPPER-BOUND PROXY, not a matched point estimate -- see this file's own module header
    # ("Surgery hazard: sourced, but a genuinely imperfect proxy") before using this value.
    surgery_hazard_multiplier_cumulative = surgery_refractory / surgery_naive,
    surgery_cumulative_n_cycles = 26  # ~1 year (52 weeks) = 26 native 2-week cycles
  )
}

#' Apply the induction multipliers to one therapy's induction row (the
#' to_moderate_severe/to_moderate_severe_responder/to_mild/to_remission partition
#' `load_published_induction()` returns, R/00_derive_transition_probs.R). Returns a row with the
#' same four columns, still summing to 1 (up to floating point).
#'
#' `remission` shrinks by `remission_multiplier` directly. `total_response` (mild + responder +
#' remission) shrinks by `response_multiplier`, then is floor-clamped at the new remission value
#' (total response cannot be less than remission itself -- remission is a subset of response).
#' The response-but-not-remission remainder is split between mild/moderate-severe-responder in
#' the SAME proportion the naive row used (both are 0.0377/0.0377, i.e. exactly 50/50, for UST;
#' kept general here rather than hardcoding that split for IFX/ADA, whose naive rows are also
#' 50/50 but shouldn't be assumed so implicitly).
apply_refractory_multiplier_induction <- function(induction_row, response_multiplier, remission_multiplier) {
  stopifnot(nrow(induction_row) == 1)
  mild <- induction_row$to_mild
  msr <- induction_row$to_moderate_severe_responder
  remission <- induction_row$to_remission
  total_response <- mild + msr + remission

  new_remission <- remission * remission_multiplier
  new_total_response <- max(total_response * response_multiplier, new_remission)
  response_only_new <- new_total_response - new_remission

  mild_frac <- if ((mild + msr) > 0) mild / (mild + msr) else 0.5
  new_mild <- response_only_new * mild_frac
  new_msr <- response_only_new * (1 - mild_frac)

  out <- induction_row
  out$to_remission <- new_remission
  out$to_mild <- new_mild
  out$to_moderate_severe_responder <- new_msr
  out$to_moderate_severe <- 1 - new_total_response
  out
}

#' Apply the maintenance remission multiplier to one square maintenance transition matrix
#' (build_transition_matrix()'s output, states x states with "Remission" and "Mild" among them),
#' per cycle, for `n_cycles_horizon` cycles of exposure to the refractory effect. Converts the
#' sourced CUMULATIVE (36-week) multiplier to a per-cycle one via
#' `cumulative_multiplier ^ (1 / cumulative_n_cycles)` -- the same "one multiplicative step per
#' Markov cycle" design as age_adjust_matrix() (R/utils/transition_matrix.R), so a lifetime
#' horizon's hundreds of cycles compound this correctly instead of applying a 36-week ratio as if
#' it held over 40+ years.
#'
#' Every row with a nonzero Remission entry has that entry scaled down by the per-cycle
#' multiplier; the removed mass moves to Mild (the state immediately below Remission on the
#' severity axis in every row this project's matrices carry a Remission entry for, including the
#' Surgery row's 0.868 -> Remission). Rows without a Remission entry (e.g. Death's absorbing
#' self-loop) are untouched. Row sums are preserved exactly (mass moved column-to-column within
#' the same row, nothing added or removed).
#'
#' `surgery_hazard_multiplier_cumulative`/`surgery_cumulative_n_cycles = NULL` (default): no
#' surgery adjustment, byte-identical to this function's pre-2026-08-05 behaviour. Passing both
#' additionally elevates each row's Surgery-column probability, run AFTER the Remission/Mild step
#' above, using the SAME cumulative-to-per-cycle conversion and age_adjust_matrix()'s own
#' proportional-rescale mechanic (R/utils/transition_matrix.R) -- every other column in the row,
#' including the just-adjusted Remission/Mild and Death, is shrunk by
#' `(1 - new_surgery) / (1 - old_surgery)` so the row still sums to 1. Rows with a zero baseline
#' Surgery probability are left at zero (0 * multiplier = 0), matching the Remission step's own
#' "nothing to adjust" convention. **This multiplier is a deliberate upper-bound/conservable
#' proxy, not a matched point estimate** -- see this file's own module header
#' ("Surgery hazard: sourced, but a genuinely imperfect proxy") before passing a non-NULL value.
apply_refractory_multiplier_maintenance <- function(m, cumulative_multiplier, cumulative_n_cycles,
                                                     surgery_hazard_multiplier_cumulative = NULL,
                                                     surgery_cumulative_n_cycles = NULL) {
  stopifnot("Remission" %in% colnames(m), "Mild" %in% colnames(m))
  per_cycle_multiplier <- cumulative_multiplier ^ (1 / cumulative_n_cycles)
  out <- m
  for (i in seq_len(nrow(m))) {
    old_rem <- m[i, "Remission"]
    if (old_rem <= 0) next
    new_rem <- old_rem * per_cycle_multiplier
    out[i, "Remission"] <- new_rem
    out[i, "Mild"] <- out[i, "Mild"] + (old_rem - new_rem)
  }

  if (!is.null(surgery_hazard_multiplier_cumulative)) {
    stopifnot(!is.null(surgery_cumulative_n_cycles), "Surgery" %in% colnames(m))
    per_cycle_surgery_multiplier <- surgery_hazard_multiplier_cumulative ^ (1 / surgery_cumulative_n_cycles)
    for (i in seq_len(nrow(out))) {
      old_surgery <- out[i, "Surgery"]
      if (old_surgery <= 0 || old_surgery >= 1) next
      new_surgery <- min(old_surgery * per_cycle_surgery_multiplier, 1)
      scale <- (1 - new_surgery) / (1 - old_surgery)
      out[i, ] <- out[i, ] * scale
      out[i, "Surgery"] <- new_surgery
    }
  }

  out
}
