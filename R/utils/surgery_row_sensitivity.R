# R2 sensitivity check (peer review 2026-08-05; independently re-verified the same day, see
# `treg-cd-peer-review-2026-08-05` project memory and `docs/model_audit_v6.md` A17's own follow-up
# note): Aliyev et al. 2019 reports TWO Surgery-state transition sources that do not reconcile.
# `data/raw/aliyev2019_appendixS1_table2_parameters.csv`'s 8-week "Out of Surgery Transition
# Probability" rows (self-loop 0.34) are dormant data, never wired into this model.
# `data/raw/aliyev2019_appendixS2_table4_maintenance_transition_probabilities.csv`'s 2-week
# Surgery row (self-loop 0.0976, 86.8% one-step return to Remission) IS what every arm in this
# model actually runs on. Compounding one to the other's cycle length doesn't reconcile them:
# Table 4's 2-week self-loop, compounded 4x to 8 weeks, implies a self-loop of roughly 0.0976^4 =
# 9.1e-5 -- nowhere near Table 2's own 0.34. Conversely, treating Table 2's 0.34 as authoritative
# and converting it DOWN to 2 weeks (this file's own `build_table2_surgery_row()`, below) implies
# a 2-week self-loop around 0.76 -- an order of magnitude above Table 4's 0.0976. In plain terms:
# Table 2 describes a Surgery state patients stay in much LONGER than the one actually wired into
# every result this project has ever produced.
#
# Why this matters beyond just "the two tables disagree": S12's own non-cured hazard-ratio scenario
# (`R/03_cure_fraction_module.R`'s `apply_non_cured_hazard_ratio()`, `docs/analysis_plan.md` §6.2)
# found a genuinely counterintuitive result -- strengthening the "advantaged" hazard ratio very
# slightly LOWERS aggregate QALYs, traced to Table 4's Surgery row having an unusually FAST 86.8%
# one-step return to Remission (faster than the ordinary Moderate-Severe -> M-SR -> Mild path), so
# an HR<1 that shrinks entries INTO Surgery removes a QALY-boosting fast lane, not just an
# expensive one (NMB still rises monotonically either way -- the cost saving dominates regardless).
# That finding is only as trustworthy as Table 4's own Surgery dynamics. If Table 2's much STICKIER
# Surgery state is the more credible figure, the "fast lane" that drives the counterintuitive sign
# may not exist, and the sign could flip. This file builds a Table-2-derived, native-2-week-
# equivalent Surgery row as an alternative, explicitly-caveated SENSITIVITY input --
# NOT a replacement for the Table-4-sourced row used everywhere else in this model, which remains
# the base case throughout -- and `R/09_scenarios.R`'s `run_scenario_r2_surgery_sensitivity()`
# re-runs S12's own HR grid against it, reporting whether the sign holds.
#
# ---- Necessary assumptions, stated explicitly (this is a sensitivity check, not a new source) ----
#
# Table 2's own extract reports only THREE of the Surgery row's entries as named PSA parameters:
# "8-Week Surgery to M-S Probability" (0.060), "8-Week Surgery to Mild/M-SR Probability" (0.075 --
# Mild AND M-SR Responder COMBINED, not split), and "8-Week Surgery to Surgery Probability" (0.34,
# the self-loop). These sum to 0.475, not 1 -- the remaining 0.525 of 8-week mass is not reported
# anywhere in this repo's sourced extract (nor, as far as this project's own sourcing has found,
# anywhere else). Three explicit assumptions fill that gap, none of them independently sourced
# beyond the reasoning stated -- if a co-author locates the missing figures in Aliyev's own
# appendix, replace these outright rather than treating them as settled:
#   1. **All unreported residual 8-week mass (0.525) is assigned to Remission.** Table 4's own
#      Surgery row is completely dominated by a Remission destination (86.8%) -- assuming Table 2's
#      own un-extracted residual follows the same "most surgical patients return directly to
#      Remission" pattern is the most defensible reading available given what this repo actually
#      has, not a neutral default (e.g. splitting the residual evenly across every destination
#      would be LESS consistent with the one piece of corroborating evidence on file).
#   2. **The combined "Mild/M-SR" 0.075 figure is split evenly, 50/50.** Table 4's own row splits
#      these two destinations exactly evenly (0.010 / 0.010) -- applied here for consistency with
#      that pattern, not because Table 2 itself supports an even split specifically.
#   3. **Death is assigned 0 probability from Surgery under this reconstruction.** Table 2 reports
#      no Surgery-to-Death figure at any cycle length; Table 4's own Death entry is already
#      negligible (0.00006), so approximating it as exactly 0 here is a small, stated
#      simplification -- inventing a nonzero figure with no source at all would be worse than
#      omitting it.
#
# ---- Cycle-length conversion method ---------------------------------------------------------------
#
# `convert_cycle_length_probability()` (R/utils/transition_matrix.R, added alongside this file) is
# applied to the AGGREGATE "leaves Surgery within 8 weeks" probability (1 - 0.34 = 0.66), not to
# each destination independently -- treating "leaving Surgery" as the single event being
# rate-converted (exactly the method Aliyev's own Appendix S2 worked example uses, B1 fix,
# `docs/model_audit_v6.md` A17), then re-splitting the converted 2-week "leaves" mass across
# destinations in the SAME relative proportions Table 2's own three entries (plus assumption 1's
# Remission residual) imply. This avoids the joint-competing-risks complexity a full matrix
# fractional-power conversion would require, while staying faithful to the one cycle-length
# conversion method this project has already validated against Aliyev's own published numbers.

if (!exists("convert_cycle_length_probability")) source("R/utils/transition_matrix.R")

#' Table 2's own three named Surgery-row entries at their native 8-week cycle (module header
#' above) -- transcribed directly from
#' `data/raw/aliyev2019_appendixS1_table2_parameters.csv`, not re-derived.
TABLE2_SURGERY_8WK <- c(
  to_moderate_severe = 0.060,
  to_mild_msr_combined = 0.075,
  self_loop = 0.34
)

#' Build a Table-2-derived, native-2-week-equivalent Surgery row -- see this file's own header for
#' the full derivation and every assumption involved. Returns a named numeric vector over
#' `MAINTENANCE_STATES`, summing to exactly 1, ready to substitute into the "Surgery" row of any
#' existing therapy matrix (`build_transition_matrix()`'s output).
build_table2_surgery_row <- function() {
  self_loop_8wk <- TABLE2_SURGERY_8WK[["self_loop"]]
  leave_8wk <- 1 - self_loop_8wk
  leave_2wk <- convert_cycle_length_probability(leave_8wk, weeks_old = 8, weeks_new = 2)

  residual_to_remission_8wk <- leave_8wk - TABLE2_SURGERY_8WK[["to_moderate_severe"]] -
    TABLE2_SURGERY_8WK[["to_mild_msr_combined"]]
  stopifnot(residual_to_remission_8wk > 0)  # sanity: Table 2's own three entries must not exceed 1

  # Relative shares of the 8-week "leaves Surgery" mass (assumption 1: residual -> Remission).
  # Sums to 1 by construction (each share is a component of leave_8wk, divided by leave_8wk).
  share_ms <- TABLE2_SURGERY_8WK[["to_moderate_severe"]] / leave_8wk
  share_mild_msr <- TABLE2_SURGERY_8WK[["to_mild_msr_combined"]] / leave_8wk
  share_remission <- residual_to_remission_8wk / leave_8wk

  row <- c(
    "Moderate-Severe" = leave_2wk * share_ms,
    "Moderate-Severe Responder" = leave_2wk * share_mild_msr / 2,  # assumption 2: even split
    "Mild" = leave_2wk * share_mild_msr / 2,                       # assumption 2: even split
    "Remission" = leave_2wk * share_remission,
    "Surgery" = 1 - leave_2wk,
    "Death" = 0                                                    # assumption 3
  )
  stopifnot(isTRUE(all.equal(sum(row), 1, tolerance = 1e-9)))
  row[MAINTENANCE_STATES]  # enforce canonical column order
}

#' Return a copy of therapy matrix `m` with its "Surgery" row replaced by
#' `build_table2_surgery_row()`'s alternative -- every other row (Moderate-Severe,
#' Moderate-Severe Responder, Mild, Remission, Death) is left untouched. Used only for the R2
#' sensitivity check (this file's own header) -- never applied to a real comparator arm's own
#' matrix, only to a local copy `run_scenario_r2_surgery_sensitivity()` (R/09_scenarios.R)
#' constructs for that purpose.
apply_table2_surgery_row <- function(m) {
  stopifnot("Surgery" %in% rownames(m), "Surgery" %in% colnames(m))
  out <- m
  out["Surgery", ] <- build_table2_surgery_row()
  out
}
