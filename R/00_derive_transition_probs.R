# Re-derivation of transition probabilities from published trial endpoints (analysis_plan.md
# §7.3), for the rebuilt engine (R/02_markov_engine.R onward) to consume directly.
#
# SCOPE CHANGE, 2026-08-04: earlier versions of this script tried to reconcile re-derived
# values against the original manuscript / v6 workbook figures, and treated a tight match as
# the Gate 1 pass condition. That target has been dropped: this is a from-scratch rebuild, the
# original figures are known to carry several audited defects (docs/model_audit_v6.md), and in
# at least one case (UST induction) the only candidate method that numerically matched the
# legacy value did so by implicitly using a 2-week window inconsistent with this study's own
# stated 8-week induction design -- i.e. the legacy figure looks like an artifact of Aliyev's
# per-2-week maintenance cycle carried over unchanged, not a value actually re-derived for this
# study's induction window. See DERIVATION_NOTES.md for the full writeup.
#
# This script instead computes ONE documented, defensible set of values per therapy and reports
# the legacy figures alongside only as non-blocking reference context (useful for the
# manuscript's discussion section -- "our re-derived estimates differ from Aliyev's because...").
# Where a genuine data gap remains (no directly observed trial endpoint, or an endpoint
# structured too differently to convert), that is flagged explicitly rather than papered over.

# ---- Core conversions --------------------------------------------------------

#' Convert a cumulative probability observed over `t_from` time units to the
#' probability implied over `t_to` time units, assuming a constant hazard
#' (DEALE: analysis_plan.md §7.3 worked example — UST week-6 remission 0.349
#' -> annual rate 3.720 -> 2-week probability 0.133).
deale_convert <- function(p, t_from, t_to) {
  stopifnot(all(p >= 0 & p < 1), all(t_from > 0), all(t_to > 0))
  rate <- -log(1 - p) / t_from
  1 - exp(-rate * t_to)
}

#' Back out the per-cycle "stay in state" probability implied by a cumulative
#' multi-cycle endpoint (e.g. "remission at week 44 given remission at
#' baseline"), assuming a stationary per-cycle Markov process:
#' p_cumulative = p_per_cycle ^ n_cycles.
retention_root <- function(p_cumulative, n_cycles) {
  stopifnot(all(p_cumulative >= 0 & p_cumulative <= 1), all(n_cycles > 0))
  p_cumulative^(1 / n_cycles)
}

get_param <- function(params, name) {
  hit <- params[params$parameter == name, "base_case_value"]
  if (length(hit) != 1) stop("Expected exactly one match for '", name, "', got ", length(hit))
  hit
}

# ---- Induction ----------------------------------------------------------

#' Re-derive induction destination probabilities for a therapy with a directly
#' observed trial-endpoint response probability and remitter:responder ratio,
#' extrapolated to `target_week`.
#'
#' Method: split the observed response proportion into its two mutually
#' exclusive destination-specific cumulative proportions -- remission
#' (response x ratio) and responder-only (response x (1 - ratio), further
#' split M-SR/Mild via `mild_split`) -- BEFORE DEALE conversion, then convert
#' each independently using its own implied hazard rate.
#'
#' This ordering matters and is the one methodological point worth being
#' explicit about: DEALE assumes a single constant hazard for a single
#' well-defined event. "Remission by week W" and "response-but-not-remission
#' by week W" are two different competing events with two different
#' implied rates -- pooling them into one response-probability hazard and
#' only splitting by the (unconverted, or separately-converted) ratio
#' afterward mixes two timescales and has no clean single-hazard
#' interpretation. Splitting first, converting each piece on its own, does.
derive_induction <- function(therapy, week, p_response, remitter_ratio, mild_split,
                              target_week, confidence = "direct", note = "") {
  stopifnot(p_response > 0, p_response < 1, remitter_ratio >= 0, remitter_ratio <= 1)

  p_remission_obs <- p_response * remitter_ratio
  p_resp_only_obs <- p_response * (1 - remitter_ratio)

  p_remission <- deale_convert(p_remission_obs, week, target_week)
  p_resp_only <- deale_convert(p_resp_only_obs, week, target_week)

  # Converting the two destination-specific proportions independently (see the ordering
  # rationale above) means they do not automatically sum to <= 1 at every horizon: each
  # individually approaches 1 as target_week grows, so their sum approaches 2 and
  # to_moderate_severe approaches -1. It stays valid near the trial window (by construction,
  # it's exactly 1 at target_week == week) but the safe range is finite and therapy-specific --
  # for ADA (week 4 -> target_week 8 is already a 2x extrapolation) it turns negative above
  # target_week ~= 10.2. Fail loudly here rather than let a future scenario change (a longer
  # induction window, a different comparator with a shorter trial week) silently feed a
  # negative transition probability into the Markov engine.
  to_moderate_severe <- 1 - p_remission - p_resp_only
  if (to_moderate_severe < 0 || to_moderate_severe > 1) {
    stop(sprintf(
      paste(
        "derive_induction(%s): extrapolating from week %g to target_week %g produced an",
        "invalid non-response probability (%.4f). The two destination-specific hazards no",
        "longer sum to a valid probability at this horizon -- reduce target_week or revisit",
        "the split-then-convert method for this therapy."
      ),
      therapy, week, target_week, to_moderate_severe
    ))
  }

  data.frame(
    therapy = therapy,
    target_week = target_week,
    confidence = confidence,
    to_moderate_severe = to_moderate_severe,
    to_moderate_severe_responder = p_resp_only * (1 - mild_split),
    to_mild = p_resp_only * mild_split,
    to_remission = p_remission,
    note = note,
    stringsAsFactors = FALSE
  )
}

#' IFX has no directly observed induction trial endpoint in data/raw/ -- only
#' an indirect adjustment of ADA's endpoint via the workbook's "IFX Induction
#' Remission/Responder Adjustment Ratio" parameters. The semantics of those
#' ratios (odds ratio? risk ratio? applied to the raw proportion or to a
#' DEALE-converted probability?) are not documented in this repository's
#' materials (Aliyev Appendix S2 is not present). This is a genuine data gap,
#' not a methodological choice this script can resolve -- flagged low
#' confidence throughout, independent of the target-week question above.
derive_ifx_induction <- function(params, mild_split, target_week) {
  ada_week <- 4
  p_response_ada <- get_param(params, "ADA Week 4 Response Probability")
  ratio_ada <- get_param(params, "ADA Week 4 Remitter:Responder Ratio")
  responder_adj <- get_param(params, "IFX Induction Responder Adjustment Ratio")
  remission_adj <- get_param(params, "IFX Induction Remission Adjustment Ratio")
  efficacy_ratio <- get_param(params, "IFX:ADA Efficacy Ratio")

  p_response_ifx <- p_response_ada * responder_adj * efficacy_ratio
  ratio_ifx <- min(ratio_ada * remission_adj, 1)

  derive_induction(
    "IFX",
    week = ada_week, p_response = p_response_ifx, remitter_ratio = ratio_ifx,
    mild_split = mild_split, target_week = target_week,
    confidence = "indirect_low_confidence",
    note = paste(
      "Derived from ADA's week-4 endpoint via the IFX Induction Responder/Remission",
      "Adjustment Ratios and the IFX:ADA Efficacy Ratio, applied as direct",
      "multiplicative adjustments to the raw ADA proportions before DEALE conversion,",
      "and anchored to ADA's week-4 timepoint rather than IFX's own ~6-week trial",
      "window (no directly observed IFX endpoint is available to anchor to instead).",
      "The adjustment-ratio semantics are undocumented in data/raw/ -- this is an open",
      "item for co-author review, not a validated re-derivation."
    )
  )
}

# ---- Maintenance ----------------------------------------------------------

#' Re-derive the one maintenance endpoint structure that maps cleanly onto a
#' single matrix cell: "remission probability at week W, given remission at
#' baseline" (UST week 44, CT/placebo week 44), via retention_root(). IFX and
#' ADA maintenance endpoints in data/raw/ are structured differently (an
#' overall week-52 response probability + remitter ratio, not a
#' baseline-conditional remission probability) and are not derived here --
#' see DERIVATION_NOTES.md.
derive_maintenance <- function(params, legacy_maintenance_df, cycle_weeks = 8) {
  specs <- list(
    UST = list(param = "UST Week 44 Remission Probability in Patients in Remission at Baseline", week = 44),
    CT  = list(param = "CT Week 44 Remission Probability in Patients in Remission at Baseline", week = 44)
  )
  rows <- lapply(names(specs), function(therapy) {
    spec <- specs[[therapy]]
    p_cumulative <- get_param(params, spec$param)
    n_cycles <- spec$week / cycle_weeks
    p_derived <- retention_root(p_cumulative, n_cycles)
    p_legacy <- legacy_maintenance_df[
      legacy_maintenance_df$therapy == therapy & legacy_maintenance_df$scenario == "base" &
        legacy_maintenance_df$from_state == "Remission" & legacy_maintenance_df$to_state == "Remission",
      "probability_per_8wk_cycle"
    ]
    data.frame(
      therapy = therapy, from_state = "Remission", to_state = "Remission",
      trial_endpoint = spec$param, trial_endpoint_value = p_cumulative, trial_endpoint_week = spec$week,
      derived_probability_per_8wk_cycle = p_derived,
      legacy_probability_per_8wk_cycle = p_legacy,
      note = "Legacy value retained for reference only, not a validation target (2026-08-04).",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# ---- Orchestration -------------------------------------------------------

run_derivation <- function(raw_dir = "data/raw", proc_dir = "data/processed",
                            induction_target_week = 8, write_output = TRUE) {
  params <- utils::read.csv(file.path(raw_dir, "aliyev2019_appendixS1_table2_parameters.csv"), stringsAsFactors = FALSE)
  legacy_induction <- utils::read.csv(file.path(raw_dir, "manuscript_supplement1_induction_transition_probabilities.csv"), stringsAsFactors = FALSE)
  legacy_maintenance <- utils::read.csv(file.path(proc_dir, "model_maintenance_transition_probabilities.csv"), stringsAsFactors = FALSE)

  mild_split <- get_param(params, "Proportion of Responders that Transition to Mild")

  ust <- derive_induction(
    "UST", week = 6,
    p_response = get_param(params, "UST Week 6 Response Probability"),
    remitter_ratio = get_param(params, "UST Week 6 Remitter:Responder Ratio"),
    mild_split = mild_split, target_week = induction_target_week
  )
  ada <- derive_induction(
    "ADA", week = 4,
    p_response = get_param(params, "ADA Week 4 Response Probability"),
    remitter_ratio = get_param(params, "ADA Week 4 Remitter:Responder Ratio"),
    mild_split = mild_split, target_week = induction_target_week
  )
  ifx <- derive_ifx_induction(params, mild_split, induction_target_week)

  induction <- rbind(ust, ada, ifx)
  induction$row_sum_check <- with(
    induction,
    to_moderate_severe + to_moderate_severe_responder + to_mild + to_remission
  )

  # Legacy manuscript values, carried through as reference context only -- NOT a
  # validation target (2026-08-04 scope change). Kept so the two can be diffed for
  # the manuscript's discussion section.
  therapy_map <- c(Ustekinumab = "UST", Infliximab = "IFX", "T-regulatory" = "TREG")
  legacy_tagged <- data.frame(
    therapy = therapy_map[legacy_induction$treatment],
    target_week = NA_real_, confidence = "legacy_reference_only",
    to_moderate_severe = legacy_induction$to_moderate_severe,
    to_moderate_severe_responder = legacy_induction$to_moderate_severe_responder,
    to_mild = legacy_induction$to_mild,
    to_remission = legacy_induction$to_remission,
    note = "Original manuscript/v6-workbook value (Aliyev Supplementary Table 3 for UST/IFX). Reference only.",
    row_sum_check = legacy_induction$to_moderate_severe + legacy_induction$to_moderate_severe_responder +
      legacy_induction$to_mild + legacy_induction$to_remission,
    stringsAsFactors = FALSE
  )
  legacy_tagged <- legacy_tagged[!is.na(legacy_tagged$therapy), ]

  induction_out <- rbind(induction, legacy_tagged)
  induction_out <- induction_out[order(induction_out$therapy, induction_out$confidence), ]

  maintenance_out <- derive_maintenance(params, legacy_maintenance, cycle_weeks = 8)

  if (write_output) {
    dir.create(proc_dir, showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(induction_out, file.path(proc_dir, "derived_induction_transition_probabilities.csv"), row.names = FALSE)
    utils::write.csv(maintenance_out, file.path(proc_dir, "derived_maintenance_probabilities.csv"), row.names = FALSE)
    write_derivation_notes(induction_out, maintenance_out, proc_dir, induction_target_week)
  }

  list(induction = induction_out, maintenance = maintenance_out)
}

write_derivation_notes <- function(induction_out, maintenance_out, proc_dir, induction_target_week) {
  derived <- induction_out[induction_out$confidence != "legacy_reference_only", ]
  legacy <- induction_out[induction_out$confidence == "legacy_reference_only", ]

  compare_lines <- vapply(intersect(derived$therapy, legacy$therapy), function(th) {
    a <- derived[derived$therapy == th, ]
    r <- legacy[legacy$therapy == th, ]
    sprintf(
      "- **%s**: re-derived to_remission %.4f vs. legacy %.4f; to_moderate_severe %.4f vs. legacy %.4f (legacy is reference only, not a target)",
      th, a$to_remission, r$to_remission, a$to_moderate_severe, r$to_moderate_severe
    )
  }, character(1))

  maint_lines <- sprintf(
    "- **%s** Remission->Remission (8-week cycle): re-derived %.4f vs. legacy %.4f, from '%s' (%d-week endpoint = %.3f)",
    maintenance_out$therapy, maintenance_out$derived_probability_per_8wk_cycle,
    maintenance_out$legacy_probability_per_8wk_cycle, maintenance_out$trial_endpoint,
    maintenance_out$trial_endpoint_week, maintenance_out$trial_endpoint_value
  )

  notes <- c(
    "# Gate 1 transition-probability derivation",
    "",
    sprintf("Generated by `R/00_derive_transition_probs.R` on %s.", format(Sys.Date())),
    "",
    "**2026-08-04 scope change:** this rebuild no longer targets reconciliation with the",
    "original manuscript/v6-workbook transition probabilities as a Gate 1 pass condition.",
    "Those values are carried in the output CSVs as `legacy_reference_only` rows, useful for",
    "the manuscript's discussion section, but not a check this script tries to pass. The rows",
    "tagged `direct` / `indirect_low_confidence` below are what the rebuilt engine",
    "(`R/02_markov_engine.R` onward) actually uses.",
    "",
    "## Method",
    "",
    sprintf(paste(
      "**Induction** (target window: %d weeks, matching this study's stated 8-week",
      "induction alignment across biologics -- configurable via `induction_target_week`):"
    ), induction_target_week),
    "each therapy's trial-endpoint response proportion is split into two mutually exclusive",
    "destination-specific cumulative proportions -- remission (response x remitter:responder",
    "ratio) and responder-only (response x (1 - ratio), further split M-SR/Mild via the",
    "workbook's 50/50 mild-split assumption) -- *before* DEALE conversion, and each is then",
    "converted independently using its own implied hazard rate. An earlier version of this",
    "script converted the pooled response probability first and split by the ratio afterward;",
    "that ordering mixes two different events' timescales into one hazard and does not",
    "reproduce a defensible probability. The values here are extrapolated forward from each",
    "trial's own (shorter) observation window under DEALE's constant-hazard assumption --",
    "worth a plausibility check against clinical judgment, since that assumption is more",
    "tested near the observed window than several weeks beyond it.",
    "",
    "**Maintenance:** only the Remission -> Remission cell has a directly comparable trial",
    "endpoint (\"remission probability at week 44, given remission at baseline\"); derived via a",
    "per-cycle retention-root conversion, not DEALE, since it's a multi-cycle cumulative",
    "retention rather than a single-transition hazard. Unaffected by the induction scope change.",
    "",
    "## Re-derived values vs. legacy (context only, not a check)",
    "",
    "### Induction",
    paste(compare_lines, collapse = "\n"),
    "",
    "### Maintenance",
    paste(maint_lines, collapse = "\n"),
    "",
    "## Open items requiring co-author input before Gate 2",
    "",
    "1. **IFX has no directly observed induction trial endpoint** in `data/raw/`; derived",
    "   indirectly from ADA's endpoint via adjustment ratios whose semantics (odds ratio vs.",
    "   risk ratio; point of application relative to DEALE conversion; anchored to ADA's",
    "   week-4 timepoint rather than IFX's own ~6-week trial window) are not documented here.",
    "   Flagged `indirect_low_confidence` -- needs either Aliyev Appendix S2 or an explicit",
    "   documented assumption. This is a data gap, not a target-week or ordering question.",
    "2. **ADA and IFX maintenance trial endpoints** are structured as an overall week-52",
    "   response probability + remitter ratio rather than a baseline-conditional remission",
    "   probability, so they are not yet derivable by `derive_maintenance()`. Needed before the",
    "   full maintenance matrix (not just the Remission cell) can be independently derived for",
    "   all four therapies.",
    "",
    "## Gate 1 status",
    "",
    "**Open**, narrowed to the two items above. UST and ADA induction, and the UST/CT",
    "maintenance Remission cell, are now derived on a single documented method with no",
    "remaining structural ambiguity -- they no longer need co-author disambiguation between",
    "candidate methods, just a plausibility read given the magnitude shift from the legacy",
    "figures shown above. IFX induction and the IFX/ADA maintenance matrix remain genuine data",
    "gaps requiring either Appendix S2 or an explicit documented assumption before Gate 2."
  )
  writeLines(notes, file.path(proc_dir, "DERIVATION_NOTES.md"))
}

if (sys.nframe() == 0L) {
  result <- run_derivation()
  cat("Wrote data/processed/derived_induction_transition_probabilities.csv,",
      "data/processed/derived_maintenance_probabilities.csv,",
      "data/processed/DERIVATION_NOTES.md\n")
}
