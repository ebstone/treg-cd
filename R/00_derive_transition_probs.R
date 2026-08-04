# DEALE re-derivation of transition probabilities from published trial endpoints (analysis_plan.md §7.3).
# Run first — Gate 1 (§14) requires re-derived probabilities to reconcile with the v6 workbook, or
# discrepancies to be explained in writing. Output: data/processed/ CSV, versioned.
#
# This does not assume Aliyev's Appendix S2 (not present in this repository) is fully
# reconstructable. Where the exact apportionment method is not recoverable from the materials in
# `data/raw/`, this script computes multiple documented candidate conversions and reports them
# side by side, rather than silently picking one that happens to match. See
# `data/processed/DERIVATION_NOTES.md` for the write-up Gate 1 requires.

# ---- DEALE conversion -------------------------------------------------------

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

# ---- Induction candidates ----------------------------------------------------

#' Split a DEALE-converted response probability into the model's four
#' induction destination states, using the workbook's "Proportion of
#' Responders that Transition to Mild" assumption to divide non-remitter
#' responders between Mild and Mod/Sev Responder (M-SR).
split_induction_destinations <- function(p_response_cyc, ratio_cyc, mild_split) {
  p_remission <- p_response_cyc * ratio_cyc
  p_resp_only <- p_response_cyc * (1 - ratio_cyc)
  data.frame(
    to_moderate_severe = 1 - p_response_cyc,
    to_moderate_severe_responder = p_resp_only * (1 - mild_split),
    to_mild = p_resp_only * mild_split,
    to_remission = p_remission
  )
}

#' Build candidate induction transition probabilities for a therapy with a
#' directly observed trial-endpoint response probability and remitter:responder
#' ratio (UST, ADA). Two methods are computed because it is not documented in
#' the materials available here whether the remitter:responder ratio itself
#' should be DEALE-converted over the trial window (nested_hazard) or applied
#' unconverted at the target cycle length (response_hazard_only).
build_direct_induction_candidates <- function(therapy, week, p_response, remitter_ratio,
                                               mild_split, cycle_weeks_grid = c(2, 8)) {
  do.call(rbind, lapply(cycle_weeks_grid, function(cw) {
    do.call(rbind, lapply(c("response_hazard_only", "nested_hazard"), function(method) {
      p_response_cyc <- deale_convert(p_response, week, cw)
      ratio_cyc <- if (method == "nested_hazard") deale_convert(remitter_ratio, week, cw) else remitter_ratio
      dest <- split_induction_destinations(p_response_cyc, ratio_cyc, mild_split)
      cbind(
        therapy = therapy, method = method, cycle_weeks = cw, confidence = "direct",
        dest,
        note = ""
      )
    }))
  }))
}

#' Best-effort candidate for IFX, which has no directly observed induction
#' trial endpoint in `data/raw/` — only an indirect adjustment of ADA's
#' endpoint via "IFX Induction Remission/Responder Adjustment Ratio". The
#' semantics of those adjustment ratios (odds ratio? risk ratio? applied to
#' the raw proportion or to a DEALE-converted probability?) are not documented
#' in this repository's materials (Aliyev Appendix S2 is not present). This
#' candidate treats them as direct multiplicative adjustments to the raw ADA
#' proportions before DEALE conversion — flagged low confidence throughout.
build_ifx_indirect_induction_candidate <- function(params, mild_split, cycle_weeks_grid = c(2, 8)) {
  ada_week <- 4
  p_response_ada <- get_param(params, "ADA Week 4 Response Probability")
  ratio_ada <- get_param(params, "ADA Week 4 Remitter:Responder Ratio")
  responder_adj <- get_param(params, "IFX Induction Responder Adjustment Ratio")
  remission_adj <- get_param(params, "IFX Induction Remission Adjustment Ratio")
  efficacy_ratio <- get_param(params, "IFX:ADA Efficacy Ratio")

  p_response_ifx <- p_response_ada * responder_adj * efficacy_ratio
  ratio_ifx <- min(ratio_ada * remission_adj, 1)

  do.call(rbind, lapply(cycle_weeks_grid, function(cw) {
    p_response_cyc <- deale_convert(p_response_ifx, ada_week, cw)
    dest <- split_induction_destinations(p_response_cyc, ratio_ifx, mild_split)
    cbind(
      therapy = "IFX", method = "indirect_adjustment", cycle_weeks = cw, confidence = "indirect_low_confidence",
      dest,
      note = "Derived from ADA endpoint via IFX Induction Responder/Remission Adjustment Ratios and IFX:ADA Efficacy Ratio; assumed applied at ADA's week-4 timepoint. Adjustment-ratio semantics are not documented in data/raw/ (Aliyev Appendix S2 not present in this repository) — treat as a starting point for co-author review, not a validated re-derivation."
    )
  }))
}

# ---- Maintenance reconciliation ----------------------------------------------

#' Reconcile the one maintenance endpoint structure that maps cleanly onto a
#' single matrix cell: "remission probability at week W, given remission at
#' baseline" (UST week 44, CT/placebo week 44) against the model's
#' Remission -> Remission per-8-week-cycle probability, via retention_root().
#' IFX and ADA maintenance endpoints in data/raw/ are structured differently
#' (an overall week-52 response probability + remitter ratio, not a
#' baseline-conditional remission probability) and are not reconciled here —
#' see DERIVATION_NOTES.md.
build_maintenance_reconciliation <- function(params, maintenance_df, cycle_weeks = 8) {
  specs <- list(
    UST = list(param = "UST Week 44 Remission Probability in Patients in Remission at Baseline", week = 44),
    CT  = list(param = "CT Week 44 Remission Probability in Patients in Remission at Baseline", week = 44)
  )
  rows <- lapply(names(specs), function(therapy) {
    spec <- specs[[therapy]]
    p_cumulative <- get_param(params, spec$param)
    n_cycles <- spec$week / cycle_weeks
    p_derived <- retention_root(p_cumulative, n_cycles)
    p_model <- maintenance_df[
      maintenance_df$therapy == therapy & maintenance_df$scenario == "base" &
        maintenance_df$from_state == "Remission" & maintenance_df$to_state == "Remission",
      "probability_per_8wk_cycle"
    ]
    data.frame(
      therapy = therapy, from_state = "Remission", to_state = "Remission",
      trial_endpoint = spec$param, trial_endpoint_value = p_cumulative, trial_endpoint_week = spec$week,
      derived_probability_per_8wk_cycle = p_derived,
      model_probability_per_8wk_cycle = p_model,
      abs_diff = abs(p_derived - p_model)
    )
  })
  do.call(rbind, rows)
}

# ---- Orchestration -------------------------------------------------------

run_deale_rederivation <- function(raw_dir = "data/raw", proc_dir = "data/processed", write_output = TRUE) {
  params <- utils::read.csv(file.path(raw_dir, "aliyev2019_appendixS1_table2_parameters.csv"), stringsAsFactors = FALSE)
  current_induction <- utils::read.csv(file.path(raw_dir, "manuscript_supplement1_induction_transition_probabilities.csv"), stringsAsFactors = FALSE)
  maintenance_df <- utils::read.csv(file.path(proc_dir, "model_maintenance_transition_probabilities.csv"), stringsAsFactors = FALSE)

  mild_split <- get_param(params, "Proportion of Responders that Transition to Mild")

  ust_candidates <- build_direct_induction_candidates(
    "UST", week = 6,
    p_response = get_param(params, "UST Week 6 Response Probability"),
    remitter_ratio = get_param(params, "UST Week 6 Remitter:Responder Ratio"),
    mild_split = mild_split
  )
  ada_candidates <- build_direct_induction_candidates(
    "ADA", week = 4,
    p_response = get_param(params, "ADA Week 4 Response Probability"),
    remitter_ratio = get_param(params, "ADA Week 4 Remitter:Responder Ratio"),
    mild_split = mild_split
  )
  ifx_candidates <- build_ifx_indirect_induction_candidate(params, mild_split)

  induction_candidates <- rbind(ust_candidates, ada_candidates, ifx_candidates)
  induction_candidates$row_sum_check <- with(
    induction_candidates,
    to_moderate_severe + to_moderate_severe_responder + to_mild + to_remission
  )

  # Current manuscript values, tagged the same way so they diff cleanly against candidates.
  therapy_map <- c(Ustekinumab = "UST", Infliximab = "IFX", "T-regulatory" = "TREG")
  current_tagged <- data.frame(
    therapy = therapy_map[current_induction$treatment],
    method = "manuscript_current", cycle_weeks = NA_real_, confidence = "reference",
    to_moderate_severe = current_induction$to_moderate_severe,
    to_moderate_severe_responder = current_induction$to_moderate_severe_responder,
    to_mild = current_induction$to_mild,
    to_remission = current_induction$to_remission,
    note = "Current value in data/raw/manuscript_supplement1_induction_transition_probabilities.csv",
    row_sum_check = current_induction$to_moderate_severe + current_induction$to_moderate_severe_responder +
      current_induction$to_mild + current_induction$to_remission
  )
  current_tagged <- current_tagged[!is.na(current_tagged$therapy), ]

  induction_out <- rbind(induction_candidates, current_tagged)
  induction_out <- induction_out[order(induction_out$therapy, induction_out$cycle_weeks, induction_out$method), ]

  maintenance_out <- build_maintenance_reconciliation(params, maintenance_df)

  if (write_output) {
    dir.create(proc_dir, showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(induction_out, file.path(proc_dir, "derived_induction_transition_probabilities.csv"), row.names = FALSE)
    utils::write.csv(maintenance_out, file.path(proc_dir, "derived_maintenance_reconciliation.csv"), row.names = FALSE)
    write_derivation_notes(induction_out, maintenance_out, proc_dir)
  }

  list(induction = induction_out, maintenance = maintenance_out)
}

write_derivation_notes <- function(induction_out, maintenance_out, proc_dir) {
  at8 <- induction_out[induction_out$cycle_weeks %in% 8 & induction_out$method == "nested_hazard", ]
  ref <- induction_out[induction_out$method == "manuscript_current", ]
  diff_lines <- vapply(intersect(at8$therapy, ref$therapy), function(th) {
    a <- at8[at8$therapy == th, ]
    r <- ref[ref$therapy == th, ]
    sprintf(
      "- **%s** (nested_hazard @ 8-week cycle vs. current manuscript value): to_remission %.4f vs %.4f (diff %.4f); to_moderate_severe %.4f vs %.4f (diff %.4f)",
      th, a$to_remission, r$to_remission, abs(a$to_remission - r$to_remission),
      a$to_moderate_severe, r$to_moderate_severe, abs(a$to_moderate_severe - r$to_moderate_severe)
    )
  }, character(1))

  maint_lines <- sprintf(
    "- **%s** Remission->Remission (8-week cycle): re-derived %.4f vs. model %.4f (diff %.4f), from '%s' (%d-week endpoint = %.3f)",
    maintenance_out$therapy, maintenance_out$derived_probability_per_8wk_cycle,
    maintenance_out$model_probability_per_8wk_cycle, maintenance_out$abs_diff,
    maintenance_out$trial_endpoint, maintenance_out$trial_endpoint_week, maintenance_out$trial_endpoint_value
  )

  notes <- c(
    "# Gate 1 DEALE re-derivation — findings and open gaps",
    "",
    sprintf("Generated by `R/00_derive_transition_probs.R` on %s.", format(Sys.Date())),
    "",
    "## Induction transition probabilities",
    "",
    "No single re-derivation method reproduces the current manuscript values within a tight",
    "tolerance for every therapy and destination state. Two structural ambiguities could not be",
    "resolved from the materials in this repository (Aliyev et al. 2019 Appendix S2, which",
    "documents the DEALE conversion methodology, is not included in `data/raw/`):",
    "",
    "1. **Cycle length at induction.** The worked example in `docs/analysis_plan.md` §7.3 converts",
    "   to a 2-week probability; `docs/analysis_plan.md` §6.1 describes an 8-week induction period.",
    "   Both are computed here (`cycle_weeks` column) since the correct one for reproducing the",
    "   current manuscript values is not stated in either the workbook or the plan.",
    "2. **Whether the remitter:responder ratio itself is DEALE-converted** over the trial window",
    "   before being applied (`nested_hazard`), or applied unconverted at the target cycle length",
    "   (`response_hazard_only`). Both are computed.",
    "",
    "Comparison against `data/raw/manuscript_supplement1_induction_transition_probabilities.csv`",
    "(nested_hazard method, 8-week cycle — the closer of the two methods to the induction period",
    "the plan describes):",
    "",
    paste(diff_lines, collapse = "\n"),
    "",
    "**IFX has no directly observed induction trial endpoint** in `data/raw/`; its candidate is",
    "derived indirectly from ADA's endpoint via the workbook's adjustment ratios, whose exact",
    "semantics (odds ratio vs. risk ratio; applied before or after DEALE conversion) are also not",
    "documented here. Treat the IFX induction candidate as low-confidence.",
    "",
    "**ADA has no current manuscript value to diff against** — it is not yet in the induction data",
    "(Decision 2, recorded 2026-08-04, adds it as a comparator; the candidate values here are the",
    "first ADA induction transition probabilities computed for this project).",
    "",
    "## Maintenance transition probabilities",
    "",
    "Only the Remission -> Remission cell reconciles cleanly against a directly comparable trial",
    "endpoint (\"remission probability at week 44, given remission at baseline\"), using a",
    "per-cycle retention-root conversion rather than DEALE (this endpoint is a multi-cycle",
    "cumulative retention probability, not a single-transition hazard):",
    "",
    paste(maint_lines, collapse = "\n"),
    "",
    "IFX and ADA maintenance trial endpoints in `data/raw/` are structured differently (an overall",
    "week-52 response probability + remitter ratio, not a baseline-conditional remission",
    "probability) and are not reconciled by this script. Reconciling them, and the rest of the",
    "maintenance matrix (Mild/M-SR/M-S/Surgery transitions) beyond the single Remission ->",
    "Remission cell, requires either the original Appendix S2 text or co-author input on the",
    "intended method — flagged here rather than guessed at.",
    "",
    "## Gate 1 status",
    "",
    "**Not closed.** Per `docs/analysis_plan.md` §14, Gate 1 requires re-derived probabilities to",
    "reconcile with the current values, *or* discrepancies to be explained in writing — this file",
    "is that explanation. Recommend review by the model lead before Gate 2 (engine build) begins:",
    "either source Appendix S2 to resolve the two structural ambiguities above, or make an",
    "explicit documented choice between the candidate methods computed here."
  )
  writeLines(notes, file.path(proc_dir, "DERIVATION_NOTES.md"))
}

if (sys.nframe() == 0L) {
  result <- run_deale_rederivation()
  cat("Wrote data/processed/derived_induction_transition_probabilities.csv,",
      "data/processed/derived_maintenance_reconciliation.csv,",
      "data/processed/DERIVATION_NOTES.md\n")
}
