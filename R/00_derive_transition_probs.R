# Load, validate, and cycle-convert the published induction and maintenance transition
# probabilities from Aliyev, Hay & Hwang 2019 (Pharmacotherapy), Appendix S2, Supplementary
# Tables 3 and 4 -- the primary source this study's biologic-arm parameterisation rests on.
#
# HISTORY, 2026-08-04: earlier versions of this script tried to *re-derive* these values from
# scratch via DEALE conversion of the raw trial endpoints in
# data/raw/aliyev2019_appendixS1_table2_parameters.csv (response probability x
# remitter:responder ratio, per therapy). That attempt is retired: Appendix S2's own worked
# example ("Adjustment Ratio Calculation") shows Aliyev's actual method for IFX/ADA involves an
# indirect, placebo-anchored cross-trial adjustment -- converting each trial's PLACEBO-arm
# endpoint to a yearly rate, taking the ratio between trials, and applying it to the drug arm's
# own rate -- which requires each trial's separate placebo-arm endpoint. That endpoint is not
# present in data/raw/ (only the active-drug-arm endpoints are). Confirmed by direct comparison:
# neither prior candidate re-derivation reproduced Aliyev's published Table 3 (off by ~0.35 on
# UST to_remission, ~0.39 on ADA) -- not a rounding-level gap, a genuine missing input. See git
# history for the retired derive_induction()/derive_ifx_induction() functions if that
# from-scratch attempt is useful for the manuscript's discussion section.
#
# With Appendix S2 now available (data/raw/aliyev2019_appendixS2_table3_induction_transition_probabilities.csv,
# ..._table4_maintenance_transition_probabilities.csv -- transcribed from the appendix and
# verified by direct visual comparison against its table images, one row group at a time), the
# well-defined task left for this script is not re-derivation but VALIDATION and one genuine
# CONVERSION:
#
#   - Table 3 (induction) needs no conversion. Its structure is a single-step absorbing
#     partition -- once a patient reaches Mild/M-SR/Remission the matrix holds them there with
#     probability 1 -- showing it is already the terminal end-of-induction split, used as-is
#     regardless of the 2-week-equivalent hazard framing Appendix S2 uses to compute it.
#   - Table 4 (maintenance) is reported at a 2-week cycle length (Appendix S2's DEALE worked
#     example; "the same method was used ... in induction and maintenance phases"), but this
#     study's redesigned maintenance Markov runs on 8-week cycles (analysis_plan.md §6.1).
#     Converting a verified 2-week transition matrix to an 8-week one is standard, exact
#     Markov-chain math (Chapman-Kolmogorov): M_8wk = M_2wk %*% M_2wk %*% M_2wk %*% M_2wk. No
#     guessing required, unlike the per-endpoint DEALE tricks the retired approach relied on.

# ---- Validation -----------------------------------------------------------

#' Check that a long-format transition table (from_state, to_state, probability, grouped by
#' `by`) has, for every group and from_state, probabilities in [0,1] summing to 1.
validate_row_sums <- function(df, by, tol = 1e-6) {
  key <- do.call(paste, c(df[c(by, "from_state")], sep = ""))
  sums <- tapply(df$probability, key, sum)
  bad <- sums[abs(sums - 1) > tol]
  if (length(bad) > 0) {
    stop("Row sums not equal to 1 for: ", paste(names(bad), "=", round(bad, 6), collapse = "; "))
  }
  if (any(df$probability < 0 | df$probability > 1)) {
    stop("Probabilities outside [0,1] found in transition table")
  }
  invisible(TRUE)
}

#' Build a square transition matrix (states x states) from a long-format
#' (from_state, to_state, probability) data frame for one therapy. Any from_state present in
#' `states` but absent from the data (biologic arms have no "Moderate-Severe" row: patients who
#' deteriorate to Moderate-Severe exit to the CT track rather than continuing on this matrix,
#' by this study's own model design, analysis_plan.md §6.1) is padded as an absorbing
#' self-loop -- the correct convention for computing "this matrix's own n-cycle behaviour",
#' since the real cross-matrix CT-switch is the Markov engine's job (Gate 2), not this script's.
build_transition_matrix <- function(df, states) {
  m <- matrix(0, nrow = length(states), ncol = length(states), dimnames = list(states, states))
  for (i in seq_len(nrow(df))) {
    m[df$from_state[i], df$to_state[i]] <- df$probability[i]
  }
  missing_from <- states[rowSums(m) == 0]
  for (s in missing_from) m[s, s] <- 1
  m
}

#' Convert a square transition matrix to its n-cycle equivalent (Chapman-Kolmogorov).
matrix_power <- function(m, n) {
  stopifnot(n >= 1, nrow(m) == ncol(m))
  result <- m
  if (n > 1) for (i in seq_len(n - 1)) result <- result %*% m
  result
}

#' Long-format (from_state, to_state, probability) rows from a square matrix.
matrix_to_long <- function(m, therapy) {
  states <- rownames(m)
  do.call(rbind, lapply(states, function(from) {
    data.frame(therapy = therapy, from_state = from, to_state = states,
               probability = unname(m[from, ]), stringsAsFactors = FALSE)
  }))
}

# ---- Induction: Table 3, used as-is ----------------------------------------

load_published_induction <- function(raw_dir) {
  df <- utils::read.csv(
    file.path(raw_dir, "aliyev2019_appendixS2_table3_induction_transition_probabilities.csv"),
    stringsAsFactors = FALSE
  )
  # Published, rounded-to-3-4-significant-figures values (e.g. 0.791 + 0.0377 + 0.0377 + 0.133 =
  # 0.9994) -- a looser tolerance than the default is expected and correct here, not a bug.
  validate_row_sums(df, by = "treatment", tol = 0.001)

  # Only the Moderate-Severe row is a live induction split; the other rows (Moderate-Severe
  # Responder/Mild/Remission -> self, probability 1) are the absorbing placeholders that make
  # the published table square, not additional induction-phase information.
  out <- df[df$from_state == "Moderate-Severe", c("treatment", "to_state", "probability")]
  wide <- reshape(out, idvar = "treatment", timevar = "to_state", direction = "wide")
  names(wide) <- sub("^probability\\.", "", names(wide))
  data.frame(
    therapy = wide$treatment,
    confidence = "published_source",
    to_moderate_severe = wide[["Moderate-Severe"]],
    to_moderate_severe_responder = wide[["Moderate-Severe Responder"]],
    to_mild = wide[["Mild"]],
    to_remission = wide[["Remission"]],
    note = "Aliyev et al. 2019 Appendix S2, Supplementary Table 3 (verified against the table image; PHAR2208 supplementary materials). Used as published, no conversion applied.",
    stringsAsFactors = FALSE
  )
}

# ---- Maintenance: Table 4, 2-week -> 8-week by matrix power ----------------

MAINTENANCE_STATES <- c("Moderate-Severe", "Moderate-Severe Responder", "Mild",
                         "Remission", "Surgery", "Death")

load_published_maintenance_8wk <- function(raw_dir, cycles_per_8wk = 4) {
  df <- utils::read.csv(
    file.path(raw_dir, "aliyev2019_appendixS2_table4_maintenance_transition_probabilities.csv"),
    stringsAsFactors = FALSE
  )
  # Biologic arms are genuinely missing a Moderate-Severe row in the source (see
  # build_transition_matrix() for why); validate only the rows Aliyev actually reported. Same
  # published-rounding tolerance as the induction table.
  validate_row_sums(df, by = "treatment", tol = 0.001)

  therapies <- unique(df$treatment)
  rows <- lapply(therapies, function(tx) {
    m_2wk <- build_transition_matrix(df[df$treatment == tx, ], MAINTENANCE_STATES)
    m_8wk <- matrix_power(m_2wk, cycles_per_8wk)
    matrix_to_long(m_8wk, tx)
  })
  out <- do.call(rbind, rows)
  # Rounding error in the published 2-week values compounds slightly across four
  # matrix multiplications; still well within a sane tolerance for published-data rounding.
  validate_row_sums(out, by = "therapy", tol = 0.005)
  out$note <- "8-week probability computed via M_2wk^4 (Chapman-Kolmogorov) from Aliyev et al. 2019 Appendix S2, Supplementary Table 4 (2-week cycle, verified against the table image). Moderate-Severe treated as absorbing for biologic arms (they exit to the CT track in this study's model, analysis_plan.md §6.1; not a value Aliyev reports)."
  out
}

#' Compare the one maintenance cell (Remission -> Remission) that means the same thing in both
#' this study's redesigned model and the existing v6-workbook snapshot -- the workbook uses a
#' 7-state structure (an extra "Mod/Sev Non-Resp" destination not in Aliyev's 6-state Table 4),
#' so a full-matrix comparison isn't attempted; Remission -> Remission is unambiguous in both.
compare_remission_retention <- function(published_8wk, workbook_snapshot) {
  therapies <- intersect(unique(published_8wk$therapy), unique(workbook_snapshot$therapy))
  rows <- lapply(therapies, function(tx) {
    derived <- published_8wk[
      published_8wk$therapy == tx & published_8wk$from_state == "Remission" &
        published_8wk$to_state == "Remission",
      "probability"
    ]
    workbook <- workbook_snapshot[
      workbook_snapshot$therapy == tx & workbook_snapshot$scenario == "base" &
        workbook_snapshot$from_state == "Remission" & workbook_snapshot$to_state == "Remission",
      "probability_per_8wk_cycle"
    ]
    if (length(workbook) == 0) return(NULL)
    data.frame(
      therapy = tx, from_state = "Remission", to_state = "Remission",
      derived_probability_per_8wk_cycle = derived,
      workbook_snapshot_probability_per_8wk_cycle = workbook,
      abs_diff = abs(derived - workbook),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

# ---- Orchestration ----------------------------------------------------------

run_derivation <- function(raw_dir = "data/raw", proc_dir = "data/processed", write_output = TRUE) {
  induction_out <- load_published_induction(raw_dir)

  maintenance_8wk <- load_published_maintenance_8wk(raw_dir)
  workbook_snapshot <- utils::read.csv(
    file.path(proc_dir, "model_maintenance_transition_probabilities.csv"),
    stringsAsFactors = FALSE
  )
  maintenance_comparison <- compare_remission_retention(maintenance_8wk, workbook_snapshot)

  if (write_output) {
    dir.create(proc_dir, showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(induction_out, file.path(proc_dir, "derived_induction_transition_probabilities.csv"), row.names = FALSE)
    utils::write.csv(maintenance_8wk, file.path(proc_dir, "derived_maintenance_probabilities.csv"), row.names = FALSE)
    write_derivation_notes(induction_out, maintenance_comparison, proc_dir)
  }

  list(induction = induction_out, maintenance_8wk = maintenance_8wk, maintenance_comparison = maintenance_comparison)
}

write_derivation_notes <- function(induction_out, maintenance_comparison, proc_dir) {
  induction_lines <- sprintf(
    "- **%s**: to_moderate_severe %.4f, to_moderate_severe_responder %.4f, to_mild %.4f, to_remission %.4f",
    induction_out$therapy, induction_out$to_moderate_severe, induction_out$to_moderate_severe_responder,
    induction_out$to_mild, induction_out$to_remission
  )

  maint_lines <- sprintf(
    "- **%s** Remission->Remission (8-week cycle): matrix-power-derived %.4f vs. v6-workbook snapshot %.4f (diff %.4f)",
    maintenance_comparison$therapy, maintenance_comparison$derived_probability_per_8wk_cycle,
    maintenance_comparison$workbook_snapshot_probability_per_8wk_cycle, maintenance_comparison$abs_diff
  )

  notes <- c(
    "# Gate 1 transition probabilities: published source, validated and cycle-converted",
    "",
    sprintf("Generated by `R/00_derive_transition_probs.R` on %s.", format(Sys.Date())),
    "",
    "**2026-08-04, second revision:** with Aliyev et al. 2019 Appendix S2 now available",
    "(`data/raw/aliyev2019_appendixS2_table3_induction_transition_probabilities.csv`,",
    "`..._table4_maintenance_transition_probabilities.csv`, transcribed from the appendix and",
    "verified by direct visual comparison against its table images), this script no longer",
    "attempts a from-scratch re-derivation from raw trial endpoints. Appendix S2's own worked",
    "example shows Aliyev's actual method for IFX/ADA requires each trial's separate",
    "placebo-arm endpoint (not present in `data/raw/aliyev2019_appendixS1_table2_parameters.csv`,",
    "which only has the active-drug-arm endpoints) -- a genuine missing input, not a modelling",
    "choice this project can resolve by picking a method. The values below are Aliyev's own",
    "published, peer-reviewed numbers: validated for internal consistency (row sums = 1,",
    "probabilities in [0,1]) and, for maintenance, cycle-converted from the published 2-week",
    "matrices to this study's 8-week cycle via exact Markov-chain matrix power -- not guessed.",
    "",
    "## Induction (Aliyev Supplementary Table 3, used as published)",
    "",
    "No conversion needed -- Table 3 is already the terminal end-of-induction split (see header",
    "comment for why). ADA and IFX are published as numerically identical: Aliyev's base case",
    "sets the IFX:ADA efficacy ratio to 1.00 (varied 0.8-1.2 only in PSA), not a transcription",
    "error -- independently confirmed against the table image.",
    "",
    paste(induction_lines, collapse = "\n"),
    "",
    "## Maintenance (Aliyev Supplementary Table 4, 2-week -> 8-week via matrix power)",
    "",
    "Compared against the existing `data/processed/model_maintenance_transition_probabilities.csv`",
    "(a snapshot of the v6 workbook's own, separately hardcoded 8-week probabilities, whose",
    "conversion-from-2-week method is undocumented in the workbook -- `data/data_dictionary.md`).",
    "Only Remission -> Remission is compared: the workbook uses a 7-state structure with an",
    "extra `Mod/Sev Non-Resp` destination not present in Aliyev's 6-state Table 4, so a full-matrix",
    "diff isn't well-defined; Remission -> Remission means the same thing in both.",
    "",
    paste(maint_lines, collapse = "\n"),
    "",
    "The gap for UST/IFX/CT (the only three the v6 workbook has -- ADA is not in it; workbook",
    "consistently higher than the matrix-power-derived value)",
    "is a real, now-precisely-quantified discrepancy between the v6 workbook's undocumented",
    "8-week conversion and the mathematically exact conversion of Aliyev's own published 2-week",
    "matrix. Worth a model-lead read before Gate 2 decides which maintenance matrix the rebuilt",
    "engine consumes -- this script does not decide that; it only makes the two comparable.",
    "",
    "## What remains open before Gate 2",
    "",
    "None of the previous \"data gap\" items (IFX induction, IFX/ADA maintenance) are open any",
    "longer -- Appendix S2 supplied all four therapies for both phases directly. What remains is",
    "a decision, not a data gap: whether the rebuilt engine (`R/02_markov_engine.R` onward)",
    "consumes these matrix-power-derived 8-week maintenance probabilities, or the v6 workbook's",
    "existing 8-week snapshot, given the gap quantified above.",
    "",
    "## Gate 1 status",
    "",
    "**Transition probabilities: closed.** Both phases, all four biologic/CT therapies, sourced",
    "directly from Aliyev et al. 2019 Appendix S2, validated for internal consistency, and",
    "maintenance cycle-converted by exact Markov-chain math. The remaining maintenance-matrix",
    "choice above is a Gate 2 engine-design decision, not a Gate 1 blocker."
  )
  writeLines(notes, file.path(proc_dir, "DERIVATION_NOTES.md"))
}

if (sys.nframe() == 0L) {
  result <- run_derivation()
  cat("Wrote data/processed/derived_induction_transition_probabilities.csv,",
      "data/processed/derived_maintenance_probabilities.csv,",
      "data/processed/DERIVATION_NOTES.md\n")
}
