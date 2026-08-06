# Load and validate the published induction and maintenance transition probabilities from
# Aliyev, Hay & Hwang 2019 (Pharmacotherapy), Appendix S2, Supplementary Tables 3 and 4 -- the
# primary source this study's biologic-arm parameterisation rests on.
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
# HISTORY, 2026-08-04 (second revision): with Aliyev et al. 2019 Appendix S2 now available
# (data/raw/aliyev2019_appendixS2_table3_induction_transition_probabilities.csv,
# ..._table4_maintenance_transition_probabilities.csv -- transcribed from the appendix and
# verified by direct visual comparison against its table images, one row group at a time), this
# script's job became loading and validating those tables, plus one genuine conversion: Table 4
# is reported at a 2-week cycle length, but this study's original design (analysis_plan.md §4.1)
# picked an 8-week maintenance cycle to align with UST/IFX's q8w dosing, so this script converted
# Table 4's 2-week matrix to 8-week via exact Markov-chain matrix power.
#
# HISTORY, 2026-08-04 (third revision, this one): that 8-week alignment is dropped. The
# rationale for it -- matching the transition-probability cycle to biologic dosing cadence --
# doesn't buy anything once drug-administration costs are attached at whatever cycles match real
# dosing regardless of the underlying transition-probability cycle length, and Treg's own cost is
# a decision-tree-level one-time/two-dose event, not a recurring per-cycle charge. This study now
# runs on Aliyev's native 2-week cycle throughout: Table 4 is used as published, with no
# conversion, same as Table 3. This also means the engine (R/02_markov_engine.R) runs Aliyev's
# actual matrices unmodified, which directly strengthens Aim 5 (external validation against
# Aliyev's own published results).
#
# HISTORY, 2026-08-06 (B1 fix, peer review 2026-08-05): the third revision's own claim above --
# "Table 3 is already the terminal end-of-induction split... a single-step absorbing partition,
# not a repeated-cycle process" -- is WRONG, and every result in this repository before this fix
# is superseded (docs/model_audit_v6.md A17). Appendix S2's own worked example is explicit that
# Table 3's Moderate-Severe row is a PER-2-WEEK-CYCLE transition probability derived FROM a
# week-6/week-4 trial endpoint via a rate conversion, not the endpoint itself:
#   "Ustekinumab week 6 remission probability -> rate = -ln(1-0.349)/(6/52) = 3.720/yr ->
#    Probability (2 weeks) = 1 - exp(-3.720/26) = 0.133 -- a 2-week Moderate-Severe to Remission
#    transition probability."
# 0.133 is exactly Table 3's published UST Moderate-Severe -> Remission entry -- confirming it
# must be applied REPEATEDLY (3 times for UST's week-6 endpoint, 2 times for IFX/ADA's week-4
# endpoint), not read once, to reproduce the trial result it was derived from:
# 1 - (1-0.133)^3 = 0.348 vs. published UNITI week-6 remission 0.349 (see
# tests/testthat/test-external-validity.R). The "absorbing" self-loops on Mild/M-SR/Remission that
# make the published table square (Induction Assumption #3: "any patient that transitions out of
# the Moderate-Severe state during induction remains in the new state until the end of induction")
# are exactly what a MULTI-cycle induction with within-phase absorbing states looks like -- they
# are not evidence the table is a one-step process. `load_published_induction()` below now runs
# each therapy's full induction matrix through `simulate_cohort()` (R/utils/transition_matrix.R)
# for `INDUCTION_CYCLES[[therapy]]` cycles and reads the terminal row, instead of reading the
# Moderate-Severe row's raw published entries directly as the end-of-induction split.

# Repo-root-relative, matching this project's convention of running scripts from the repo root
# (README.md: `source("analysis/run_full_analysis.R")`). Guarded so this file can also be
# source()'d from tests/testthat/ (a different working directory) after the caller has already
# loaded R/utils/transition_matrix.R itself via a repo-root-relative path -- see
# tests/testthat/test-derive-transition-probs.R.
if (!exists("validate_row_sums")) source("R/utils/transition_matrix.R")

# ---- Induction: Table 3, applied for the number of 2-week cycles that reproduces its endpoint --

#' Number of native 2-week cycles Table 3's Moderate-Severe row must be applied for to reach each
#' therapy's own trial induction endpoint (module header, B1 fix, 2026-08-06). UST's induction
#' trial (UNITI) reads out at week 6 = 3 cycles; IFX and ADA's own induction trials (Aliyev's
#' IFX:ADA efficacy ratio = 1.00 base case, R/00's third-revision note above) read out at week 4 =
#' 2 cycles -- both directly against the endpoint labels in
#' data/raw/aliyev2019_appendixS1_table2_parameters.csv ("UST Week 6 Response Probability", "ADA
#' Week 4 Response Probability"). Verified by two independent reconstructions of the published
#' endpoint, not just this table's own worked example (see test-external-validity.R): UST,
#' `1-(1-0.133)^3 = 0.348` vs. published UNITI week-6 remission 0.349; IFX/ADA,
#' `1-(1-0.203)^2 = 0.365` vs. an ADA week-4 remission of 0.50 response x 0.72 remitter:responder =
#' 0.360 computed from the same appendix's raw trial-endpoint table.
INDUCTION_CYCLES <- c(UST = 3, IFX = 2, ADA = 2)

#' The four induction states Table 3 actually reports rows for -- deliberately not the full
#' 6-state MAINTENANCE_STATES universe (Table 3 has no Surgery/Death entries; induction has no
#' route to either, R/01_decision_tree.R's own module header) -- so build_transition_matrix()'s
#' padding never has to invent an induction-phase Surgery/Death row that doesn't exist in Aliyev's
#' own table.
INDUCTION_STATES <- c("Moderate-Severe", "Moderate-Severe Responder", "Mild", "Remission")

load_published_induction <- function(raw_dir) {
  df <- utils::read.csv(
    file.path(raw_dir, "aliyev2019_appendixS2_table3_induction_transition_probabilities.csv"),
    stringsAsFactors = FALSE
  )
  # Published, rounded-to-3-4-significant-figures values (e.g. 0.791 + 0.0377 + 0.0377 + 0.133 =
  # 0.9994) -- a looser tolerance than the default is expected and correct here, not a bug.
  validate_row_sums(df, by = "treatment", tol = 0.001)

  # B1 fix (module header, 2026-08-06): Table 3's Moderate-Severe row is a PER-2-WEEK-CYCLE
  # transition probability, not the terminal end-of-induction split -- it must be run through
  # simulate_cohort() for INDUCTION_CYCLES[[therapy]] cycles, starting from 100% Moderate-Severe,
  # to reach the split this project's induction decision tree (R/01_decision_tree.R) actually
  # needs. The other three rows (Moderate-Severe Responder/Mild/Remission -> self, probability 1)
  # are the absorbing within-induction states Aliyev's own Induction Assumption #3 describes
  # ("any patient that transitions out of the Moderate-Severe state during induction remains in
  # the new state until the end of induction") -- they make the repeated multiplication correct
  # (once a patient responds, further induction cycles leave them exactly where they are), not
  # evidence the whole table is a single step (the third-revision comment this fix retracts,
  # module header).
  therapies <- unique(df$treatment)
  stopifnot(all(therapies %in% names(INDUCTION_CYCLES)))
  rows <- lapply(therapies, function(tx) {
    tx_df <- df[df$treatment == tx, ]
    m <- build_transition_matrix(tx_df, INDUCTION_STATES)
    n_cycles <- INDUCTION_CYCLES[[tx]]
    initial <- stats::setNames(c(1, 0, 0, 0), INDUCTION_STATES)
    trace <- simulate_cohort(m, initial, n_cycles)
    terminal <- trace[n_cycles + 1, ]

    data.frame(
      therapy = tx,
      confidence = "published_source",
      to_moderate_severe = unname(terminal[["Moderate-Severe"]]),
      to_moderate_severe_responder = unname(terminal[["Moderate-Severe Responder"]]),
      to_mild = unname(terminal[["Mild"]]),
      to_remission = unname(terminal[["Remission"]]),
      induction_cycles = n_cycles,
      note = sprintf(
        paste(
          "Aliyev et al. 2019 Appendix S2, Supplementary Table 3 (verified against the table",
          "image; PHAR2208 supplementary materials): published 2-week Moderate-Severe transition",
          "probabilities applied for %d native 2-week cycles (B1 fix, 2026-08-06,",
          "docs/model_audit_v6.md A17) to reach the published week-%d induction endpoint -- not",
          "read directly as a 1-cycle terminal split."
        ),
        n_cycles, n_cycles * 2
      ),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# ---- Maintenance: Table 4, used as-is at Aliyev's native 2-week cycle ------

MAINTENANCE_STATES <- c("Moderate-Severe", "Moderate-Severe Responder", "Mild",
                         "Remission", "Surgery", "Death")

load_published_maintenance <- function(raw_dir) {
  df <- utils::read.csv(
    file.path(raw_dir, "aliyev2019_appendixS2_table4_maintenance_transition_probabilities.csv"),
    stringsAsFactors = FALSE
  )
  # Biologic arms are genuinely missing a Moderate-Severe row in the source (see
  # build_transition_matrix() for why); validate only the rows Aliyev actually reported. Same
  # published-rounding tolerance as the induction table.
  validate_row_sums(df, by = "treatment", tol = 0.001)

  out <- data.frame(
    therapy = df$treatment, from_state = df$from_state, to_state = df$to_state,
    probability = df$probability,
    note = "Aliyev et al. 2019 Appendix S2, Supplementary Table 4 (verified against the table image; PHAR2208 supplementary materials). Used as published at its native 2-week cycle, no conversion applied.",
    stringsAsFactors = FALSE
  )
  out
}

# ---- Orchestration ----------------------------------------------------------

run_derivation <- function(raw_dir = "data/raw", proc_dir = "data/processed", write_output = TRUE) {
  induction_out <- load_published_induction(raw_dir)
  maintenance_out <- load_published_maintenance(raw_dir)

  if (write_output) {
    dir.create(proc_dir, showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(induction_out, file.path(proc_dir, "derived_induction_transition_probabilities.csv"), row.names = FALSE)
    utils::write.csv(maintenance_out, file.path(proc_dir, "derived_maintenance_probabilities.csv"), row.names = FALSE)
    write_derivation_notes(induction_out, maintenance_out, proc_dir)
  }

  list(induction = induction_out, maintenance = maintenance_out)
}

write_derivation_notes <- function(induction_out, maintenance_out, proc_dir) {
  induction_lines <- sprintf(
    "- **%s** (%d induction cycles, week %d endpoint): to_moderate_severe %.4f, to_moderate_severe_responder %.4f, to_mild %.4f, to_remission %.4f",
    induction_out$therapy, induction_out$induction_cycles, induction_out$induction_cycles * 2,
    induction_out$to_moderate_severe, induction_out$to_moderate_severe_responder,
    induction_out$to_mild, induction_out$to_remission
  )

  maintenance_therapies <- unique(maintenance_out$therapy)
  remission_lines <- sprintf(
    "- **%s** Remission -> Remission (2-week cycle, as published): %.4f",
    maintenance_therapies,
    vapply(maintenance_therapies, function(tx) {
      maintenance_out[
        maintenance_out$therapy == tx & maintenance_out$from_state == "Remission" &
          maintenance_out$to_state == "Remission",
        "probability"
      ]
    }, numeric(1))
  )

  notes <- c(
    "# Gate 1 transition probabilities: Aliyev et al. 2019 Appendix S2, native 2-week cycle",
    "",
    sprintf("Generated by `R/00_derive_transition_probs.R` on %s.", format(Sys.Date())),
    "",
    "**2026-08-06, B1 fix (peer review 2026-08-05, docs/model_audit_v6.md A17):** the third",
    "revision's own claim below that Table 3 is already the terminal end-of-induction split was",
    "WRONG -- every induction number produced before this fix is superseded. Table 3's",
    "Moderate-Severe row is a PER-2-WEEK-CYCLE transition probability derived FROM a week-6",
    "(UST)/week-4 (IFX/ADA) trial endpoint (Appendix S2's own worked example), not the endpoint",
    "itself, and must be applied repeatedly -- 3 cycles for UST, 2 for IFX/ADA -- to reach the",
    "actual end-of-induction split. `load_published_induction()` now does this via",
    "`simulate_cohort()`; the figures below are the corrected terminal values, not the raw",
    "published Moderate-Severe row.",
    "",
    "**2026-08-04, third revision:** this study now runs on Aliyev's native 2-week cycle",
    "throughout, dropping the earlier 8-week-cycle design (analysis_plan.md §4.1) and the",
    "matrix-power conversion it required. The rationale for 8-week -- aligning the",
    "transition-probability cycle with UST/IFX's q8w dosing cadence -- doesn't buy anything once",
    "drug-administration costs are attached at whatever cycles match real dosing regardless of",
    "the transition-probability cycle length, and Treg's cost is a decision-tree-level",
    "one-time/two-dose event rather than a recurring per-cycle charge. Both Table 3 (induction)",
    "and Table 4 (maintenance) are now used exactly as Aliyev published them, with no cycle-length",
    "conversion anywhere in this script. This also means `R/02_markov_engine.R` runs Aliyev's",
    "actual matrices unmodified, directly strengthening Aim 5 (external validation against",
    "Aliyev's own published results).",
    "",
    "The comparison against `data/processed/model_maintenance_transition_probabilities.csv`",
    "(the v6 workbook's own hardcoded, separately-derived 8-week snapshot) that a previous",
    "revision of this file computed is retired along with the 8-week design: this project now",
    "diverges from the workbook's cycle length by construction, not by approximation error, so a",
    "numeric diff against it is no longer informative. The workbook snapshot remains available",
    "for historical reference (`data/data_dictionary.md`).",
    "",
    "## Induction (Aliyev Supplementary Table 3, applied for induction_cycles native 2-week cycles)",
    "",
    "Table 3's Moderate-Severe row is Aliyev's own PER-2-WEEK-CYCLE transition probability, derived",
    "from a week-6 (UST) / week-4 (IFX/ADA) trial endpoint (Appendix S2's worked example) -- applied",
    "here for `induction_cycles` cycles (3 for UST, 2 for IFX/ADA) via `simulate_cohort()`, not read",
    "directly as a 1-cycle terminal split (B1 fix, 2026-08-06, docs/model_audit_v6.md A17,",
    "superseding this file's own third-revision claim to the contrary). Once a patient reaches",
    "Mild/M-SR/Remission the published matrix holds them there with probability 1 (Aliyev's own",
    "Induction Assumption #3) -- an absorbing state WITHIN a multi-cycle induction, not evidence the",
    "table is a single step. ADA and IFX are published as numerically identical: Aliyev's base case",
    "sets the IFX:ADA efficacy ratio to 1.00 (varied 0.8-1.2 only in PSA), not a transcription",
    "error -- independently confirmed against the table image.",
    "",
    paste(induction_lines, collapse = "\n"),
    "",
    "## Maintenance (Aliyev Supplementary Table 4, used as published at its native 2-week cycle)",
    "",
    "Biologic arms (UST/IFX/ADA) have no `Moderate-Severe` row in the published table by design:",
    "patients who deteriorate to Moderate-Severe exit to the CT track (analysis_plan.md §6.1)",
    "rather than continuing on the biologic's own matrix. `R/02_markov_engine.R` implements that",
    "switch explicitly; `build_transition_matrix()` (`R/utils/transition_matrix.R`) pads the",
    "missing row as an absorbing self-loop only for the purpose of representing each matrix in",
    "isolation, not as a modelling claim about what happens to those patients.",
    "",
    "Remission -> Remission per therapy, for reference:",
    "",
    paste(remission_lines, collapse = "\n"),
    "",
    "## Gate 1 status",
    "",
    "**Closed.** Both phases, all four biologic/CT therapies, sourced directly from Aliyev et al.",
    "2019 Appendix S2, validated for internal consistency, used at Aliyev's native 2-week cycle",
    "with no conversion. No open data gaps remain before Gate 2 (`R/02_markov_engine.R` onward)."
  )
  writeLines(notes, file.path(proc_dir, "DERIVATION_NOTES.md"))
}

if (sys.nframe() == 0L) {
  result <- run_derivation()
  cat("Wrote data/processed/derived_induction_transition_probabilities.csv,",
      "data/processed/derived_maintenance_probabilities.csv,",
      "data/processed/DERIVATION_NOTES.md\n")
}
