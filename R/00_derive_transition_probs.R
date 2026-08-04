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

# Repo-root-relative, matching this project's convention of running scripts from the repo root
# (README.md: `source("analysis/run_full_analysis.R")`). Guarded so this file can also be
# source()'d from tests/testthat/ (a different working directory) after the caller has already
# loaded R/utils/transition_matrix.R itself via a repo-root-relative path -- see
# tests/testthat/test-derive-transition-probs.R.
if (!exists("validate_row_sums")) source("R/utils/transition_matrix.R")

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
    "- **%s**: to_moderate_severe %.4f, to_moderate_severe_responder %.4f, to_mild %.4f, to_remission %.4f",
    induction_out$therapy, induction_out$to_moderate_severe, induction_out$to_moderate_severe_responder,
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
    "## Induction (Aliyev Supplementary Table 3, used as published)",
    "",
    "No conversion needed -- Table 3 is already the terminal end-of-induction split: once a",
    "patient reaches Mild/M-SR/Remission the published matrix holds them there with probability",
    "1, showing this is a single-step absorbing partition, not a repeated-cycle process. ADA and",
    "IFX are published as numerically identical: Aliyev's base case sets the IFX:ADA efficacy",
    "ratio to 1.00 (varied 0.8-1.2 only in PSA), not a transcription error -- independently",
    "confirmed against the table image.",
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
