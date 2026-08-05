# Parameter-provenance acceptance test (§1.4, docs/decision_resolutions_2026-08-05.md; A15/A16,
# docs/model_audit_v6.md). `data/processed/` was extracted from a live Excel workbook snapshot,
# and not everything in it is a parameter -- some of it is *state*: a live PSA draw or a scenario
# multiplier that happened to be sitting in a cell when the extraction ran. A15
# (model_dose_costs_and_psa_ranges.csv's `value_usd` column) and A16
# (model_health_utilities.csv's Remission figure) are both instances of exactly this bug, and both
# went undetected until someone happened to cross-check two files against each other by hand. This
# test makes that check permanent and mechanical rather than relying on it happening again by luck:
#
#   1. Every data/processed/ file that a deterministic base-case code path (R/04, R/05) actually
#      reads from disk must appear in PROCESSED_FILE_PROVENANCE below, with a one-line note of
#      what data/raw/ source or literature citation justifies treating its values as parameters
#      rather than snapshot state. Adding a new deterministic read.csv()/readLines() call against
#      `proc_dir` without adding (and actually checking) a matching registry entry fails this test.
#   2. Files already found to be snapshot artefacts (A15's `value_usd` column, A16's
#      model_health_utilities.csv) are explicitly retired -- this test fails if deterministic code
#      starts reading them again.

repo_root_relative <- function(...) file.path("..", "..", ...)
source(repo_root_relative("R", "utils", "transition_matrix.R"), local = TRUE)
source(repo_root_relative("R", "utils", "life_table.R"), local = TRUE)
source(repo_root_relative("R", "00_derive_transition_probs.R"), local = TRUE)
source(repo_root_relative("R", "02_markov_engine.R"), local = TRUE)
source(repo_root_relative("R", "03_cure_fraction_module.R"), local = TRUE)
source(repo_root_relative("R", "04_costs_utilities.R"), local = TRUE)

# ---- Registry: every data/processed/ file a deterministic code path reads, and why it's trusted --
PROCESSED_FILE_PROVENANCE <- list(
  "model_psa_parameter_distributions.csv" = paste(
    "Base values (0.82 Remission, 0.89/0.89/0.87/0.87 ratios) trace to Aliyev / NICE TA352",
    "(analysis_plan.md §7.1 items 19-20) and reproduce the manuscript's own published utility",
    "table exactly when chained -- A16 resolution §1.1(d), docs/model_audit_v6.md."
  ),
  "model_tenham_derived_treg_dose_cost.csv" = paste(
    "Full ten Ham et al. 2020 manufacturing-cost build-up (data/raw/",
    "tenham2020_manufacturing_cost_casestudies.csv) reproduced end to end -- analysis_plan.md",
    "§7.1 item 12, status 'Derived'."
  )
)

# Files a previous audit finding determined ARE snapshot state, not parameters -- must never be
# read by a deterministic (non-PSA-sampling) code path again.
RETIRED_PROCESSED_FILES <- c(
  "model_health_utilities.csv" # A16: Remission figure was a stray live PSA draw, not a parameter.
)

#' Every `proc_dir`-relative CSV path referenced by a deterministic-results source file, extracted
#' from the file's own text rather than by executing it -- this test is about what the code reads,
#' not what one particular call happens to pass.
processed_csvs_referenced_in <- function(path) {
  src <- readLines(path, warn = FALSE)
  hits <- regmatches(src, regexpr('proc_dir,\\s*"([^"]+)"', src))
  gsub('.*"([^"]+)".*', "\\1", hits[nzchar(hits)])
}

test_that("every data/processed/ file read deterministically by R/04 or R/05 is in the provenance registry", {
  files_to_scan <- c(
    repo_root_relative("R", "04_costs_utilities.R"),
    repo_root_relative("R", "05_deterministic_results.R")
  )
  referenced <- unique(unlist(lapply(files_to_scan, processed_csvs_referenced_in)))
  expect_true(length(referenced) > 0)  # sanity check the extraction itself isn't silently finding nothing

  for (f in referenced) {
    expect_true(
      f %in% names(PROCESSED_FILE_PROVENANCE),
      info = paste0(
        f, " is read as a deterministic input from data/processed/ but has no provenance ",
        "registry entry above -- add one (after actually verifying the provenance, not as a ",
        "placeholder) before this can pass."
      )
    )
    expect_false(
      f %in% RETIRED_PROCESSED_FILES,
      info = paste0(f, " was retired as a snapshot artefact and must not be read deterministically.")
    )
  }
})

test_that("load_health_state_utilities() no longer reads the retired A16 snapshot file", {
  src <- paste(readLines(repo_root_relative("R", "04_costs_utilities.R"), warn = FALSE), collapse = "\n")
  fn_body <- sub(".*load_health_state_utilities <- function[^{]*\\{(.*?)\\n}\\n.*", "\\1", src)
  expect_false(grepl("model_health_utilities.csv", fn_body, fixed = TRUE))
})

test_that("the deterministic utility vector is the same central chain R/06_psa.R's PSA samples around", {
  proc_dir <- repo_root_relative("data", "processed")
  df <- utils::read.csv(file.path(proc_dir, "model_psa_parameter_distributions.csv"), stringsAsFactors = FALSE)
  value_for <- function(v) df$value[df$variable == v]

  remission <- value_for("Remission")
  mild <- remission * value_for("Mild:Remission Ratio")
  msr <- mild * value_for("M-SR:Mild Ratio")
  ms <- msr * value_for("M-S:M-SR Ratio")
  surgery <- msr * value_for("Surgery:M-SR Ratio")

  utilities <- load_health_state_utilities(proc_dir)
  expect_equal(unname(utilities[["Remission"]]), remission)
  expect_equal(unname(utilities[["Mild"]]), mild)
  expect_equal(unname(utilities[["Moderate-Severe Responder"]]), msr)
  expect_equal(unname(utilities[["Moderate-Severe"]]), ms)
  expect_equal(unname(utilities[["Surgery"]]), surgery)

  # Regression guard against reverting to the retired snapshot figure (A16) -- e.g. by "fixing" a
  # future failing test back to the old value instead of recognising the base case changed.
  expect_false(isTRUE(all.equal(unname(utilities[["Remission"]]), 0.9554396355934719)))
})

test_that("load_treg_dose_acquisition_cost's psa_base-vs-value_usd choice (A15) stays on the traceable side", {
  # A15's resolution (docs/model_audit_v6.md): psa_base/the ten Ham derivation chain is
  # independently traceable; model_dose_costs_and_psa_ranges.csv's value_usd column is not, and is
  # never read by any deterministic function. This test doesn't re-derive the ten Ham chain (see
  # test-costs-utilities.R for that) -- it only guards that the untraceable column stays unused.
  files_to_scan <- c(
    repo_root_relative("R", "04_costs_utilities.R"),
    repo_root_relative("R", "05_deterministic_results.R")
  )
  referenced <- unique(unlist(lapply(files_to_scan, processed_csvs_referenced_in)))
  expect_false("model_dose_costs_and_psa_ranges.csv" %in% referenced)
})
