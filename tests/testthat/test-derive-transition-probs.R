# testthat::test_dir() runs with the working directory set to tests/testthat/, but
# R/00_derive_transition_probs.R's top-level source("R/utils/transition_matrix.R") and
# run_derivation()'s default raw_dir/proc_dir are repo-root-relative. Load utils first via a
# repo-root-relative path so R/00's own `if (!exists("validate_row_sums"))` guard skips its
# internal (otherwise-broken-from-here) source() call, then source R/00 itself.
repo_root_relative <- function(...) file.path("..", "..", ...)
source(repo_root_relative("R", "utils", "transition_matrix.R"), local = TRUE)
source(repo_root_relative("R", "00_derive_transition_probs.R"), local = TRUE)

run_from_repo_root <- function(...) {
  run_derivation(
    raw_dir = repo_root_relative("data", "raw"),
    proc_dir = repo_root_relative("data", "processed"),
    ...
  )
}

test_that("load_published_induction reproduces the verified Table 3 endpoint (B1 fix) and sums to 1", {
  out <- load_published_induction(repo_root_relative("data", "raw"))
  expect_setequal(out$therapy, c("UST", "IFX", "ADA"))

  # B1 fix (docs/model_audit_v6.md A17): Table 3's Moderate-Severe row is a per-2-week-cycle
  # transition probability, applied INDUCTION_CYCLES[[therapy]] times, not read directly as a
  # 1-cycle terminal split -- the raw published row entries (0.791/0.133 for UST) are no longer
  # what this function returns.
  ust <- out[out$therapy == "UST", ]
  expect_equal(ust$induction_cycles, 3)
  expect_false(isTRUE(all.equal(ust$to_moderate_severe, 0.791)))
  expect_false(isTRUE(all.equal(ust$to_remission, 0.133)))

  ifx <- out[out$therapy == "IFX", ]
  expect_equal(ifx$induction_cycles, 2)

  # ADA and IFX are published as numerically identical (IFX:ADA efficacy ratio = 1.00 base case),
  # confirmed against the appendix table image -- not a transcription error. Still true of the
  # corrected multi-cycle terminal split, since both share the same raw Table 3 row and cycle count.
  ada <- out[out$therapy == "ADA", ]
  expect_equal(ada$induction_cycles, ifx$induction_cycles)
  expect_equal(ada$to_remission, ifx$to_remission)

  # Published values are rounded to 3-4 significant figures, so a single 2-week-cycle row sums to
  # ~1, not exactly 1 -- but repeated multiplication through a row-normalised matrix
  # (build_transition_matrix()) conserves cohort mass exactly, so the corrected terminal split
  # sums to 1 far more tightly than the raw table's own rounding tolerance.
  totals <- out$to_moderate_severe + out$to_moderate_severe_responder + out$to_mild + out$to_remission
  expect_equal(totals, rep(1, nrow(out)), tolerance = 1e-9)
})

test_that("load_published_induction's terminal split independently reproduces the closed-form repeated-multiplication result", {
  # Re-derives the expected terminal split from Table 3's own raw (un-normalised) Moderate-Severe
  # row using nothing but base-R vector arithmetic -- deliberately NOT calling
  # build_transition_matrix()/simulate_cohort() (the functions under test) -- so this is a real
  # independent check, not a tautology. Absorbing self-loops on Mild/M-SR/Remission mean the
  # closed form is exactly geometric: after n cycles, mass remaining in Moderate-Severe is a^n
  # (a = the row-normalised Moderate-Severe self-loop), and mass having flowed into state X is
  # x * (a^(n-1) + ... + a + 1), x = X's own row-normalised one-cycle transition probability.
  closed_form_terminal <- function(raw_row, n_cycles) {
    total <- sum(raw_row)
    a <- raw_row[["Moderate-Severe"]] / total
    geometric_sum <- sum(a^seq(0, n_cycles - 1))
    c(
      "Moderate-Severe" = a^n_cycles,
      "Moderate-Severe Responder" = (raw_row[["Moderate-Severe Responder"]] / total) * geometric_sum,
      "Mild" = (raw_row[["Mild"]] / total) * geometric_sum,
      "Remission" = (raw_row[["Remission"]] / total) * geometric_sum
    )
  }

  out <- load_published_induction(repo_root_relative("data", "raw"))
  raw_rows <- list(
    UST = c("Moderate-Severe" = 0.791, "Moderate-Severe Responder" = 0.0377, "Mild" = 0.0377, "Remission" = 0.133),
    IFX = c("Moderate-Severe" = 0.756, "Moderate-Severe Responder" = 0.0208, "Mild" = 0.0208, "Remission" = 0.203),
    ADA = c("Moderate-Severe" = 0.756, "Moderate-Severe Responder" = 0.0208, "Mild" = 0.0208, "Remission" = 0.203)
  )

  for (tx in names(raw_rows)) {
    row <- out[out$therapy == tx, ]
    expected <- closed_form_terminal(raw_rows[[tx]], row$induction_cycles)
    expect_equal(row$to_moderate_severe, unname(expected[["Moderate-Severe"]]), tolerance = 1e-9, info = tx)
    expect_equal(row$to_moderate_severe_responder, unname(expected[["Moderate-Severe Responder"]]), tolerance = 1e-9, info = tx)
    expect_equal(row$to_mild, unname(expected[["Mild"]]), tolerance = 1e-9, info = tx)
    expect_equal(row$to_remission, unname(expected[["Remission"]]), tolerance = 1e-9, info = tx)
  }
})

test_that("load_published_maintenance reproduces Table 4 unmodified at its native 2-week cycle", {
  out <- load_published_maintenance(repo_root_relative("data", "raw"))
  expect_setequal(out$therapy, c("UST", "IFX", "ADA", "CT"))
  expect_true(all(out$probability >= 0 & out$probability <= 1))
  # No cycle-length conversion applied: published rounding tolerance only, not the wider
  # compounded tolerance a matrix-power conversion would need.
  expect_true(validate_row_sums(out, by = "therapy", tol = 0.001))

  # Spot check against the table image (Supplementary Table 4, UST panel).
  ust_remission_remission <- out[
    out$therapy == "UST" & out$from_state == "Remission" & out$to_state == "Remission", "probability"
  ]
  expect_equal(ust_remission_remission, 0.982)
})

test_that("load_published_maintenance: biologic arms genuinely have no Moderate-Severe row, CT does", {
  # As published -- patients who deteriorate to Moderate-Severe exit to the CT track
  # (analysis_plan.md §6.1), so Aliyev's own Table 4 has no Moderate-Severe row for UST/IFX/ADA.
  # build_transition_matrix()'s absorbing-padding for this case is tested in
  # test-transition-matrix.R; this test just confirms the raw published data has the gap.
  out <- load_published_maintenance(repo_root_relative("data", "raw"))
  expect_false("Moderate-Severe" %in% out[out$therapy == "UST", "from_state"])
  expect_true("Moderate-Severe" %in% out[out$therapy == "CT", "from_state"])
})

test_that("run_derivation produces deterministic output", {
  first <- run_from_repo_root(write_output = FALSE)
  second <- run_from_repo_root(write_output = FALSE)
  expect_equal(first, second)
})
