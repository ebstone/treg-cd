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

test_that("load_published_induction reproduces the verified Table 3 values and sums to 1", {
  out <- load_published_induction(repo_root_relative("data", "raw"))
  expect_setequal(out$therapy, c("UST", "IFX", "ADA"))

  ust <- out[out$therapy == "UST", ]
  expect_equal(ust$to_moderate_severe, 0.791)
  expect_equal(ust$to_remission, 0.133)

  # ADA and IFX are published as numerically identical (IFX:ADA efficacy ratio = 1.00 base case),
  # confirmed against the appendix table image -- not a transcription error.
  ada <- out[out$therapy == "ADA", ]
  ifx <- out[out$therapy == "IFX", ]
  expect_equal(ada$to_remission, ifx$to_remission)

  # Published values are rounded to 3-4 significant figures, so rows sum to ~1, not exactly 1.
  totals <- out$to_moderate_severe + out$to_moderate_severe_responder + out$to_mild + out$to_remission
  expect_true(all(abs(totals - 1) < 0.001))
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
