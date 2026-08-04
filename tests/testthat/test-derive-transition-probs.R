source(file.path("..", "..", "R", "00_derive_transition_probs.R"), local = TRUE)

# testthat::test_dir() runs with the working directory set to tests/testthat/, but
# run_deale_rederivation()'s default raw_dir/proc_dir are repo-root-relative — pass explicit
# repo-root-relative paths here rather than relying on the defaults.
repo_root_relative <- function(...) file.path("..", "..", ...)
run_from_repo_root <- function(...) {
  run_deale_rederivation(
    raw_dir = repo_root_relative("data", "raw"),
    proc_dir = repo_root_relative("data", "processed"),
    ...
  )
}

test_that("deale_convert reproduces the analysis_plan.md §7.3 worked example", {
  # UST week-6 remission 0.349 -> annual rate 3.720 -> 2-week probability 0.133
  expect_equal(deale_convert(0.349, 6, 2), 0.133, tolerance = 0.001)
})

test_that("deale_convert is the identity when t_from == t_to", {
  expect_equal(deale_convert(0.4, 8, 8), 0.4, tolerance = 1e-9)
})

test_that("deale_convert is monotonically increasing in t_to", {
  short <- deale_convert(0.5, 10, 2)
  long <- deale_convert(0.5, 10, 12)
  expect_lt(short, long)
})

test_that("deale_convert rejects out-of-range probabilities", {
  expect_error(deale_convert(-0.1, 6, 2))
  expect_error(deale_convert(1, 6, 2))
})

test_that("retention_root inverts a stationary per-cycle Markov process", {
  p_per_cycle <- 0.93
  n_cycles <- 5.5
  p_cumulative <- p_per_cycle^n_cycles
  expect_equal(retention_root(p_cumulative, n_cycles), p_per_cycle, tolerance = 1e-9)
})

test_that("run_deale_rederivation produces deterministic output", {
  first <- run_from_repo_root(write_output = FALSE)
  second <- run_from_repo_root(write_output = FALSE)
  expect_equal(first, second)
})

test_that("every computed induction candidate's destination probabilities sum to 1 and lie in [0,1]", {
  result <- run_from_repo_root(write_output = FALSE)
  computed <- result$induction[result$induction$confidence != "reference", ]

  expect_true(all(computed$row_sum_check >= 1 - 1e-9 & computed$row_sum_check <= 1 + 1e-9))

  dest_cols <- c("to_moderate_severe", "to_moderate_severe_responder", "to_mild", "to_remission")
  for (col in dest_cols) {
    expect_true(all(computed[[col]] >= 0 & computed[[col]] <= 1), info = col)
  }
})

test_that("the current manuscript reference rows are carried through unchanged", {
  result <- run_from_repo_root(write_output = FALSE)
  reference <- result$induction[result$induction$confidence == "reference", ]

  expect_setequal(reference$therapy, c("UST", "IFX", "TREG"))
  ust <- reference[reference$therapy == "UST", ]
  expect_equal(ust$to_remission, 0.133)
})

test_that("maintenance reconciliation is limited to UST and CT, with finite differences", {
  result <- run_from_repo_root(write_output = FALSE)
  expect_setequal(result$maintenance$therapy, c("UST", "CT"))
  expect_true(all(is.finite(result$maintenance$abs_diff)))
  expect_true(all(result$maintenance$derived_probability_per_8wk_cycle >= 0 &
    result$maintenance$derived_probability_per_8wk_cycle <= 1))
})
