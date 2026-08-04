source(file.path("..", "..", "R", "00_derive_transition_probs.R"), local = TRUE)

# testthat::test_dir() runs with the working directory set to tests/testthat/, but
# run_derivation()'s default raw_dir/proc_dir are repo-root-relative — pass explicit
# repo-root-relative paths here rather than relying on the defaults.
repo_root_relative <- function(...) file.path("..", "..", ...)
run_from_repo_root <- function(...) {
  run_derivation(
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

test_that("derive_induction is the identity decomposition when target_week == week (no extrapolation)", {
  out <- derive_induction("TEST", week = 6, p_response = 0.56, remitter_ratio = 0.63,
                           mild_split = 0.5, target_week = 6)
  expect_equal(out$to_remission, 0.56 * 0.63, tolerance = 1e-9)
  expect_equal(
    out$to_moderate_severe + out$to_moderate_severe_responder + out$to_mild + out$to_remission,
    1, tolerance = 1e-9
  )
})

test_that("derive_induction destination probabilities stay in [0,1] and sum to 1 when extrapolating forward", {
  out <- derive_induction("TEST", week = 6, p_response = 0.56, remitter_ratio = 0.63,
                           mild_split = 0.5, target_week = 8)
  total <- out$to_moderate_severe + out$to_moderate_severe_responder + out$to_mild + out$to_remission
  expect_equal(total, 1, tolerance = 1e-9)
  for (col in c("to_moderate_severe", "to_moderate_severe_responder", "to_mild", "to_remission")) {
    expect_true(out[[col]] >= 0 && out[[col]] <= 1, info = col)
  }
})

test_that("derive_induction fails loudly rather than returning a negative probability when over-extrapolated", {
  # ADA-shaped inputs (week 4, response 0.5, remitter ratio 0.72) are the tightest margin among
  # the therapies this script derives -- to_moderate_severe goes negative above target_week
  # ~10.2. This must error, not silently return an invalid probability that could corrupt a
  # downstream Markov engine.
  expect_error(
    derive_induction("ADA", week = 4, p_response = 0.5, remitter_ratio = 0.72,
                      mild_split = 0.5, target_week = 12),
    "invalid non-response probability"
  )
  # Comfortably within range at the actual base-case target_week (8) used by run_derivation().
  expect_no_error(
    derive_induction("ADA", week = 4, p_response = 0.5, remitter_ratio = 0.72,
                      mild_split = 0.5, target_week = 8)
  )
})

test_that("run_derivation produces deterministic output", {
  first <- run_from_repo_root(write_output = FALSE)
  second <- run_from_repo_root(write_output = FALSE)
  expect_equal(first, second)
})

test_that("every derived (non-legacy) induction row sums to 1 and lies in [0,1]", {
  result <- run_from_repo_root(write_output = FALSE)
  computed <- result$induction[result$induction$confidence != "legacy_reference_only", ]

  expect_true(all(computed$row_sum_check >= 1 - 1e-9 & computed$row_sum_check <= 1 + 1e-9))

  dest_cols <- c("to_moderate_severe", "to_moderate_severe_responder", "to_mild", "to_remission")
  for (col in dest_cols) {
    expect_true(all(computed[[col]] >= 0 & computed[[col]] <= 1), info = col)
  }
})

test_that("legacy rows are carried through unchanged as reference context", {
  result <- run_from_repo_root(write_output = FALSE)
  legacy <- result$induction[result$induction$confidence == "legacy_reference_only", ]

  expect_setequal(legacy$therapy, c("UST", "IFX", "TREG"))
  ust_legacy <- legacy[legacy$therapy == "UST", ]
  expect_equal(ust_legacy$to_remission, 0.133)
})

test_that("derived UST induction is an independent re-derivation, not a copy of the legacy value", {
  result <- run_from_repo_root(write_output = FALSE)
  derived_ust <- result$induction[result$induction$therapy == "UST" & result$induction$confidence == "direct", ]
  expect_false(isTRUE(all.equal(derived_ust$to_remission, 0.133, tolerance = 1e-6)))
})

test_that("IFX induction is flagged low-confidence (data gap, not resolved by this script)", {
  result <- run_from_repo_root(write_output = FALSE)
  ifx <- result$induction[result$induction$therapy == "IFX" & result$induction$confidence != "legacy_reference_only", ]
  expect_equal(ifx$confidence, "indirect_low_confidence")
})

test_that("maintenance derivation covers UST and CT with finite, in-range probabilities", {
  result <- run_from_repo_root(write_output = FALSE)
  expect_setequal(result$maintenance$therapy, c("UST", "CT"))
  expect_true(all(is.finite(result$maintenance$derived_probability_per_8wk_cycle)))
  expect_true(all(result$maintenance$derived_probability_per_8wk_cycle >= 0 &
    result$maintenance$derived_probability_per_8wk_cycle <= 1))
})
