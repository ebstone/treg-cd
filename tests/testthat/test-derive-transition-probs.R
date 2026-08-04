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

test_that("validate_row_sums passes a well-formed transition table", {
  df <- data.frame(
    treatment = c("A", "A", "B", "B"),
    from_state = c("X", "X", "X", "X"),
    to_state = c("X", "Y", "X", "Y"),
    probability = c(0.4, 0.6, 0.9, 0.1)
  )
  expect_true(validate_row_sums(df, by = "treatment"))
})

test_that("validate_row_sums catches rows that don't sum to 1", {
  df <- data.frame(
    treatment = c("A", "A"),
    from_state = c("X", "X"),
    to_state = c("X", "Y"),
    probability = c(0.4, 0.4)
  )
  expect_error(validate_row_sums(df, by = "treatment"), "Row sums not equal to 1")
})

test_that("validate_row_sums catches out-of-range probabilities", {
  df <- data.frame(
    treatment = "A", from_state = "X", to_state = c("X", "Y"), probability = c(1.5, -0.5)
  )
  expect_error(validate_row_sums(df, by = "treatment"))
})

test_that("build_transition_matrix pads a missing from_state as an absorbing self-loop", {
  df <- data.frame(
    from_state = c("A", "A"), to_state = c("A", "B"), probability = c(0.3, 0.7),
    stringsAsFactors = FALSE
  )
  m <- build_transition_matrix(df, states = c("A", "B"))
  expect_equal(unname(m["A", ]), c(0.3, 0.7))
  expect_equal(unname(m["B", ]), c(0, 1))
})

test_that("matrix_power(m, 1) is the identity operation", {
  m <- matrix(c(0.5, 0.5, 0.2, 0.8), nrow = 2, byrow = TRUE, dimnames = list(c("A", "B"), c("A", "B")))
  expect_equal(matrix_power(m, 1), m)
})

test_that("matrix_power correctly composes a two-state chain (hand-computable)", {
  # Simple 2-state chain: from A, 0.5 stay / 0.5 go to B; from B, always stay.
  m <- matrix(c(0.5, 0.5, 0, 1), nrow = 2, byrow = TRUE, dimnames = list(c("A", "B"), c("A", "B")))
  m2 <- matrix_power(m, 2)
  # P(A->A after 2 steps) = 0.5 * 0.5 = 0.25; P(A->B after 2 steps) = 1 - 0.25 = 0.75
  expect_equal(m2["A", "A"], 0.25, tolerance = 1e-9)
  expect_equal(m2["A", "B"], 0.75, tolerance = 1e-9)
  expect_equal(m2["B", "B"], 1, tolerance = 1e-9)
})

test_that("matrix_to_long round-trips through build_transition_matrix", {
  m <- matrix(c(0.3, 0.7, 0.1, 0.9), nrow = 2, byrow = TRUE, dimnames = list(c("A", "B"), c("A", "B")))
  long <- matrix_to_long(m, therapy = "TEST")
  expect_setequal(long$therapy, "TEST")
  rebuilt <- build_transition_matrix(
    data.frame(from_state = long$from_state, to_state = long$to_state, probability = long$probability),
    states = c("A", "B")
  )
  expect_equal(rebuilt, m)
})

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

test_that("load_published_maintenance_8wk produces valid, row-sum-1 matrices for all four therapies", {
  out <- load_published_maintenance_8wk(repo_root_relative("data", "raw"))
  expect_setequal(out$therapy, c("UST", "IFX", "ADA", "CT"))
  # Rounding in the published 2-week values compounds slightly across four matrix
  # multiplications (up to ~0.2% off), same allowance run_derivation() uses.
  expect_true(all(out$probability >= -1e-9 & out$probability <= 1 + 1e-9))
  expect_true(validate_row_sums(out, by = "therapy", tol = 0.005))
})

test_that("load_published_maintenance_8wk: Moderate-Severe is absorbing for biologic arms, not for CT", {
  out <- load_published_maintenance_8wk(repo_root_relative("data", "raw"))
  ust_ms_to_ms <- out[out$therapy == "UST" & out$from_state == "Moderate-Severe" & out$to_state == "Moderate-Severe", "probability"]
  expect_equal(ust_ms_to_ms, 1)

  ct_ms_to_ms <- out[out$therapy == "CT" & out$from_state == "Moderate-Severe" & out$to_state == "Moderate-Severe", "probability"]
  expect_true(ct_ms_to_ms < 1)
})

test_that("run_derivation produces deterministic output", {
  first <- run_from_repo_root(write_output = FALSE)
  second <- run_from_repo_root(write_output = FALSE)
  expect_equal(first, second)
})

test_that("maintenance comparison against the workbook snapshot covers UST, IFX, CT (ADA is not in the v6 workbook)", {
  result <- run_from_repo_root(write_output = FALSE)
  expect_setequal(result$maintenance_comparison$therapy, c("UST", "IFX", "CT"))
  expect_true(all(is.finite(result$maintenance_comparison$abs_diff)))
})
