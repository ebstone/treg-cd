source(file.path("..", "..", "R", "utils", "transition_matrix.R"), local = TRUE)

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

# ---- age_adjust_matrix (lifetime horizon, 2026-08-05) -----------------------------------------

test_that("age_adjust_matrix replaces Death probability and rescales the rest to still sum to 1", {
  m <- matrix(
    c(0.5, 0.3, 0.2, 0, 1, 0, 0, 0, 1),
    nrow = 3, byrow = TRUE, dimnames = list(c("A", "B", "Death"), c("A", "B", "Death"))
  )
  adj <- age_adjust_matrix(m, 0.01)
  expect_equal(unname(rowSums(adj)), rep(1, 3), tolerance = 1e-12)
  expect_equal(unname(adj["A", "Death"]), 0.01)
  # Ratio of A's non-Death entries is preserved (0.5:0.3 = 5:3), only the Death share moves.
  expect_equal(unname(adj["A", "A"] / adj["A", "B"]), 0.5 / 0.3, tolerance = 1e-12)
})

test_that("age_adjust_matrix installs a Death probability into a row that had none (a padded absorbing self-loop)", {
  # B has no Death entry at all (old_death = 0) -- e.g. build_transition_matrix()'s padding for a
  # biologic arm's missing Moderate-Severe row.
  m <- matrix(
    c(0.5, 0.3, 0.2, 0, 1, 0, 0, 0, 1),
    nrow = 3, byrow = TRUE, dimnames = list(c("A", "B", "Death"), c("A", "B", "Death"))
  )
  adj <- age_adjust_matrix(m, 0.02)
  expect_equal(unname(adj["B", ]), c(0, 0.98, 0.02))
})

test_that("age_adjust_matrix(m, 0) zeroes every row's Death probability and rescales the rest up to still sum to 1", {
  # NOT a no-op even though the target probability is 0 -- a row that already HAD a nonzero
  # Death entry (A: 0.2) has that risk entirely removed and its mass redistributed
  # proportionally to the survival states, same mechanic as any other target probability.
  m <- matrix(
    c(0.5, 0.3, 0.2, 0, 1, 0, 0, 0, 1),
    nrow = 3, byrow = TRUE, dimnames = list(c("A", "B", "Death"), c("A", "B", "Death"))
  )
  adj <- age_adjust_matrix(m, 0)
  expect_equal(unname(adj["A", "Death"]), 0)
  expect_equal(unname(rowSums(adj)), rep(1, 3), tolerance = 1e-12)
  expect_equal(unname(adj["A", "A"] / adj["A", "B"]), 0.5 / 0.3, tolerance = 1e-12)
})

test_that("age_adjust_matrix is a true no-op when death_prob already equals every row's own Death entry", {
  m <- matrix(
    c(0.5, 0.3, 0.2, 0.4, 0.4, 0.2, 0, 0, 1),
    nrow = 3, byrow = TRUE, dimnames = list(c("A", "B", "Death"), c("A", "B", "Death"))
  )
  expect_equal(age_adjust_matrix(m, 0.2), m)
})

test_that("age_adjust_matrix leaves the Death row itself untouched", {
  m <- matrix(
    c(0.5, 0.3, 0.2, 0, 1, 0, 0, 0, 1),
    nrow = 3, byrow = TRUE, dimnames = list(c("A", "B", "Death"), c("A", "B", "Death"))
  )
  adj <- age_adjust_matrix(m, 0.5)
  expect_equal(unname(adj["Death", ]), c(0, 0, 1))
})

test_that("age_adjust_matrix leaves a row whose Death probability is already 1 alone (no division by zero)", {
  # C is not the "Death" row itself, but is already 100% bound for it (old_death = 1) -- the
  # denominator (1 - old_death) in age_adjust_matrix()'s rescale formula would be 0 here.
  m <- matrix(
    c(0.5, 0.3, 0.2, 0, 0, 1, 0, 0, 1),
    nrow = 3, byrow = TRUE, dimnames = list(c("A", "C", "Death"), c("A", "C", "Death"))
  )
  adj <- age_adjust_matrix(m, 0.1)
  expect_equal(unname(adj["C", ]), c(0, 0, 1))
})
