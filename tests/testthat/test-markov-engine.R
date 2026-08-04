# testthat::test_dir() runs with the working directory set to tests/testthat/, but this
# project's scripts source() each other with repo-root-relative paths (README.md convention).
# Load utils first via a repo-root-relative path so R/00's and R/02's own
# `if (!exists("validate_row_sums"))` guards skip their internal (otherwise-broken-from-here)
# source() calls.
repo_root_relative <- function(...) file.path("..", "..", ...)
source(repo_root_relative("R", "utils", "transition_matrix.R"), local = TRUE)
source(repo_root_relative("R", "00_derive_transition_probs.R"), local = TRUE)
source(repo_root_relative("R", "02_markov_engine.R"), local = TRUE)

test_that("transition matrix rows sum to 1", {
  maintenance <- load_published_maintenance(repo_root_relative("data", "raw"))
  for (tx in unique(maintenance$therapy)) {
    m <- build_transition_matrix(maintenance[maintenance$therapy == tx, ], MAINTENANCE_STATES)
    expect_equal(unname(rowSums(m)), rep(1, length(MAINTENANCE_STATES)), tolerance = 1e-9, info = tx)
  }
})

test_that("cohort size is conserved across cycles", {
  maintenance <- load_published_maintenance(repo_root_relative("data", "raw"))
  ust_m <- build_transition_matrix(maintenance[maintenance$therapy == "UST", ], MAINTENANCE_STATES)
  ct_m <- build_transition_matrix(maintenance[maintenance$therapy == "CT", ], MAINTENANCE_STATES)
  init <- stats::setNames(rep(0, length(MAINTENANCE_STATES)), MAINTENANCE_STATES)
  init["Remission"] <- 1

  # Long horizon (300 cycles = ~11.5 years at 2-week cycles) so any per-cycle rounding drift
  # would be visible -- build_transition_matrix()'s row-normalisation (R/utils/transition_matrix.R)
  # exists specifically to prevent this from compounding.
  for (apply_cap in c(TRUE, FALSE)) {
    res <- run_maintenance_arm(ust_m, ct_m, init, n_cycles = 300, cap_cycle = 52, apply_cap = apply_cap)
    totals <- rowSums(res$total)
    expect_equal(totals, rep(1, length(totals)), tolerance = 1e-9, info = paste("apply_cap =", apply_cap))
  }
})

test_that("QALY/cost discounting matches a closed-form check", {
  expect_equal(discount_factor(0), 1)
  # Cycle 26 at a 2-week cycle length = 52 weeks = 1 year.
  expect_equal(discount_factor(26, cycle_weeks = 2, annual_rate = 0.03), 1 / 1.03, tolerance = 1e-9)
  expect_equal(discount_factor(52, cycle_weeks = 2, annual_rate = 0.03), 1 / 1.03^2, tolerance = 1e-9)
  # 0% discount rate is always 1, regardless of elapsed time.
  expect_equal(discount_factor(100, annual_rate = 0), 1)
})

test_that("simulate_cohort reproduces a fixed, hand-computed expected output", {
  states <- c("A", "B")
  # From A: 60% stay, 40% to B. From B: absorbing.
  m <- matrix(c(0.6, 0.4, 0, 1), nrow = 2, byrow = TRUE, dimnames = list(states, states))
  trace <- simulate_cohort(m, initial_state = c(A = 1, B = 0), n_cycles = 3)

  # Hand-computed: P(A) after n cycles = 0.6^n.
  expect_equal(trace[, "A"], c(1, 0.6, 0.36, 0.216), tolerance = 1e-9)
  expect_equal(trace[, "B"], c(0, 0.4, 0.64, 0.784), tolerance = 1e-9)
})

test_that("run_maintenance_arm: Moderate-Severe switches to CT every cycle (hand-computed, no cap)", {
  states <- c("Moderate-Severe", "Remission")
  # Biologic: from M-S, 50/50 stay/go to Remission; Remission absorbing.
  bio_m <- matrix(c(0.5, 0.5, 0, 1), nrow = 2, byrow = TRUE, dimnames = list(states, states))
  # CT: both states absorbing (kept simple -- this test is about the switch mechanics, not CT's
  # own dynamics).
  ct_m <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE, dimnames = list(states, states))
  init <- c(`Moderate-Severe` = 1, Remission = 0)

  res <- run_maintenance_arm(bio_m, ct_m, init, n_cycles = 2, cap_cycle = 52, apply_cap = FALSE)

  # Cycle 1: bio (1,0) -> (0.5,0.5); M-S mass (0.5) exits to CT; bio left holding only the
  # Remission mass it produced this cycle: (0, 0.5). CT: (0.5, 0).
  expect_equal(unname(res$on_biologic[2, ]), c(0, 0.5), tolerance = 1e-9)
  expect_equal(unname(res$on_ct[2, ]), c(0.5, 0), tolerance = 1e-9)
  # Cycle 2: bio's M-S component is already 0, so nothing further to switch -- steady state.
  expect_equal(unname(res$on_biologic[3, ]), c(0, 0.5), tolerance = 1e-9)
  expect_equal(unname(res$on_ct[3, ]), c(0.5, 0), tolerance = 1e-9)
  expect_equal(rowSums(res$total), rep(1, 3), tolerance = 1e-9)
})

test_that("run_maintenance_arm: 2-year cap moves ALL remaining biologic mass to CT at cap_cycle (hand-computed)", {
  states <- c("Moderate-Severe", "Remission")
  bio_m <- matrix(c(0.5, 0.5, 0, 1), nrow = 2, byrow = TRUE, dimnames = list(states, states))
  ct_m <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE, dimnames = list(states, states))
  init <- c(`Moderate-Severe` = 1, Remission = 0)

  res <- run_maintenance_arm(bio_m, ct_m, init, n_cycles = 2, cap_cycle = 1, apply_cap = TRUE)

  # Cycle 1 = cap_cycle: after the M-S switch (bio -> (0, 0.5), ct -> (0.5, 0)), the cap then
  # moves ALL remaining bio mass (the 0.5 in Remission) to CT too: bio -> (0,0), ct -> (0.5, 0.5).
  expect_equal(unname(res$on_biologic[2, ]), c(0, 0), tolerance = 1e-9)
  expect_equal(unname(res$on_ct[2, ]), c(0.5, 0.5), tolerance = 1e-9)
  # Cycle 2: biologic track stays at zero for the rest of the horizon.
  expect_equal(unname(res$on_biologic[3, ]), c(0, 0), tolerance = 1e-9)
  expect_equal(unname(res$on_ct[3, ]), c(0.5, 0.5), tolerance = 1e-9)
})
