# testthat::test_dir() runs with the working directory set to tests/testthat/, but this
# project's scripts source() each other with repo-root-relative paths (README.md convention).
# Load utils first via a repo-root-relative path so R/00's and R/02's own
# `if (!exists("validate_row_sums"))` guards skip their internal (otherwise-broken-from-here)
# source() calls.
repo_root_relative <- function(...) file.path("..", "..", ...)
source(repo_root_relative("R", "utils", "transition_matrix.R"), local = TRUE)
source(repo_root_relative("R", "utils", "life_table.R"), local = TRUE)
source(repo_root_relative("R", "00_derive_transition_probs.R"), local = TRUE)
source(repo_root_relative("R", "02_markov_engine.R"), local = TRUE)

RAW_DIR <- repo_root_relative("data", "raw")

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

# ---- Lifetime horizon: age- and sex-specific background mortality (2026-08-05) ----------------

test_that("run_maintenance_arm's death_prob_schedule defaults to NULL and is fully backward compatible", {
  maintenance <- load_published_maintenance(RAW_DIR)
  ust_m <- build_transition_matrix(maintenance[maintenance$therapy == "UST", ], MAINTENANCE_STATES)
  ct_m <- build_transition_matrix(maintenance[maintenance$therapy == "CT", ], MAINTENANCE_STATES)
  init <- stats::setNames(rep(0, length(MAINTENANCE_STATES)), MAINTENANCE_STATES)
  init["Remission"] <- 1

  without_arg <- run_maintenance_arm(ust_m, ct_m, init, n_cycles = 50, cap_cycle = 52)
  with_explicit_null <- run_maintenance_arm(ust_m, ct_m, init, n_cycles = 50, cap_cycle = 52,
                                             death_prob_schedule = NULL)
  expect_equal(without_arg, with_explicit_null)
})

test_that("run_maintenance_arm rejects a death_prob_schedule of the wrong length", {
  maintenance <- load_published_maintenance(RAW_DIR)
  ust_m <- build_transition_matrix(maintenance[maintenance$therapy == "UST", ], MAINTENANCE_STATES)
  ct_m <- build_transition_matrix(maintenance[maintenance$therapy == "CT", ], MAINTENANCE_STATES)
  init <- stats::setNames(rep(0, length(MAINTENANCE_STATES)), MAINTENANCE_STATES)
  init["Remission"] <- 1
  expect_error(run_maintenance_arm(ust_m, ct_m, init, n_cycles = 10, death_prob_schedule = c(0.01, 0.02)))
})

test_that("run_maintenance_arm with a death_prob_schedule still conserves total cohort mass at 1", {
  maintenance <- load_published_maintenance(RAW_DIR)
  ust_m <- build_transition_matrix(maintenance[maintenance$therapy == "UST", ], MAINTENANCE_STATES)
  ct_m <- build_transition_matrix(maintenance[maintenance$therapy == "CT", ], MAINTENANCE_STATES)
  init <- stats::setNames(rep(0, length(MAINTENANCE_STATES)), MAINTENANCE_STATES)
  init["Remission"] <- 1

  lt <- load_life_table(RAW_DIR)
  schedule <- death_prob_schedule(80, "male", n_cycles = 200, life_table = lt)  # old cohort -> high mortality
  res <- run_maintenance_arm(ust_m, ct_m, init, n_cycles = 200, death_prob_schedule = schedule)
  expect_equal(rowSums(res$total), rep(1, 201), tolerance = 1e-9)
  # A meaningfully older cohort should end up with substantial Death-state mass, not near-zero.
  expect_true(res$total[201, "Death"] > 0.3)
})

test_that("age-adjusted mortality strictly increases Death occupancy relative to Aliyev's own tiny fixed rate, for an old cohort", {
  maintenance <- load_published_maintenance(RAW_DIR)
  ust_m <- build_transition_matrix(maintenance[maintenance$therapy == "UST", ], MAINTENANCE_STATES)
  ct_m <- build_transition_matrix(maintenance[maintenance$therapy == "CT", ], MAINTENANCE_STATES)
  init <- stats::setNames(rep(0, length(MAINTENANCE_STATES)), MAINTENANCE_STATES)
  init["Remission"] <- 1

  lt <- load_life_table(RAW_DIR)
  schedule <- death_prob_schedule(85, "male", n_cycles = 100, life_table = lt)
  with_mortality <- run_maintenance_arm(ust_m, ct_m, init, n_cycles = 100, death_prob_schedule = schedule)
  without <- run_maintenance_arm(ust_m, ct_m, init, n_cycles = 100)  # Aliyev's own flat 0.00006/cycle

  expect_true(with_mortality$total[101, "Death"] > without$total[101, "Death"])
})

test_that("run_maintenance_arm_with_mortality splits 50/50 by sex and sums to the same total mass as a single-sex run at double weight", {
  maintenance <- load_published_maintenance(RAW_DIR)
  ust_m <- build_transition_matrix(maintenance[maintenance$therapy == "UST", ], MAINTENANCE_STATES)
  ct_m <- build_transition_matrix(maintenance[maintenance$therapy == "CT", ], MAINTENANCE_STATES)
  init <- stats::setNames(rep(0, length(MAINTENANCE_STATES)), MAINTENANCE_STATES)
  init["Remission"] <- 1
  lt <- load_life_table(RAW_DIR)

  combined <- run_maintenance_arm_with_mortality(ust_m, ct_m, init, n_cycles = 60, baseline_age = 60,
                                                  life_table = lt)
  expect_equal(rowSums(combined$total), rep(1, 61), tolerance = 1e-9)

  # Sanity: male-only and female-only runs (each starting the FULL cohort, not half) should
  # differ (male mortality exceeds female at every adult age in this table -- test-life-table.R),
  # and the combined 50/50 run's Death share should land strictly between them.
  male_schedule <- death_prob_schedule(60, "male", n_cycles = 60, life_table = lt)
  female_schedule <- death_prob_schedule(60, "female", n_cycles = 60, life_table = lt)
  male_only <- run_maintenance_arm(ust_m, ct_m, init, n_cycles = 60, death_prob_schedule = male_schedule)
  female_only <- run_maintenance_arm(ust_m, ct_m, init, n_cycles = 60, death_prob_schedule = female_schedule)

  expect_true(male_only$total[61, "Death"] > combined$total[61, "Death"])
  expect_true(combined$total[61, "Death"] > female_only$total[61, "Death"])
  # And the combined run is exactly the average (each sub-cohort weighted 0.5) of the two.
  expect_equal(combined$total[61, "Death"], 0.5 * male_only$total[61, "Death"] + 0.5 * female_only$total[61, "Death"],
               tolerance = 1e-9)
})
