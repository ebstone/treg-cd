# testthat::test_dir() runs with the working directory set to tests/testthat/, but this
# project's scripts source() each other with repo-root-relative paths (README.md convention).
# Load in dependency order via repo-root-relative paths so each script's own
# `if (!exists(...))` guard skips its internal (otherwise-broken-from-here) source() call.
repo_root_relative <- function(...) file.path("..", "..", ...)
source(repo_root_relative("R", "utils", "transition_matrix.R"), local = TRUE)
source(repo_root_relative("R", "00_derive_transition_probs.R"), local = TRUE)
source(repo_root_relative("R", "01_decision_tree.R"), local = TRUE)
source(repo_root_relative("R", "02_markov_engine.R"), local = TRUE)
source(repo_root_relative("R", "03_cure_fraction_module.R"), local = TRUE)

test_that("hazard_to_cycle_probability matches a closed-form check", {
  expect_equal(hazard_to_cycle_probability(0), 0)
  # h = 0.1/year at a 2-week cycle: 1 - exp(-0.1 * 2/52).
  expect_equal(hazard_to_cycle_probability(0.1, cycle_weeks = 2), 1 - exp(-0.1 * 2 / 52), tolerance = 1e-12)
  expect_equal(hazard_to_cycle_probability(0.1, cycle_weeks = 8), 1 - exp(-0.1 * 8 / 52), tolerance = 1e-12)
})

test_that("pi_sdr = 0 exactly reproduces run_maintenance_arm()'s plain output -- the null case falls out for free", {
  states <- MAINTENANCE_STATES
  # Small synthetic matrices are enough here -- this test is about the null case being a genuine
  # no-op, not about UST's specific numbers (that's covered by the real-data test below).
  bio_m <- build_transition_matrix(
    data.frame(from_state = "Remission", to_state = c("Remission", "Mild"), probability = c(0.9, 0.1)),
    states
  )
  ct_m <- build_transition_matrix(
    data.frame(from_state = "Remission", to_state = c("Remission", "Mild"), probability = c(0.8, 0.2)),
    states
  )
  init_bio <- stats::setNames(rep(0, length(states)), states); init_bio["Remission"] <- 1
  init_ct <- stats::setNames(rep(0, length(states)), states)

  treg <- run_treg_arm(bio_m, ct_m, init_bio, init_ct, n_cycles = 60,
                        pi_sdr = 0, relapse_hazard_annual = 0.2, landmark_cycle = 28, cap_cycle = 52)
  plain <- run_maintenance_arm(bio_m, ct_m, init_bio, n_cycles = 60, initial_on_ct = init_ct)

  expect_equal(treg$on_biologic, plain$on_biologic)
  expect_equal(treg$on_ct, plain$on_ct)
  expect_true(all(treg$on_sdr == 0))
})

test_that("horizon at or before the landmark returns the plain trace with on_sdr all zero (cure never manifests)", {
  states <- MAINTENANCE_STATES
  bio_m <- build_transition_matrix(data.frame(from_state = "Remission", to_state = "Remission", probability = 1), states)
  ct_m <- build_transition_matrix(data.frame(from_state = "Remission", to_state = "Remission", probability = 1), states)
  init_bio <- stats::setNames(rep(0, length(states)), states); init_bio["Remission"] <- 1
  init_ct <- stats::setNames(rep(0, length(states)), states)

  treg <- run_treg_arm(bio_m, ct_m, init_bio, init_ct, n_cycles = 10,
                        pi_sdr = 0.9, relapse_hazard_annual = 0.2, landmark_cycle = 28, cap_cycle = 52)
  expect_true(all(treg$on_sdr == 0))
  expect_equal(nrow(treg$on_biologic), 11)
})

test_that("landmark split, SDR decay, and cap-boundary relapse routing (hand-computed)", {
  # Both M-S and Remission absorbing under both matrices -- isolates the split/relapse/cap
  # mechanics from ordinary Markov drift, so every number is hand-traceable.
  states <- c("Moderate-Severe", "Remission")
  bio_m <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE, dimnames = list(states, states))
  ct_m <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE, dimnames = list(states, states))
  init_bio <- c(`Moderate-Severe` = 0, Remission = 1)
  init_ct <- c(`Moderate-Severe` = 0, Remission = 0)

  # landmark_cycle=1, cap_cycle=4: relapses at cycles 2-3 are pre-cap (-> biologic, then
  # immediately switched to CT the *following* cycle by the ordinary Moderate-Severe switch,
  # same as any other biologic-track Moderate-Severe mass); relapses at cycle 4 onward are
  # post-cap (-> CT directly, same cycle).
  h <- 3 # deliberately large so relapse mass is easy to see against tolerance
  res <- run_treg_arm(bio_m, ct_m, init_bio, init_ct, n_cycles = 6, pi_sdr = 1,
                       relapse_hazard_annual = h, landmark_cycle = 1, cap_cycle = 4,
                       apply_cap = TRUE, relapse_destination = "Moderate-Severe")

  p <- hazard_to_cycle_probability(h) # per-cycle relapse probability

  # Landmark (row 2 = cycle 1): all Remission mass (1) moves to SDR, since pi_sdr = 1.
  expect_equal(res$on_sdr[2], 1)
  expect_equal(unname(res$on_biologic[2, "Remission"]), 0)

  # Cycle 2 (row 3): relapse = on_sdr[cycle1] * p, lands in on_biologic (2 < cap_cycle=4).
  relapse_c2 <- 1 * p
  expect_equal(res$on_sdr[3], 1 - relapse_c2, tolerance = 1e-12)
  expect_equal(unname(res$on_biologic[3, "Moderate-Severe"]), relapse_c2, tolerance = 1e-12)
  expect_equal(unname(res$on_ct[3, "Moderate-Severe"]), 0, tolerance = 1e-12)

  # Cycle 3 (row 4): the cycle-2 relapse mass switches from biologic to CT (ordinary M-S switch,
  # every cycle); a fresh relapse lands in on_biologic again (3 < 4).
  relapse_c3 <- res$on_sdr[3] * p
  expect_equal(unname(res$on_ct[4, "Moderate-Severe"]), relapse_c2, tolerance = 1e-12)
  expect_equal(unname(res$on_biologic[4, "Moderate-Severe"]), relapse_c3, tolerance = 1e-12)

  # Cycle 4 == cap_cycle (row 5): the cap sweeps the cycle-3 relapse mass (still on the
  # biologic) to CT; a fresh relapse this cycle lands directly in CT too (4 >= cap_cycle).
  relapse_c4 <- res$on_sdr[4] * p
  expect_equal(unname(res$on_biologic[5, "Moderate-Severe"]), 0, tolerance = 1e-12)
  expect_equal(unname(res$on_ct[5, "Moderate-Severe"]), relapse_c2 + relapse_c3 + relapse_c4, tolerance = 1e-12)

  # Cohort conservation throughout.
  expect_equal(res$total, rep(1, 7), tolerance = 1e-12)
})

test_that("end-to-end: real UST/CT matrices over a full lifetime horizon conserve cohort and exercise post-cap relapse routing", {
  maintenance <- load_published_maintenance(repo_root_relative("data", "raw"))
  ust_m <- build_transition_matrix(maintenance[maintenance$therapy == "UST", ], MAINTENANCE_STATES)
  ct_m <- build_transition_matrix(maintenance[maintenance$therapy == "CT", ], MAINTENANCE_STATES)
  induction <- load_published_induction(repo_root_relative("data", "raw"))
  split <- run_decision_tree(induction[induction$therapy == "UST", ])

  for (pi_sdr in c(0.25, 0.75)) {
    res <- run_treg_arm(
      ust_m, ct_m, split$initial_on_biologic, split$initial_on_ct, n_cycles = 1300,
      pi_sdr = pi_sdr, relapse_hazard_annual = 0.05, landmark_cycle = 28, cap_cycle = 52,
      apply_cap = TRUE
    )
    expect_equal(res$total, rep(1, 1301), tolerance = 1e-9, info = paste("pi_sdr =", pi_sdr))
    # A non-trivial relapse hazard over a 1300-cycle (50-year) horizon should produce some
    # relapse mass that lands after the cap -- confirms the post-cap routing path is actually
    # exercised here, not just reachable in principle.
    expect_true(res$on_ct[53, "Mild"] > 0, info = paste("pi_sdr =", pi_sdr))
    # SDR mass should be strictly decaying (monotonic) after the landmark.
    post_landmark_sdr <- res$on_sdr[29:1301]
    expect_true(all(diff(post_landmark_sdr) <= 1e-12), info = paste("pi_sdr =", pi_sdr))
  }
})
