# testthat::test_dir() runs with the working directory set to tests/testthat/, but this
# project's scripts source() each other with repo-root-relative paths (README.md convention).
# Same loading pattern as test-markov-engine.R.
repo_root_relative <- function(...) file.path("..", "..", ...)
source(repo_root_relative("R", "utils", "transition_matrix.R"), local = TRUE)
source(repo_root_relative("R", "00_derive_transition_probs.R"), local = TRUE)
source(repo_root_relative("R", "02_markov_engine.R"), local = TRUE)
source(repo_root_relative("R", "03_cure_fraction_module.R"), local = TRUE)
source(repo_root_relative("R", "04_costs_utilities.R"), local = TRUE)

test_that("health_state_monitoring_costs reproduces the cited 2017->2025 inflation exactly", {
  costs <- health_state_monitoring_costs()
  expect_equal(unname(costs["Moderate-Severe"]), 217 * 1.3035)
  expect_equal(unname(costs["Moderate-Severe Responder"]), 217 * 1.3035)
  expect_equal(unname(costs["Mild"]), 91 * 1.3035)
  expect_equal(unname(costs["Remission"]), 10 * 1.3035)
  expect_equal(unname(costs["Death"]), 0)
  # Cross-check against the A9 fix already committed to data/processed/model_health_state_costs.csv
  # (data/data_dictionary.md, resolved 2026-08-04) -- this module computes Surgery's cost
  # independently, from Aliyev's own $884 figure, and the two must agree.
  expect_equal(unname(costs["Surgery"]), 1152.29, tolerance = 1e-2)
})

test_that("ct_drug_cost is Aliyev's $67 (2017) inflated by the same factor as every other line", {
  expect_equal(ct_drug_cost(), 67 * 1.3035)
})

test_that("load_health_state_utilities maps the workbook's abbreviated state names onto MAINTENANCE_STATES", {
  utilities <- load_health_state_utilities(repo_root_relative("data", "processed"))
  expect_setequal(names(utilities), MAINTENANCE_STATES)
  expect_equal(unname(utilities["Death"]), 0)
  # Remission should be the best health state, Death/Surgery among the worst -- a sanity check on
  # the mapping direction (catches an accidental name swap), not a numeric assertion.
  expect_gt(utilities[["Remission"]], utilities[["Mild"]])
  expect_gt(utilities[["Mild"]], utilities[["Moderate-Severe"]])
})

test_that("trace_qalys reproduces a hand-computed 2-state, 2-cycle example", {
  states <- c("A", "B")
  trace <- matrix(c(1, 0, 0.5, 0.5, 0.25, 0.75), nrow = 3, byrow = TRUE, dimnames = list(NULL, states))
  utilities <- c(A = 1, B = 0.5)

  qalys <- trace_qalys(trace, utilities, cycle_weeks = 2, annual_rate = 0)
  # 0% discount rate -> discount_factor is always 1, so this is just utility-mass x cycle-years.
  expected_utility_mass <- c(1, 0.5 * 1 + 0.5 * 0.5, 0.25 * 1 + 0.75 * 0.5)
  expect_equal(qalys, expected_utility_mass * 2 / 52, tolerance = 1e-9)
})

test_that("trace_costs applies the CT drug cost only when add_ct_drug_cost = TRUE", {
  states <- c("Mild", "Death")
  trace <- matrix(c(1, 0, 1, 0), nrow = 2, byrow = TRUE, dimnames = list(NULL, states))
  monitoring_costs <- c(Mild = 100, Death = 0)

  no_drug <- trace_costs(trace, monitoring_costs, annual_rate = 0, add_ct_drug_cost = FALSE)
  with_drug <- trace_costs(trace, monitoring_costs, annual_rate = 0, add_ct_drug_cost = TRUE)

  expect_equal(no_drug, c(100, 100))
  expect_equal(with_drug, c(100, 100) + ct_drug_cost())
})

test_that("attach_sdr_costs_utilities halves the Remission monitoring cost after the cap boundary", {
  on_sdr <- c(1, 1, 1)  # cycles 0, 1, 2; halve_after_cycle = 1 -> cycle 2 is halved
  utilities <- c(Remission = 0.9)
  monitoring_costs <- c(Remission = 20)

  res <- attach_sdr_costs_utilities(on_sdr, utilities, monitoring_costs, annual_rate = 0,
                                     halve_after_cycle = 1)

  expect_equal(res$non_drug_cost_by_cycle, c(20, 20, 10))
  expect_equal(res$qalys_by_cycle, rep(0.9 * 2 / 52, 3))
  # No drug cost is ever added for SDR -- nothing in this module's cost path references a drug
  # cost for SDR, which this test documents as intentional rather than assuming.
})

test_that("attach_maintenance_costs_utilities: CT-track occupancy gets the drug cost, biologic-track doesn't", {
  states <- c("Remission", "Death")
  on_biologic <- matrix(c(1, 0, 1, 0), nrow = 2, byrow = TRUE, dimnames = list(NULL, states))
  on_ct <- matrix(c(0, 0, 0, 0), nrow = 2, byrow = TRUE, dimnames = list(NULL, states))
  arm_result <- list(on_biologic = on_biologic, on_ct = on_ct)
  utilities <- c(Remission = 0.9, Death = 0)
  monitoring_costs <- c(Remission = 20, Death = 0)

  res <- attach_maintenance_costs_utilities(arm_result, utilities, monitoring_costs, annual_rate = 0)
  expect_equal(res$non_drug_cost_by_cycle, c(20, 20))  # no CT drug cost: everyone's on_biologic

  # Swap: everyone on CT instead -- now the drug cost should show up.
  arm_result2 <- list(on_biologic = on_ct, on_ct = on_biologic)
  res2 <- attach_maintenance_costs_utilities(arm_result2, utilities, monitoring_costs, annual_rate = 0)
  expect_equal(res2$non_drug_cost_by_cycle, c(20, 20) + ct_drug_cost())
})

test_that("attach_treg_costs_utilities integrates cleanly against a real, small run_treg_arm() output", {
  maintenance <- load_published_maintenance(repo_root_relative("data", "raw"))
  ust_m <- build_transition_matrix(maintenance[maintenance$therapy == "UST", ], MAINTENANCE_STATES)
  ct_m <- build_transition_matrix(maintenance[maintenance$therapy == "CT", ], MAINTENANCE_STATES)
  init <- stats::setNames(rep(0, length(MAINTENANCE_STATES)), MAINTENANCE_STATES)
  init["Remission"] <- 1

  arm <- run_treg_arm(ust_m, ct_m, init, rep(0, length(MAINTENANCE_STATES)), n_cycles = 60,
                       pi_sdr = 0.3, relapse_hazard_annual = 0.05, landmark_cycle = 28, cap_cycle = 52)

  utilities <- load_health_state_utilities(repo_root_relative("data", "processed"))
  monitoring_costs <- health_state_monitoring_costs()
  res <- attach_treg_costs_utilities(arm, utilities, monitoring_costs)

  # No NAs/negative values, and lifetime QALYs for a cohort of size 1 over 60 cycles (~2.3 years)
  # should be well under the max possible (1 QALY/year x 2.3 years), consistent with occupancy
  # summing to <=1 throughout and utilities <=1.
  expect_false(anyNA(res$qalys_by_cycle))
  expect_false(anyNA(res$non_drug_cost_by_cycle))
  expect_true(all(res$non_drug_cost_by_cycle >= 0))
  expect_true(sum(res$qalys_by_cycle) < 60 * 2 / 52)
})
