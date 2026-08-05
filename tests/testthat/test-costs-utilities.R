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

# ---- UST/IFX drug-cost layer (added 2026-08-04) ------------------------------------------------

test_that("load_dosing_schedule and load_drug_prices read the real sourced files without error", {
  schedule <- load_dosing_schedule(repo_root_relative("data", "raw"))
  prices <- load_drug_prices(repo_root_relative("data", "raw"))
  expect_setequal(schedule$therapy, c("UST", "IFX", "ADA"))
  expect_equal(prices$ust_induction_usd_per_mg, 12.808)
  expect_equal(prices$ust_maintenance_usd_per_mg, 155.883)
  expect_equal(prices$ifx_usd_per_mg, 3.109)
  expect_equal(prices$ada_usd_per_mg, 84.1678)
  expect_equal(prices$iv_administration_usd, 57.9)
})

test_that("ADA's sourced dosing schedule is correct", {
  schedule <- load_dosing_schedule(repo_root_relative("data", "raw"))
  ada <- schedule[schedule$therapy == "ADA", ]

  induction <- ada[ada$phase == "Induction", ]
  expect_equal(nrow(induction), 2)  # two UNEQUAL sequential doses, not one repeated dose like IFX
  expect_setequal(induction$dose_amount, c(160, 80))
  expect_equal(sum(induction$dose_amount * induction$n_doses), 240)  # total induction mg

  maintenance <- ada[ada$phase == "Maintenance", ]
  expect_equal(nrow(maintenance), 1)
  expect_equal(maintenance$dose_amount, 40)
  # Every-other-week = every single native 2-week cycle, NOT every 4th like UST/IFX's q8w --
  # this is the one most likely to silently break if ADA ever reuses UST/IFX's cadence logic.
  expect_equal(maintenance$maintenance_interval_cycles_native_2wk, 1)
  expect_equal(maintenance$first_maintenance_cycle_native_2wk, 2)

  # Unlike UST (induction) and IFX (all doses), ADA is entirely SC -- no row should be IV.
  expect_true(all(ada$route == "SC"))
  # Unlike UST induction (weight-tiered) and IFX (mg/kg), ADA's doses are flat/fixed.
  expect_true(all(is.na(ada$weight_tier_min_kg)), all(is.na(ada$weight_tier_max_kg)))
  expect_equal(unique(ada$dose_unit), "mg")  # not mg/kg
})

test_that("ust_induction_dose_mg picks the correct FDA weight tier, including exact boundaries", {
  schedule <- load_dosing_schedule(repo_root_relative("data", "raw"))
  expect_equal(ust_induction_dose_mg(50, schedule), 260)
  expect_equal(ust_induction_dose_mg(55, schedule), 260)   # boundary: "<=55kg" tier
  expect_equal(ust_induction_dose_mg(55.01, schedule), 390) # just over -> next tier
  expect_equal(ust_induction_dose_mg(71, schedule), 390)    # this project's assumed weight
  expect_equal(ust_induction_dose_mg(85, schedule), 390)    # boundary: "<=85kg" tier
  expect_equal(ust_induction_dose_mg(85.01, schedule), 520)
  expect_equal(ust_induction_dose_mg(150, schedule), 520)
})

test_that("ifx_dose_mg rounds up to the nearest whole 100mg vial, never down", {
  schedule <- load_dosing_schedule(repo_root_relative("data", "raw"))
  expect_equal(ifx_dose_mg(20, schedule), 100)   # 20*5=100, exact -> no waste
  expect_equal(ifx_dose_mg(71, schedule), 400)   # 71*5=355 -> rounds up to 400, not down to 300
  expect_equal(ifx_dose_mg(21, schedule), 200)   # 21*5=105 -> rounds up to 200
})

test_that("induction_drug_cost matches a hand-computed UST and IFX example at the assumed weight", {
  schedule <- load_dosing_schedule(repo_root_relative("data", "raw"))
  prices <- load_drug_prices(repo_root_relative("data", "raw"))

  ust_cost <- induction_drug_cost("UST", ASSUMED_PATIENT_WEIGHT_KG, schedule, prices)
  expect_equal(ust_cost, 390 * prices$ust_induction_usd_per_mg + prices$iv_administration_usd)

  ifx_cost <- induction_drug_cost("IFX", ASSUMED_PATIENT_WEIGHT_KG, schedule, prices)
  # 71kg -> 400mg per dose (vial-rounded), 3 doses, each with its own administration fee.
  expect_equal(ifx_cost, 3 * (400 * prices$ifx_usd_per_mg + prices$iv_administration_usd))
})

test_that("induction_drug_cost sums ADA's two UNEQUAL sequential doses, weight-independent, no admin fee", {
  schedule <- load_dosing_schedule(repo_root_relative("data", "raw"))
  prices <- load_drug_prices(repo_root_relative("data", "raw"))

  ada_cost <- induction_drug_cost("ADA", ASSUMED_PATIENT_WEIGHT_KG, schedule, prices)
  # 160mg + 80mg = 240mg total, at ADA's single NADAC price, no IV administration fee (SC).
  expect_equal(ada_cost, 240 * prices$ada_usd_per_mg)

  # Weight-independent: a very different weight must give the identical cost (unlike UST/IFX).
  ada_cost_other_weight <- induction_drug_cost("ADA", 45, schedule, prices)
  expect_equal(ada_cost, ada_cost_other_weight)
})

test_that("maintenance_drug_cost_by_cycle charges only on dose cycles (4, 8, 12, ...), UST has no admin fee", {
  states <- MAINTENANCE_STATES
  on_biologic <- matrix(0, nrow = 13, ncol = length(states), dimnames = list(NULL, states))
  on_biologic[, "Remission"] <- 1  # everyone on-biologic, in Remission, for all 13 cycles (0..12)

  schedule <- load_dosing_schedule(repo_root_relative("data", "raw"))
  prices <- load_drug_prices(repo_root_relative("data", "raw"))

  res <- maintenance_drug_cost_by_cycle(on_biologic, "UST", ASSUMED_PATIENT_WEIGHT_KG, schedule, prices)
  expected_per_dose <- 90 * prices$ust_maintenance_usd_per_mg  # SC, no admin fee
  is_dose_cycle <- (0:12) %in% c(4, 8, 12)
  expect_equal(res, ifelse(is_dose_cycle, expected_per_dose, 0))

  res_ifx <- maintenance_drug_cost_by_cycle(on_biologic, "IFX", ASSUMED_PATIENT_WEIGHT_KG, schedule, prices)
  expected_per_dose_ifx <- 400 * prices$ifx_usd_per_mg + prices$iv_administration_usd
  expect_equal(res_ifx, ifelse(is_dose_cycle, expected_per_dose_ifx, 0))
})

test_that("maintenance_drug_cost_by_cycle: ADA charges every SINGLE cycle (EOW), not every 4th like UST/IFX", {
  states <- MAINTENANCE_STATES
  on_biologic <- matrix(0, nrow = 9, ncol = length(states), dimnames = list(NULL, states))
  on_biologic[, "Remission"] <- 1  # everyone on-biologic for all 9 cycles (0..8)

  schedule <- load_dosing_schedule(repo_root_relative("data", "raw"))
  prices <- load_drug_prices(repo_root_relative("data", "raw"))

  res <- maintenance_drug_cost_by_cycle(on_biologic, "ADA", ASSUMED_PATIENT_WEIGHT_KG, schedule, prices)
  expected_per_dose <- 40 * prices$ada_usd_per_mg  # SC, no admin fee
  # first_maintenance_cycle_native_2wk = 2, interval = 1 -> every cycle from 2 onward, not
  # every 4th (this is the cadence trap the module header explicitly warns about).
  is_dose_cycle <- (0:8) >= 2
  expect_equal(res, ifelse(is_dose_cycle, expected_per_dose, 0))
})

test_that("maintenance_drug_cost_by_cycle stops charging once on_biologic mass has left (cap/switch)", {
  states <- MAINTENANCE_STATES
  on_biologic <- matrix(0, nrow = 9, ncol = length(states), dimnames = list(NULL, states))
  on_biologic[1:5, "Remission"] <- 1  # cycles 0-4 on biologic
  on_biologic[6:9, "Remission"] <- 0  # cycles 5-8: switched away (cap or CT-switch), mass = 0

  schedule <- load_dosing_schedule(repo_root_relative("data", "raw"))
  prices <- load_drug_prices(repo_root_relative("data", "raw"))
  res <- maintenance_drug_cost_by_cycle(on_biologic, "UST", ASSUMED_PATIENT_WEIGHT_KG, schedule, prices)

  # Only cycle 4 is both a dose cycle AND has mass still on the biologic; cycle 8 would also be a
  # dose cycle but mass is already 0 there.
  expect_equal(res, c(0, 0, 0, 0, 90 * prices$ust_maintenance_usd_per_mg, 0, 0, 0, 0))
})

test_that("attach_maintenance_costs_utilities: passing therapy activates drug_cost_by_cycle correctly", {
  states <- c("Remission", "Death")
  n <- 9
  on_biologic <- matrix(0, nrow = n, ncol = 2, dimnames = list(NULL, states))
  on_biologic[, "Remission"] <- 1
  on_ct <- matrix(0, nrow = n, ncol = 2, dimnames = list(NULL, states))
  arm_result <- list(on_biologic = on_biologic, on_ct = on_ct)
  utilities <- c(Remission = 0.9, Death = 0)
  monitoring_costs <- c(Remission = 20, Death = 0)

  # therapy = NULL (default): drug_cost_by_cycle is all zero, total = non_drug.
  res_none <- attach_maintenance_costs_utilities(arm_result, utilities, monitoring_costs, annual_rate = 0)
  expect_equal(res_none$drug_cost_by_cycle, rep(0, n))
  expect_equal(res_none$total_cost_by_cycle, res_none$non_drug_cost_by_cycle)

  # therapy = "UST": drug_cost_by_cycle picks up the maintenance dosing schedule.
  # Only need a fake schedule/prices here to keep the test hand-computable and independent of the
  # real sourced files' exact numbers (those are already checked by the loader tests above).
  fake_schedule <- data.frame(
    therapy = "UST", phase = "Maintenance", dose_amount = 90, route = "SC",
    first_maintenance_cycle_native_2wk = 2, maintenance_interval_cycles_native_2wk = 2,
    stringsAsFactors = FALSE
  )
  fake_prices <- list(ust_maintenance_usd_per_mg = 1)  # $1/mg -> per-dose cost = $90, easy to check
  res_ust <- attach_maintenance_costs_utilities(
    arm_result, utilities, monitoring_costs, annual_rate = 0,
    therapy = "UST", weight_kg = 71, schedule = fake_schedule, prices = fake_prices
  )
  is_dose_cycle <- (0:(n - 1)) %in% c(2, 4, 6, 8)
  expect_equal(res_ust$drug_cost_by_cycle, ifelse(is_dose_cycle, 90, 0))
  expect_equal(res_ust$total_cost_by_cycle, res_ust$non_drug_cost_by_cycle + res_ust$drug_cost_by_cycle)
})

test_that("attach_treg_costs_utilities never charges UST's drug cost to the non-cured track", {
  # Real integration check: even though Treg's non-cured portion literally reuses UST's own
  # maintenance matrix (R/03_cure_fraction_module.R), attach_treg_costs_utilities() must never
  # produce a nonzero drug_cost_by_cycle -- confirmed with E. Stone, 2026-08-04 (module header).
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

  expect_equal(res$drug_cost_by_cycle, rep(0, 61))
  expect_equal(res$total_cost_by_cycle, res$non_drug_cost_by_cycle)
})

test_that("summarise_arm adds a one-time induction cost undiscounted, on top of the per-cycle sums", {
  attached <- list(
    qalys_by_cycle = c(1, 2, 3),
    non_drug_cost_by_cycle = c(10, 10, 10),
    drug_cost_by_cycle = c(5, 0, 5),
    total_cost_by_cycle = c(15, 10, 15)
  )
  res <- summarise_arm(attached, induction_cost = 1000)
  expect_equal(res$qalys, 6)
  expect_equal(res$non_drug_cost, 30)
  expect_equal(res$drug_cost, 10 + 1000)
  expect_equal(res$total_cost, 40 + 1000)

  res_no_induction <- summarise_arm(attached)
  expect_equal(res_no_induction$drug_cost, 10)
  expect_equal(res_no_induction$total_cost, 40)
})
