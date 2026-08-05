# Tests for R/09_scenarios.R -- structural scenarios S1, S2, S4, S5, S10, S11 (analysis_plan.md
# §10.3), the ones needing no new sourcing (S3's own tests live in test-refractory-multipliers.R).

repo_root_relative <- function(...) file.path("..", "..", ...)
source(repo_root_relative("R", "utils", "transition_matrix.R"), local = TRUE)
source(repo_root_relative("R", "utils", "life_table.R"), local = TRUE)
source(repo_root_relative("R", "utils", "refractory_multipliers.R"), local = TRUE)
source(repo_root_relative("R", "00_derive_transition_probs.R"), local = TRUE)
source(repo_root_relative("R", "01_decision_tree.R"), local = TRUE)
source(repo_root_relative("R", "02_markov_engine.R"), local = TRUE)
source(repo_root_relative("R", "03_cure_fraction_module.R"), local = TRUE)
source(repo_root_relative("R", "04_costs_utilities.R"), local = TRUE)
source(repo_root_relative("R", "05_deterministic_results.R"), local = TRUE)
source(repo_root_relative("R", "09_scenarios.R"), local = TRUE)

RAW_DIR <- repo_root_relative("data", "raw")
PROC_DIR <- repo_root_relative("data", "processed")

# ---- S1: 2-year cap on vs. off -----------------------------------------------------------------

test_that("run_scenario_s1_cap's 'on' rows exactly reproduce run_base_case()'s own default", {
  scenario <- run_scenario_s1_cap(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  on_rows <- scenario[scenario$cap == "on", setdiff(names(scenario), "cap")]
  base <- run_base_case(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  rownames(on_rows) <- NULL
  rownames(base) <- NULL
  expect_equal(on_rows, base)
})

test_that("run_scenario_s1_cap's 'off' rows differ from 'on' and raise both QALYs and cost for the comparators", {
  scenario <- run_scenario_s1_cap(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  for (tx in COMPARATOR_THERAPIES) {
    on_row <- scenario[scenario$cap == "on" & scenario$intervention == tx, ]
    off_row <- scenario[scenario$cap == "off" & scenario$intervention == tx, ]
    # Removing the 2-year cap keeps patients on the (more effective, more expensive) biologic
    # track longer instead of switching them to CT -- both QALYs and cost rise.
    expect_true(off_row$qalys > on_row$qalys)
    expect_true(off_row$total_cost > on_row$total_cost)
  }
})

# ---- S2: ADA included vs. excluded -------------------------------------------------------------

test_that("run_scenario_s2_comparator_set drops ADA's row but leaves UST/IFX/TREG unchanged", {
  scenario <- run_scenario_s2_comparator_set(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  with_ada <- scenario[scenario$comparator_set == "UST_IFX_ADA", ]
  without_ada <- scenario[scenario$comparator_set == "UST_IFX", ]
  expect_setequal(with_ada$intervention, c("UST", "IFX", "ADA", "TREG"))
  expect_setequal(without_ada$intervention, c("UST", "IFX", "TREG"))
  shared <- c("UST", "IFX", "TREG")
  cols <- setdiff(names(scenario), "comparator_set")
  a <- with_ada[with_ada$intervention %in% shared, cols]
  b <- without_ada[without_ada$intervention %in% shared, cols]
  a <- a[order(a$intervention), ]
  b <- b[order(b$intervention), ]
  rownames(a) <- NULL; rownames(b) <- NULL
  expect_equal(a, b)
})

test_that("run_scenario_s2_headroom_at_sourced_price returns valid pi* rows for both comparator sets", {
  res <- run_scenario_s2_headroom_at_sourced_price(
    wtp_grid_usd = 150000, n_cycles = HORIZON_CYCLES_LIFETIME,
    baseline_age = ASSUMED_PATIENT_AGE_YEARS, raw_dir = RAW_DIR, proc_dir = PROC_DIR
  )
  expect_equal(nrow(res), 2)
  expect_setequal(res$comparator_set, c("UST_IFX_ADA", "UST_IFX"))
  expect_true(all(res$feasible))
  expect_true(all(res$pi_star >= 0 & res$pi_star <= 1))
})

# ---- S4: cure-fraction anchor table ------------------------------------------------------------

test_that("run_cure_fraction_anchor_table covers the full pi x duration grid with sane, monotone values", {
  res <- run_cure_fraction_anchor_table(pi_grid = c(0, 0.5, 0.9), duration_grid_years = c(2, 10),
                                        raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_equal(nrow(res), 3 * 2)
  expect_false(anyNA(res))
  # Higher pi (more durable cures) can only raise QALYs and lower cost, at fixed duration --
  # same monotonicity headroom_pi_star() itself already relies on (R/05's own docstring).
  for (d in c(2, 10)) {
    sub <- res[res$duration_years == d, ]
    sub <- sub[order(sub$pi_sdr), ]
    expect_true(all(diff(sub$qalys) >= -1e-9))
    expect_true(all(diff(sub$total_cost) <= 1e-6))
  }
  # At pi = 0, duration/relapse hazard cannot matter (no one is ever in SDR to relapse).
  pi0 <- res[res$pi_sdr == 0, ]
  expect_equal(pi0$qalys[1], pi0$qalys[2])
  expect_equal(pi0$total_cost[1], pi0$total_cost[2])
})

# ---- S5: time horizon ---------------------------------------------------------------------------

test_that("run_scenario_s5_horizon's 6.15-year rows reproduce run_base_case()'s own default, and QALYs/cost rise with horizon length", {
  scenario <- run_scenario_s5_horizon(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_setequal(scenario$horizon, c("6.15_year", "10_year", "lifetime"))

  six_year_rows <- scenario[scenario$horizon == "6.15_year", setdiff(names(scenario), "horizon")]
  base <- run_base_case(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  rownames(six_year_rows) <- NULL; rownames(base) <- NULL
  expect_equal(six_year_rows, base)

  for (tx in c(COMPARATOR_THERAPIES, "TREG")) {
    q6 <- scenario$qalys[scenario$horizon == "6.15_year" & scenario$intervention == tx]
    q10 <- scenario$qalys[scenario$horizon == "10_year" & scenario$intervention == tx]
    qlife <- scenario$qalys[scenario$horizon == "lifetime" & scenario$intervention == tx]
    expect_true(q10 > q6)
    expect_true(qlife > q10)
  }
})

# ---- S10: discount rate --------------------------------------------------------------------------

test_that("run_scenario_s10_discount_rate's 3% rows reproduce run_base_case()'s own default, and lower rates raise both QALYs and cost", {
  scenario <- run_scenario_s10_discount_rate(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_setequal(scenario$annual_discount_rate, DISCOUNT_RATE_GRID)

  base_rate_rows <- scenario[scenario$annual_discount_rate == 0.03,
                              setdiff(names(scenario), "annual_discount_rate")]
  base <- run_base_case(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  rownames(base_rate_rows) <- NULL; rownames(base) <- NULL
  expect_equal(base_rate_rows, base)

  for (tx in c(COMPARATOR_THERAPIES, "TREG")) {
    by_rate <- scenario[scenario$intervention == tx, ]
    by_rate <- by_rate[order(by_rate$annual_discount_rate), ]
    # Lower discounting weights future cycles more heavily -- QALYs and cost both rise as the
    # rate falls (0% > 1.5% > 3% > 5%), matching every other annual_rate call site in this project.
    expect_true(all(diff(by_rate$qalys) <= 1e-9))
    expect_true(all(diff(by_rate$total_cost) <= 1e-6))
  }
})

# ---- S11: SDR relapse destination -----------------------------------------------------------------

test_that("run_scenario_s11_relapse_destination shows no difference at pi=0 and a small, correctly-signed difference at pi>0", {
  res <- run_scenario_s11_relapse_destination(pi_grid = c(0, 0.5, 0.9), raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_equal(nrow(res), 3 * 2)

  mild0 <- res$qalys[res$pi_sdr == 0 & res$relapse_destination == "Mild"]
  msr0 <- res$qalys[res$pi_sdr == 0 & res$relapse_destination == "Moderate-Severe Responder"]
  expect_equal(mild0, msr0)

  for (pi_val in c(0.5, 0.9)) {
    mild <- res$qalys[res$pi_sdr == pi_val & res$relapse_destination == "Mild"]
    msr <- res$qalys[res$pi_sdr == pi_val & res$relapse_destination == "Moderate-Severe Responder"]
    # Mild is a better health state than Moderate-Severe Responder, so relapsing to Mild can only
    # be at least as good -- a small but non-negative QALY gap.
    expect_true(mild >= msr)
  }
})

test_that("run_treg_arm_lifetime's relapse_destination = 'Mild' default is byte-identical to omitting it", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  treg_price <- load_treg_dose_acquisition_cost(PROC_DIR)
  without_arg <- run_treg_arm_lifetime(HORIZON_CYCLES_6YR, pi_sdr = 0.5, relapse_hazard_annual = 0.07,
                                       price_usd = treg_price, matrices = matrices,
                                       raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  explicit_default <- run_treg_arm_lifetime(HORIZON_CYCLES_6YR, pi_sdr = 0.5, relapse_hazard_annual = 0.07,
                                             price_usd = treg_price, matrices = matrices,
                                             raw_dir = RAW_DIR, proc_dir = PROC_DIR,
                                             relapse_destination = "Mild")
  expect_equal(without_arg, explicit_default)
})

# ---- S12: non-cured Treg = UST-equivalent vs. HR-advantaged ------------------------------------

test_that("apply_non_cured_hazard_ratio(m, 1) is the exact identity transform", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  for (tx in c(COMPARATOR_THERAPIES, "CT")) {
    m <- matrices[[tx]]
    expect_identical(apply_non_cured_hazard_ratio(m, 1), m)
  }
})

test_that("apply_non_cured_hazard_ratio preserves row-stochasticity and only ever reduces M-S/Surgery, leaving Death untouched", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  for (tx in c(COMPARATOR_THERAPIES, "CT")) {
    m <- matrices[[tx]]
    for (hr in c(0.9, 0.5, 0.1)) {
      adj <- apply_non_cured_hazard_ratio(m, hr)
      expect_equal(unname(rowSums(adj)), rep(1, nrow(adj)), tolerance = 1e-9)
      expect_true(all(adj >= -1e-12 & adj <= 1 + 1e-12))
      expect_true(all(adj[, "Moderate-Severe"] <= m[, "Moderate-Severe"] + 1e-12))
      expect_true(all(adj[, "Surgery"] <= m[, "Surgery"] + 1e-12))
      expect_equal(adj[, "Death"], m[, "Death"])
    }
  }
})

test_that("apply_non_cured_hazard_ratio matches the discrete-time proportional-hazards formula exactly", {
  m <- build_all_transition_matrices(RAW_DIR)[["UST"]]
  hr <- 0.6
  adj <- apply_non_cured_hazard_ratio(m, hr)
  expect_equal(adj[, "Moderate-Severe"], 1 - (1 - m[, "Moderate-Severe"]) ^ hr)
  expect_equal(adj[, "Surgery"], 1 - (1 - m[, "Surgery"]) ^ hr)
})

test_that("apply_non_cured_hazard_ratio rejects a non-positive hazard ratio", {
  m <- build_all_transition_matrices(RAW_DIR)[["UST"]]
  expect_error(apply_non_cured_hazard_ratio(m, 0))
  expect_error(apply_non_cured_hazard_ratio(m, -0.5))
})

test_that("run_treg_arm_lifetime's non_cured_hazard_ratio = 1 default is byte-identical to omitting it", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  treg_price <- load_treg_dose_acquisition_cost(PROC_DIR)
  without_arg <- run_treg_arm_lifetime(HORIZON_CYCLES_6YR, pi_sdr = 0, relapse_hazard_annual = 0,
                                        price_usd = treg_price, matrices = matrices,
                                        raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  explicit_default <- run_treg_arm_lifetime(HORIZON_CYCLES_6YR, pi_sdr = 0, relapse_hazard_annual = 0,
                                             price_usd = treg_price, matrices = matrices,
                                             raw_dir = RAW_DIR, proc_dir = PROC_DIR,
                                             non_cured_hazard_ratio = 1)
  expect_equal(without_arg, explicit_default)
})

test_that("run_treg_arm_lifetime(non_cured_hazard_ratio < 1) never mutates the shared matrices list", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  before <- matrices[["UST"]]
  treg_price <- load_treg_dose_acquisition_cost(PROC_DIR)
  invisible(run_treg_arm_lifetime(HORIZON_CYCLES_6YR, pi_sdr = 0, relapse_hazard_annual = 0,
                                   price_usd = treg_price, matrices = matrices,
                                   raw_dir = RAW_DIR, proc_dir = PROC_DIR,
                                   non_cured_hazard_ratio = 0.5))
  expect_identical(matrices[["UST"]], before)
})

test_that("run_scenario_s12_non_cured_hr's HR=1 row exactly reproduces the identity-transform Treg outcome, and NMB rises monotonically as the advantage strengthens", {
  res <- run_scenario_s12_non_cured_hr(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_equal(nrow(res), length(NON_CURED_HR_GRID))

  matrices <- build_all_transition_matrices(RAW_DIR)
  treg_price <- load_treg_dose_acquisition_cost(PROC_DIR)
  base <- run_treg_arm_lifetime(HORIZON_CYCLES_6YR, pi_sdr = 0, relapse_hazard_annual = 0,
                                 price_usd = treg_price, matrices = matrices,
                                 raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  hr1_row <- res[res$non_cured_hazard_ratio == 1, ]
  expect_equal(hr1_row$qalys, base$qalys)
  expect_equal(hr1_row$total_cost, base$total_cost)

  # Cost falls monotonically as the advantage strengthens (Surgery/M-S are the expensive states);
  # NMB rises monotonically at every WTP threshold, dominated by that cost saving -- verified
  # (R/09_scenarios.R's own module header) that QALYs alone are NOT guaranteed monotone here,
  # a real, explained structural feature of Aliyev's own matrix (Surgery's own row has a fast
  # one-step return to Remission), not a bug -- so this test deliberately does NOT assert
  # monotone QALYs, only monotone cost and NMB, the metric that's actually decision-relevant.
  ordered <- res[order(-res$non_cured_hazard_ratio), ]
  expect_true(all(diff(ordered$total_cost) <= 1e-6))
  for (col in c("nmb_at_50000", "nmb_at_100000", "nmb_at_150000")) {
    expect_true(all(diff(ordered[[col]]) >= -1e-6))
  }
})
