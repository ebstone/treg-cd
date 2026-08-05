# testthat::test_dir() runs with the working directory set to tests/testthat/, but this
# project's scripts source() each other with repo-root-relative paths (README.md convention).
# Same loading pattern as test-costs-utilities.R.
repo_root_relative <- function(...) file.path("..", "..", ...)
source(repo_root_relative("R", "utils", "transition_matrix.R"), local = TRUE)
source(repo_root_relative("R", "utils", "life_table.R"), local = TRUE)
source(repo_root_relative("R", "00_derive_transition_probs.R"), local = TRUE)
source(repo_root_relative("R", "01_decision_tree.R"), local = TRUE)
source(repo_root_relative("R", "02_markov_engine.R"), local = TRUE)
source(repo_root_relative("R", "03_cure_fraction_module.R"), local = TRUE)
source(repo_root_relative("R", "04_costs_utilities.R"), local = TRUE)
source(repo_root_relative("R", "05_deterministic_results.R"), local = TRUE)

RAW_DIR <- repo_root_relative("data", "raw")
PROC_DIR <- repo_root_relative("data", "processed")

test_that("net_monetary_benefit is a simple, hand-computed WTP*QALY - cost", {
  s <- list(qalys = 4, total_cost = 100000)
  expect_equal(net_monetary_benefit(s, 50000), 4 * 50000 - 100000)
  expect_equal(net_monetary_benefit(s, 0), -100000)
})

test_that("best_comparator_nmb picks the arm with the highest NMB, not the lowest cost or highest QALYs alone", {
  summaries <- list(
    A = list(qalys = 4, total_cost = 100000),   # NMB@100k = 300000
    B = list(qalys = 3, total_cost = 50000),    # NMB@100k = 250000
    C = list(qalys = 5, total_cost = 460000)    # NMB@100k = 40000 -- highest QALYs, worst NMB
  )
  res <- best_comparator_nmb(summaries, 100000)
  expect_equal(res$comparator, "A")
  expect_equal(res$nmb, 300000)
})

test_that("treg_price_dependent_dose_cost's composition agrees with treg_dose_cost() when the same price is used", {
  prices <- load_drug_prices(RAW_DIR)
  sourced_acquisition <- load_treg_dose_acquisition_cost(PROC_DIR)

  from_treg_dose_cost <- treg_dose_cost(
    cyclophosphamide_dose_mg = 25, observation_stay_cost_usd = 2672.15,
    proc_dir = PROC_DIR, raw_dir = RAW_DIR, prices = prices
  )
  from_price_dependent <- treg_price_dependent_dose_cost(
    price_usd = sourced_acquisition, cyclophosphamide_dose_mg = 25,
    observation_stay_cost_usd = 2672.15, prices = prices
  )
  expect_equal(from_price_dependent, from_treg_dose_cost)
})

test_that("run_comparator_arm_lifetime runs cleanly for all three comparators with sane, positive outputs", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  for (tx in COMPARATOR_THERAPIES) {
    res <- run_comparator_arm_lifetime(tx, HORIZON_CYCLES_6YR, matrices, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
    expect_true(res$qalys > 0 && res$qalys < HORIZON_CYCLES_6YR * 2 / 52)  # can't exceed the horizon in years
    expect_true(res$total_cost > 0)
    expect_false(anyNA(unlist(res)))
  }
})

test_that("run_base_case's TREG row at pi=0 is exactly UST-equivalent on QALYs (Null-scenario check, analysis_plan.md sec 7.2)", {
  bc <- run_base_case(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_setequal(bc$intervention, c("UST", "IFX", "ADA", "TREG"))

  ust_qalys <- bc$qalys[bc$intervention == "UST"]
  treg_qalys <- bc$qalys[bc$intervention == "TREG"]
  expect_equal(treg_qalys, ust_qalys)

  # Treg's one-time dose cost differs from UST's cumulative drug cost over the horizon, so total
  # cost should NOT match even though QALYs do -- this is the whole reason it's worth reporting
  # as a separate row rather than just noting "same as UST".
  ust_cost <- bc$total_cost[bc$intervention == "UST"]
  treg_cost <- bc$total_cost[bc$intervention == "TREG"]
  expect_false(isTRUE(all.equal(ust_cost, treg_cost)))
})

test_that("run_base_case's nmb_at_* columns are named with plain integers, not scientific notation", {
  bc <- run_base_case(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  # Regression test for a real bug caught during development: paste0("nmb_at_", 100000) renders
  # as "nmb_at_1e+05", which rbind.data.frame()'s make.names() then mangles further to
  # "nmb_at_1e.05" -- silently losing the "nmb_at_100000" column a caller would actually look for.
  expect_true(all(paste0("nmb_at_", format(WTP_THRESHOLDS_USD, scientific = FALSE, trim = TRUE)) %in% names(bc)))
  expect_false(any(grepl("e[+-]", names(bc))))
})

test_that("headroom_pi_star: a very cheap price is already cost-effective at pi=0, an absurd price is infeasible even at pi=1", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  comparator_summaries <- stats::setNames(
    lapply(COMPARATOR_THERAPIES, function(tx) {
      run_comparator_arm_lifetime(tx, HORIZON_CYCLES_6YR, matrices, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
    }),
    COMPARATOR_THERAPIES
  )

  cheap <- headroom_pi_star(price_usd = 1, wtp_usd = 150000, matrices = matrices,
                             comparator_summaries = comparator_summaries, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_equal(cheap$pi_star, 0)
  expect_true(cheap$feasible)

  absurd <- headroom_pi_star(price_usd = 10000000, wtp_usd = 150000, matrices = matrices,
                              comparator_summaries = comparator_summaries, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_true(is.na(absurd$pi_star))
  expect_false(absurd$feasible)
})

test_that("headroom_pi_star's root, where one exists, actually satisfies the NMB break-even condition", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  comparator_summaries <- stats::setNames(
    lapply(COMPARATOR_THERAPIES, function(tx) {
      run_comparator_arm_lifetime(tx, HORIZON_CYCLES_6YR, matrices, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
    }),
    COMPARATOR_THERAPIES
  )
  target_nmb <- best_comparator_nmb(comparator_summaries, 150000)$nmb

  res <- headroom_pi_star(price_usd = 5000, wtp_usd = 150000, matrices = matrices,
                           comparator_summaries = comparator_summaries, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_true(res$feasible)
  expect_true(res$pi_star > 0 && res$pi_star < 1)  # a genuine interior root, not a boundary case

  treg_at_root <- run_treg_arm_lifetime(HORIZON_CYCLES_6YR, res$pi_star, 0, 5000, matrices,
                                         raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_equal(net_monetary_benefit(treg_at_root, 150000), target_nmb, tolerance = 1e-3)
})

test_that("headroom_frontier returns a monotone-non-decreasing pi_star as price rises, over a grid with a real feasible/infeasible boundary", {
  price_grid <- c(1000, 5000, 10000, 15000, 20000)
  f <- headroom_frontier(price_grid, wtp_usd = 150000, raw_dir = RAW_DIR, proc_dir = PROC_DIR)

  expect_equal(nrow(f), length(price_grid))
  expect_true(any(f$feasible))       # the grid brackets a real feasible region
  expect_true(any(!f$feasible))      # ...and a real infeasible one, not just all-or-nothing

  feasible_rows <- f[f$feasible, ]
  ord <- order(feasible_rows$price_usd)
  # pi* must not decrease as price rises: a more expensive dose can only ever require an equal or
  # higher cure fraction to break even, never a lower one.
  expect_true(all(diff(feasible_rows$pi_star[ord]) >= -1e-8))
})

test_that("a longer horizon gives Treg strictly more headroom (lower pi* for the same price) -- sanity-checks that discounting/compounding flows through correctly", {
  # price = 7000, not 10000: after the A16 fix (deterministic utility base case rebuilt on 0.82,
  # not the retired 0.9554 snapshot draw -- docs/model_audit_v6.md), every arm's QALYs are lower,
  # so a price that used to leave the 6-year horizon feasible at pi=1 no longer does (10000 is now
  # infeasible even at the 6yr horizon's own pi=1 -- correctly, not a bug: less durable benefit per
  # cure buys less headroom). 7000 keeps both horizons in the interior-root regime this test is
  # actually about.
  price <- 7000
  short <- headroom_pi_star(price, wtp_usd = 150000, n_cycles = HORIZON_CYCLES_6YR, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  long <- headroom_pi_star(price, wtp_usd = 150000, n_cycles = HORIZON_CYCLES_10YR, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_true(short$feasible && long$feasible)
  expect_true(long$pi_star < short$pi_star)
})

# ---- Lifetime horizon (2026-08-05) -------------------------------------------------------------

test_that("baseline_age = NULL (default) reproduces the exact pre-lifetime-horizon base case, unchanged", {
  without_arg <- run_base_case(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  with_explicit_null <- run_base_case(raw_dir = RAW_DIR, proc_dir = PROC_DIR, baseline_age = NULL)
  expect_equal(without_arg, with_explicit_null)
})

test_that("HORIZON_CYCLES_LIFETIME is derived to land the cohort exactly at the life table's terminal age", {
  # ASSUMED_PATIENT_AGE_YEARS (35) + HORIZON_CYCLES_LIFETIME cycles of 2 weeks each should reach
  # age 100 (data/raw/nchs_us_life_tables_2021.csv's terminal row) on the LAST transition.
  last_cycle_age <- floor(ASSUMED_PATIENT_AGE_YEARS + (HORIZON_CYCLES_LIFETIME - 1) * 2 / 52)
  expect_equal(last_cycle_age, 100)
  # One cycle earlier must NOT yet have reached it -- confirms this is the minimal such horizon,
  # not an arbitrary large number that happens to be big enough.
  second_last_cycle_age <- floor(ASSUMED_PATIENT_AGE_YEARS + (HORIZON_CYCLES_LIFETIME - 2) * 2 / 52)
  expect_true(second_last_cycle_age < 100)
})

test_that("run_base_case at the lifetime horizon produces QALYs far above the 6.15-year horizon, and the cohort is fully exhausted by the end", {
  bc_short <- run_base_case(n_cycles = HORIZON_CYCLES_6YR, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  bc_lifetime <- run_base_case(n_cycles = HORIZON_CYCLES_LIFETIME, baseline_age = ASSUMED_PATIENT_AGE_YEARS,
                                raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  # A lifetime run should report meaningfully more QALYs per arm than the 6.15-year comparability
  # scenario -- the entire point of extending the horizon (module header).
  for (arm in bc_lifetime$intervention) {
    expect_true(
      bc_lifetime$qalys[bc_lifetime$intervention == arm] > bc_short$qalys[bc_short$intervention == arm]
    )
  }
})

test_that("at the lifetime horizon, deterministic EJP at pi=1 clears Treg's sourced acquisition price -- the finding this work exists to check", {
  # At the 6.15-year horizon (documented in the run_full_analysis.R headline finding, 2026-08-05)
  # the deterministic EJP ceiling even at pi=1 was ~$8,724 at WTP=$150k -- well under the sourced
  # $19,917 acquisition cost, an artefact of the horizon being too short to capture a cure's
  # value. At the lifetime horizon this should no longer be true.
  treg_price <- load_treg_dose_acquisition_cost(PROC_DIR)
  res <- headroom_pi_star(treg_price, wtp_usd = 150000, n_cycles = HORIZON_CYCLES_LIFETIME,
                           baseline_age = ASSUMED_PATIENT_AGE_YEARS, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_true(res$feasible)
  expect_true(res$pi_star > 0 && res$pi_star < 1)  # a genuinely interior, non-trivial cure fraction
})

test_that("headroom_frontier's baseline_age argument is actually threaded through, not silently ignored", {
  # Same price, same WTP: the 6.15-year comparability scenario and the lifetime scenario should
  # disagree (a longer horizon with age-appropriate mortality changes the answer materially) --
  # if baseline_age/n_cycles were being ignored somewhere in the plumbing, these would coincide.
  price <- 15000
  short <- headroom_pi_star(price, wtp_usd = 150000, n_cycles = HORIZON_CYCLES_6YR, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  lifetime <- headroom_pi_star(price, wtp_usd = 150000, n_cycles = HORIZON_CYCLES_LIFETIME,
                                baseline_age = ASSUMED_PATIENT_AGE_YEARS, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_false(isTRUE(all.equal(short$pi_star, lifetime$pi_star)) && short$feasible == lifetime$feasible)
})
