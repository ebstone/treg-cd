# testthat::test_dir() runs with the working directory set to tests/testthat/, but this
# project's scripts source() each other with repo-root-relative paths (README.md convention).
# Same loading pattern as test-costs-utilities.R.
repo_root_relative <- function(...) file.path("..", "..", ...)
source(repo_root_relative("R", "utils", "transition_matrix.R"), local = TRUE)
source(repo_root_relative("R", "utils", "life_table.R"), local = TRUE)
source(repo_root_relative("R", "utils", "population_utility.R"), local = TRUE)
source(repo_root_relative("R", "utils", "refractory_multipliers.R"), local = TRUE)
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

  # The relapse hazard must be re-stated here to match whatever headroom_pi_star() itself used --
  # i.e. its own default, the 10-year anchor since 2026-08-06 (B3). This argument is passed by
  # NAME, not positionally as it once was: the old positional `0` in this slot silently became a
  # different assumption from the function under test the moment that default changed, which is
  # the same class of mistake B3 itself was.
  treg_at_root <- run_treg_arm_lifetime(
    HORIZON_CYCLES_6YR, res$pi_star,
    relapse_hazard_annual = duration_to_hazard(DEFAULT_RELAPSE_DURATION_YEARS),
    price_usd = 5000, matrices = matrices, raw_dir = RAW_DIR, proc_dir = PROC_DIR
  )
  expect_equal(net_monetary_benefit(treg_at_root, 150000), target_nmb, tolerance = 1e-3)
})

test_that("headroom_frontier returns a monotone-non-decreasing pi_star as price rises, over a grid with a real feasible/infeasible boundary", {
  # price_grid recalibrated 2026-08-06 (B1 fix, docs/model_audit_v6.md A17): the corrected,
  # larger comparator responder pools push comparator maintenance drug cost up enough that the
  # feasible/infeasible boundary at the 6.15-year horizon (this test's default n_cycles) moved
  # from inside [1000, 20000] to inside (20000, 30000] -- more headroom for Treg at a given price,
  # not less, since a bigger comparator responder pool raises comparator drug spend more than it
  # raises comparator QALYs. 20000 is still feasible, 30000 is not (verified directly against
  # headroom_frontier() post-fix); this grid just needs to keep bracketing a real boundary, the
  # same intent as the original grid, not a specific price level.
  price_grid <- c(1000, 5000, 10000, 20000, 30000)
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

# ---- h sweep (relapse duration T) and the pi * g(h) factorisation check (2026-08-05) -----------

test_that("headroom_pi_star's qaly_gain is 0 exactly at pi_star = 0, NA when infeasible, and positive at an interior root", {
  cheap <- headroom_pi_star(price_usd = 1, wtp_usd = 150000, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_equal(cheap$pi_star, 0)
  expect_equal(cheap$qaly_gain, 0)

  absurd <- headroom_pi_star(price_usd = 10000000, wtp_usd = 150000, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_false(absurd$feasible)
  expect_true(is.na(absurd$qaly_gain))

  interior <- headroom_pi_star(price_usd = 5000, wtp_usd = 150000, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_true(interior$feasible)
  expect_true(interior$pi_star > 0 && interior$pi_star < 1)
  expect_true(interior$qaly_gain > 0)
})

test_that("a shorter median SDR duration T (higher relapse hazard) requires a higher pi* for the same price -- the h-sweep sanity check", {
  # price = 5000, not 3000/7000 (recalibrated 2026-08-06, B1 fix, docs/model_audit_v6.md A17): the
  # corrected, larger comparator responder pools raised comparator maintenance drug cost enough
  # that 3000 is now cheap enough to be feasible at pi_star = 0 under EITHER duration (a boundary
  # case, not the interior comparison this test needs), while at the 6.15-year horizon a
  # pessimistic T = 2 relapse duration pushes 7000 out of the feasible region entirely (the same
  # failure mode the original comment already knew to avoid, just at a different price now).
  # 5000 keeps both durations in the interior-root regime post-fix (verified directly against
  # headroom_pi_star()).
  price <- 5000
  optimistic <- headroom_pi_star(price, wtp_usd = 150000, relapse_hazard_annual = duration_to_hazard(20),
                                  raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  pessimistic <- headroom_pi_star(price, wtp_usd = 150000, relapse_hazard_annual = duration_to_hazard(2),
                                   raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_true(optimistic$feasible && pessimistic$feasible)
  expect_true(pessimistic$pi_star > optimistic$pi_star)
})

test_that("headroom_pi_star's DEFAULT relapse hazard is the 10-year anchor, not 0 -- the regression test for B3's deterministic half", {
  # B3 (peer review 2026-08-05, docs/model_audit_v6.md A19) was two defects wearing one name. The
  # PSA half (h never sampled) is tested in test-psa.R. THIS is the deterministic half: the h sweep
  # tooling landed 2026-08-05 with DEFAULT_RELAPSE_DURATION_YEARS <- 10 recorded as the base-case
  # anchor, but headroom_pi_star()/headroom_frontier()/ejp_deterministic()/ejp_frontier() kept
  # defaulting to h = 0, so every headline number was computed assuming a permanent cure. Nothing
  # in the suite caught it, because every h-aware test passed h explicitly and every other test
  # only ever compared the default against itself. This test compares the default against a
  # SOURCED value and against the value it used to be.
  price <- 5000
  default_result <- headroom_pi_star(price, wtp_usd = 150000, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  at_anchor <- headroom_pi_star(price, wtp_usd = 150000,
                                 relapse_hazard_annual = duration_to_hazard(DEFAULT_RELAPSE_DURATION_YEARS),
                                 raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  at_zero <- headroom_pi_star(price, wtp_usd = 150000, relapse_hazard_annual = 0,
                               raw_dir = RAW_DIR, proc_dir = PROC_DIR)

  # The default IS the anchor...
  expect_equal(default_result$pi_star, at_anchor$pi_star)
  expect_equal(default_result$qaly_gain, at_anchor$qaly_gain)
  # ...and is NOT the old permanent-cure assumption. Asserted as a strict inequality in the
  # correct direction, not merely as inequality: a decaying cure is worth less, so it takes a
  # LARGER cure fraction to break even at the same price.
  expect_true(at_zero$feasible && default_result$feasible)
  expect_true(default_result$pi_star > at_zero$pi_star)
  expect_false(isTRUE(all.equal(default_result$pi_star, at_zero$pi_star)))

  # headroom_frontier() must default identically -- it forwards to headroom_pi_star() per price
  # point, and a divergence between the two defaults would silently make the frontier and the
  # single-point headline number two different curves.
  price_grid <- c(2000, 5000)
  frontier_default <- headroom_frontier(price_grid, wtp_usd = 150000, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_equal(unique(frontier_default$relapse_hazard_annual),
                duration_to_hazard(DEFAULT_RELAPSE_DURATION_YEARS))
  expect_equal(frontier_default$pi_star[frontier_default$price_usd == price], at_anchor$pi_star)
})

test_that("the relapse hazard is structurally inert at pi = 0, so run_base_case()'s explicit h = 0 is not a stale default", {
  # run_base_case() (and run_refractory_scenario(), and S12) still pass relapse_hazard_annual = 0
  # explicitly at pi_sdr = 0, and their comments claim that is provably harmless rather than an
  # oversight left behind by B3's fix. This asserts the claim instead of trusting the comment: with
  # pi = 0 the SDR track is seeded with remission_mass * 0 and every subsequent update multiplies
  # that mass, so no h can ever change the trace. Checked to exact floating-point equality, since
  # the argument is that the arithmetic is literally identical, not merely close.
  matrices <- build_all_transition_matrices(RAW_DIR)
  hazards <- c(0, duration_to_hazard(DEFAULT_RELAPSE_DURATION_YEARS), duration_to_hazard(2), 5)
  at_pi_zero <- lapply(hazards, function(h) {
    run_treg_arm_lifetime(HORIZON_CYCLES_6YR, pi_sdr = 0, relapse_hazard_annual = h,
                           price_usd = 15000, matrices = matrices, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  })
  for (i in seq_along(hazards)[-1]) {
    expect_identical(at_pi_zero[[i]]$qalys, at_pi_zero[[1]]$qalys, info = paste("h =", hazards[i]))
    expect_identical(at_pi_zero[[i]]$total_cost, at_pi_zero[[1]]$total_cost, info = paste("h =", hazards[i]))
  }

  # The converse, so this isn't passing because the hazard is broken everywhere: at pi > 0 the
  # very same comparison must show a real difference.
  at_pi_half <- lapply(c(0, duration_to_hazard(2)), function(h) {
    run_treg_arm_lifetime(HORIZON_CYCLES_6YR, pi_sdr = 0.5, relapse_hazard_annual = h,
                           price_usd = 15000, matrices = matrices, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  })
  expect_true(at_pi_half[[1]]$qalys - at_pi_half[[2]]$qalys > 1e-4)
})

test_that("headroom_frontier_by_duration stacks one full headroom_frontier() per duration, each tagged correctly", {
  price_grid <- c(1000, 5000, 10000)
  duration_grid <- c(5, 10, Inf)
  surf <- headroom_frontier_by_duration(price_grid, wtp_usd = 150000, duration_grid_years = duration_grid,
                                         raw_dir = RAW_DIR, proc_dir = PROC_DIR)

  expect_equal(nrow(surf), length(price_grid) * length(duration_grid))
  expect_setequal(surf$duration_years, duration_grid)
  # T = Inf (permanent remission) should reproduce an EXPLICIT h = 0 headroom_frontier() call
  # exactly -- same underlying computation, just reached through duration_to_hazard(Inf) = 0.
  # `relapse_hazard_annual = 0` has to be passed explicitly here as of 2026-08-06 (B3): the
  # function's default is no longer 0 but the 10-year anchor, so leaving it implicit would compare
  # the T = Inf row against the T = 10 base case and (correctly) fail.
  direct <- headroom_frontier(price_grid, wtp_usd = 150000, relapse_hazard_annual = 0,
                               raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  from_surface <- surf[surf$duration_years == Inf, ]
  expect_equal(from_surface$pi_star[order(from_surface$price_usd)], direct$pi_star[order(direct$price_usd)])

  # A shorter duration (higher hazard) can never make headroom EASIER to reach at the same price --
  # pi* should be non-decreasing as duration shortens, mirroring the single-duration sanity check.
  feasible <- surf[surf$feasible, ]
  by_price <- split(feasible, feasible$price_usd)
  for (rows in by_price) {
    ord <- order(rows$duration_years)  # ascending duration = descending hazard
    expect_true(all(diff(rows$pi_star[ord]) <= 1e-8))
  }
})

test_that("verify_pi_factorization: pi * g(h) is exact to floating-point precision, not merely close -- the finding this function's docstring reports", {
  # The resolutions memo hedged this as "very good but not exact." It's checked here to actually
  # be exact (run_treg_arm()'s pi-split feeds into a purely linear recursion -- R/05's own
  # docstring on verify_pi_factorization() spells out why) -- residuals should be ordinary
  # floating-point noise (~1e-10 relative or smaller), not a genuine few-percent approximation
  # gap. Checked at both a mid-range price and both the 6.15-year and lifetime horizons, since the
  # linearity argument doesn't depend on either.
  for (n_cycles_arg in list(list(n_cycles = HORIZON_CYCLES_6YR, baseline_age = NULL),
                            list(n_cycles = HORIZON_CYCLES_LIFETIME, baseline_age = ASSUMED_PATIENT_AGE_YEARS))) {
    res <- do.call(verify_pi_factorization, c(
      list(price_usd = 15000, relapse_hazard_annual = duration_to_hazard(10),
           pi_grid = seq(0, 1, by = 0.25), raw_dir = RAW_DIR, proc_dir = PROC_DIR),
      n_cycles_arg
    ))
    expect_equal(nrow(res), 5)
    # Exact by construction at the endpoints: g(h) is defined as the pi=1-vs-pi=0 gap itself.
    expect_equal(res$observed_qaly_gain[res$pi == 0], 0)
    expect_equal(res$predicted_qaly_gain[res$pi == 0], 0)
    expect_equal(res$observed_qaly_gain[res$pi == 1], res$predicted_qaly_gain[res$pi == 1], tolerance = 1e-9)
    # Interior points: floating-point-scale residuals, not a genuine approximation gap.
    expect_true(all(res$rel_error < 1e-8))
    expect_true(all(is.finite(res$abs_error)))
  }
})
