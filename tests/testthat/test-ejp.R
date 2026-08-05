# testthat::test_dir() runs with the working directory set to tests/testthat/, but this
# project's scripts source() each other with repo-root-relative paths (README.md convention).
# Same loading pattern as test-psa.R / test-evpi-evppi.R.
repo_root_relative <- function(...) file.path("..", "..", ...)
source(repo_root_relative("R", "utils", "transition_matrix.R"), local = TRUE)
source(repo_root_relative("R", "00_derive_transition_probs.R"), local = TRUE)
source(repo_root_relative("R", "01_decision_tree.R"), local = TRUE)
source(repo_root_relative("R", "02_markov_engine.R"), local = TRUE)
source(repo_root_relative("R", "03_cure_fraction_module.R"), local = TRUE)
source(repo_root_relative("R", "04_costs_utilities.R"), local = TRUE)
source(repo_root_relative("R", "05_deterministic_results.R"), local = TRUE)
source(repo_root_relative("R", "06_psa.R"), local = TRUE)
source(repo_root_relative("R", "08_ejp.R"), local = TRUE)

RAW_DIR <- repo_root_relative("data", "raw")
PROC_DIR <- repo_root_relative("data", "processed")

test_that("load_treg_cogs reads the sourced COGS figure exactly", {
  expect_equal(load_treg_cogs(PROC_DIR), 4979.187449068181)
})

test_that("gross_margin_over_cogs is a simple, hand-computed formula", {
  cogs <- load_treg_cogs(PROC_DIR)
  expect_equal(gross_margin_over_cogs(cogs * 2, PROC_DIR), 0.5)
  expect_equal(gross_margin_over_cogs(cogs, PROC_DIR), 0)
  expect_true(gross_margin_over_cogs(cogs / 2, PROC_DIR) < 0)  # priced below COGS -> negative margin
})

test_that("ejp_deterministic's p_star round-trips through R/05's independently-implemented headroom_pi_star", {
  # This is the strongest available check: two independently-derived methods (closed-form
  # algebra here vs. uniroot() in R/05_deterministic_results.R) solving the same INMB=0 equation
  # from opposite directions must agree, or one of them has a real bug -- exactly how a genuine
  # bug (using non_drug_cost instead of total_cost) was caught while building this module.
  for (pi in c(0, 0.3, 0.6, 1)) {
    res <- ejp_deterministic(pi, wtp_usd = 150000, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
    rt <- headroom_pi_star(res$p_star, wtp_usd = 150000, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
    expect_true(rt$feasible)
    expect_equal(rt$pi_star, pi, tolerance = 1e-3)
  }
})

test_that("ejp_deterministic's p_star makes INMB exactly 0 against the identified comparator, verified independently", {
  res <- ejp_deterministic(0.4, wtp_usd = 100000, raw_dir = RAW_DIR, proc_dir = PROC_DIR)

  matrices <- build_all_transition_matrices(RAW_DIR)
  treg_at_pstar <- run_treg_arm_lifetime(HORIZON_CYCLES_6YR, 0.4, 0, res$p_star, matrices,
                                          raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  comp_summary <- run_comparator_arm_lifetime(res$comparator, HORIZON_CYCLES_6YR, matrices,
                                               raw_dir = RAW_DIR, proc_dir = PROC_DIR)

  treg_nmb <- net_monetary_benefit(treg_at_pstar, 100000)
  comp_nmb <- net_monetary_benefit(comp_summary, 100000)
  expect_equal(treg_nmb, comp_nmb, tolerance = 1e-6)
})

test_that("ejp_frontier and headroom_frontier trace the same (pi, price) curve from opposite directions", {
  pi_grid <- c(0, 0.2, 0.5, 0.8, 1)
  ef <- ejp_frontier(pi_grid, wtp_usd = 150000, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_equal(nrow(ef), length(pi_grid))
  expect_true(all(diff(ef$p_star) >= 0))  # higher pi can only justify an equal or higher price

  for (i in seq_along(pi_grid)) {
    rt <- headroom_pi_star(ef$p_star[i], wtp_usd = 150000, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
    expect_equal(rt$pi_star, pi_grid[i], tolerance = 1e-3)
  }
})

test_that("ejp_probabilistic_draws' p_star is an exact algebraic root: plugging it back in recovers the comparator's NMB exactly", {
  psa <- run_psa(n_draws = 100, raw_dir = RAW_DIR, proc_dir = PROC_DIR, seed = 21)
  draws <- ejp_probabilistic_draws(psa, wtp_usd = 100000)
  expect_equal(nrow(draws), 100)

  treg <- psa[psa$intervention == "TREG", ]
  treg <- treg[match(draws$draw, treg$draw), ]
  cost_at_pstar <- treg$total_cost - treg$treg_price + draws$p_star
  nmb_at_pstar <- 100000 * treg$qalys - cost_at_pstar

  psa_key <- paste(psa$draw, psa$intervention)
  comp_rows <- psa[match(paste(draws$draw, draws$comparator), psa_key), ]
  target_nmb <- 100000 * comp_rows$qalys - comp_rows$total_cost

  expect_equal(nmb_at_pstar, target_nmb, tolerance = 1e-6)
})

test_that("ejp_probabilistic's median lies within its own 95% credible interval, and the draws frame matches its summary", {
  psa <- run_psa(n_draws = 150, raw_dir = RAW_DIR, proc_dir = PROC_DIR, seed = 22)
  summ <- ejp_probabilistic(psa, wtp_usd = 100000)
  expect_true(summ$ci_low <= summ$median && summ$median <= summ$ci_high)
  expect_equal(nrow(summ$draws), 150)
  expect_equal(unname(stats::median(summ$draws$p_star)), summ$median)
})

test_that("ejp_p50 matches a hand-computed step-function crossing on a tiny toy sample", {
  # 4 draws' worth of p_star: 10, 20, 30, 40. Survival fractions at each: 4/4, 3/4, 2/4, 1/4.
  # The first value where survival <= 0.5 is 30 (survival exactly 0.5 there).
  expect_equal(ejp_p50(c(10, 20, 30, 40)), 30)
})

test_that("ejp_p50 closely agrees with median() on a large, real PSA-derived sample (expected under near-symmetry, per analysis_plan.md sec 9.1)", {
  psa <- run_psa(n_draws = 400, raw_dir = RAW_DIR, proc_dir = PROC_DIR, seed = 23)
  draws <- ejp_probabilistic_draws(psa, wtp_usd = 100000)
  p50 <- ejp_p50(draws$p_star)
  med <- stats::median(draws$p_star)
  # Independently computed (module header: NOT derived from each other) -- expected close, not
  # required identical; a wide but real tolerance, not a tautological check.
  expect_equal(p50, med, tolerance = 0.02)
})

test_that("D_BAR is 1, matching the single-dose base case enforced elsewhere (treg_dose_cost(), run_treg_arm_lifetime())", {
  expect_equal(D_BAR, 1)
})
