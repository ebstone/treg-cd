# testthat::test_dir() runs with the working directory set to tests/testthat/, but this
# project's scripts source() each other with repo-root-relative paths (README.md convention).
# Same loading pattern as test-deterministic-results.R.
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
source(repo_root_relative("R", "06_psa.R"), local = TRUE)

RAW_DIR <- repo_root_relative("data", "raw")
PROC_DIR <- repo_root_relative("data", "processed")

test_that("gamma_shape_scale reproduces analysis_plan.md's own cited Treg-price parameterisation", {
  g <- gamma_shape_scale(19917, 5080.87)
  # analysis_plan.md sec 7.1 item 12: "Gamma (alpha ~= 15.4, scale ~= 1,297)" -- their own rounding
  # of 1296.14 to the nearest 1, so check within 1, not exact round-trip equality.
  expect_equal(round(g$shape, 1), 15.4)
  expect_equal(g$scale, 1297, tolerance = 1)
})

test_that("gamma_shape_scale is a simple, hand-computed method-of-moments formula", {
  g <- gamma_shape_scale(mean = 100, sd = 10)
  expect_equal(g$shape, 100)
  expect_equal(g$scale, 1)
  expect_equal(g$shape * g$scale, 100)                # recovers the mean
  expect_equal(sqrt(g$shape) * g$scale, 10)            # recovers the sd
})

test_that("load_psa_distributions reads the real sourced file with the expected rows", {
  psa_df <- load_psa_distributions(PROC_DIR)
  expect_setequal(
    psa_df$variable,
    c("Remission", "Mild:Remission Ratio", "M-SR:Mild Ratio", "M-S:M-SR Ratio",
      "Surgery:M-SR Ratio", "psa_costs_ifx_dose", "psa_cost_treg_dose")
  )
})

test_that("psa_uniform_bounds and psa_gamma_params extract the right row and reject the wrong distribution type", {
  psa_df <- load_psa_distributions(PROC_DIR)
  b <- psa_uniform_bounds(psa_df, "Remission")
  expect_equal(unname(b), c(0.66, 0.98))
  expect_error(psa_gamma_params(psa_df, "Remission"))  # Remission is Uniform, not Gamma

  g <- psa_gamma_params(psa_df, "psa_cost_treg_dose")
  expect_equal(round(g$shape, 1), 15.4)
  expect_error(psa_uniform_bounds(psa_df, "psa_cost_treg_dose"))  # Gamma, not Uniform
})

test_that("sample_utility_chain_draws produces a strictly-decreasing chain within sourced bounds, Death always 0", {
  set.seed(42)
  draws <- sample_utility_chain_draws(200, PROC_DIR)
  expect_equal(dim(draws), c(200, length(MAINTENANCE_STATES)))
  expect_true(all(draws[, "Death"] == 0))

  # Every ratio in the sourced file is < 1 (0.79-0.98), so the chain must be strictly decreasing:
  # Remission > Mild > Moderate-Severe Responder > Moderate-Severe, and Remission > M-SR > Surgery.
  expect_true(all(draws[, "Remission"] > draws[, "Mild"]))
  expect_true(all(draws[, "Mild"] > draws[, "Moderate-Severe Responder"]))
  expect_true(all(draws[, "Moderate-Severe Responder"] > draws[, "Moderate-Severe"]))
  expect_true(all(draws[, "Moderate-Severe Responder"] > draws[, "Surgery"]))

  # Remission itself must stay within its sourced Uniform(0.66, 0.98) bounds.
  expect_true(all(draws[, "Remission"] >= 0.66 & draws[, "Remission"] <= 0.98))
})

test_that("sample_utility_chain_draws is reproducible under a fixed seed", {
  set.seed(7)
  a <- sample_utility_chain_draws(20, PROC_DIR)
  set.seed(7)
  b <- sample_utility_chain_draws(20, PROC_DIR)
  expect_equal(a, b)
})

test_that("sample_treg_price_draws is always positive and has approximately the sourced mean over many draws", {
  set.seed(123)
  draws <- sample_treg_price_draws(20000, PROC_DIR)
  expect_true(all(draws > 0))
  expect_equal(mean(draws), 19917, tolerance = 0.02)          # within 2% over 20k draws
  expect_equal(sd(draws), 5080.87, tolerance = 0.05)          # within 5% over 20k draws
})

# ---- Post-cure relapse hazard sampling (B3, 2026-08-06) ----------------------------------------

test_that("sample_relapse_hazard_draws recovers the 10-year anchor as its MEAN, with the intended Gamma spread", {
  set.seed(1234)
  draws <- sample_relapse_hazard_draws(50000)
  target_mean <- duration_to_hazard(DEFAULT_RELAPSE_DURATION_YEARS)

  expect_true(all(draws > 0))                                  # a hazard is strictly positive
  # The mean, specifically, is the anchored quantity -- that is what makes the PSA the stochastic
  # counterpart of the deterministic base case rather than a differently-centred analysis.
  expect_equal(mean(draws), target_mean, tolerance = 0.02)     # within 2% at 50k draws
  # Gamma(shape k, mean m) has sd = m/sqrt(k), i.e. CV = 1/sqrt(k) -- re-derived here from the
  # closed form rather than read off the implementation, so a wrong scale/rate parameterisation
  # (the classic rgamma() footgun) would fail this even though the mean still came out right.
  expect_equal(sd(draws), target_mean / sqrt(RELAPSE_HAZARD_PSA_GAMMA_SHAPE), tolerance = 0.03)
})

test_that("the sampled hazard's implied SDR durations sit where sample_relapse_hazard_draws' docstring says they do, against RELAPSE_DURATION_GRID_YEARS", {
  # The percentile table in that docstring is the actual justification offered to a reviewer for
  # the chosen shape, so it is asserted rather than left as prose that could drift from the code.
  # Everything here is computed from the distribution, then translated into T-years through the
  # module's own hazard_to_duration() -- the units the sweep grid is stated in.
  set.seed(4321)
  draws <- sample_relapse_hazard_draws(100000)
  implied_T <- hazard_to_duration_vec <- vapply(draws, hazard_to_duration, numeric(1))

  # Median implied duration ~11.9 years: near the 10-year anchor, and (deliberately) a little
  # above it, since fixing the MEAN of a right-skewed hazard puts the median hazard below the mean.
  expect_equal(stats::median(implied_T), 11.9, tolerance = 0.05)
  # Interquartile range brackets the anchor and its two neighbours on RELAPSE_DURATION_GRID_YEARS.
  iqr <- stats::quantile(implied_T, c(0.25, 0.75))
  expect_true(iqr[["75%"]] > DEFAULT_RELAPSE_DURATION_YEARS)
  expect_true(iqr[["25%"]] < DEFAULT_RELAPSE_DURATION_YEARS)
  expect_true(iqr[["25%"]] > 5 && iqr[["75%"]] < 25)

  # The sweep grid's pessimistic endpoint (T = 2 years) must be IN the sampled tail but well out
  # of the body: present as a possibility, not as a central expectation. ~5 draws in 10,000.
  frac_worse_than_2yr <- mean(implied_T < 2)
  expect_true(frac_worse_than_2yr > 0 && frac_worse_than_2yr < 0.002)
  # Permanence (T = Inf, the Ovasave/CATS1-style upper bound) must have essentially no mass: with
  # shape > 1 the Gamma density vanishes at h = 0, which is the whole reason shape 1 (Exponential)
  # was rejected. Nothing beyond a human lifetime from a baseline age of 35.
  expect_true(mean(implied_T > 200) < 0.005)
})

test_that("sample_relapse_hazard_draws is reproducible under a fixed seed and rejects nonsense parameters", {
  set.seed(8)
  a <- sample_relapse_hazard_draws(50)
  set.seed(8)
  b <- sample_relapse_hazard_draws(50)
  expect_equal(a, b)
  expect_error(sample_relapse_hazard_draws(10, mean_hazard = 0))
  expect_error(sample_relapse_hazard_draws(10, mean_hazard = -0.1))
  expect_error(sample_relapse_hazard_draws(10, shape = 0))
})

test_that("sample_dirichlet_row's rows sum to 1 and its mean converges to the target probs", {
  set.seed(99)
  probs <- c(0.7, 0.2, 0.1)
  draws <- sample_dirichlet_row(20000, probs, concentration = 200)
  expect_equal(rowSums(draws), rep(1, 20000))
  expect_equal(colMeans(draws), probs, tolerance = 0.02)

  # Higher concentration must give tighter (lower-variance) draws around the same mean.
  tight <- sample_dirichlet_row(5000, probs, concentration = 10000)
  loose <- sample_dirichlet_row(5000, probs, concentration = 10)
  expect_true(var(tight[, 1]) < var(loose[, 1]))
})

test_that("sample_dirichlet_row rejects probs that don't sum to 1", {
  expect_error(sample_dirichlet_row(10, c(0.5, 0.4), concentration = 100))
})

test_that("run_psa is reproducible under a fixed seed and produces one row per (draw, arm)", {
  a <- run_psa(n_draws = 10, raw_dir = RAW_DIR, proc_dir = PROC_DIR, seed = 1)
  b <- run_psa(n_draws = 10, raw_dir = RAW_DIR, proc_dir = PROC_DIR, seed = 1)
  expect_equal(a, b)
  expect_equal(nrow(a), 10 * 4)
  expect_setequal(a$intervention, c("UST", "IFX", "ADA", "TREG"))
})

test_that("run_psa's comparator-arm hoist (2026-08-06) reproduces run_comparator_arm_lifetime() called directly, per draw", {
  # The hoist itself (this function's own docstring, "Performance note"): comparator arms are
  # simulated ONCE outside the draw loop, not once per draw, on the argument that nothing PSA
  # samples changes their occupancy trace. This test is the check that argument actually holds --
  # comparing run_psa()'s own output against run_comparator_arm_lifetime() called directly (the
  # UN-hoisted code path) for every draw, not just asserting the hoist runs without error.
  res <- run_psa(n_draws = 5, raw_dir = RAW_DIR, proc_dir = PROC_DIR, seed = 42)
  # Re-seed identically and re-draw the utility chain FIRST, exactly as run_psa() itself does
  # right after its own set.seed(42) call -- sample_utility_chain_draws() is the first
  # RNG-consuming call inside run_psa(), so this reproduces the same n_draws x 5 utility matrix.
  set.seed(42)
  utility_draws <- sample_utility_chain_draws(5, PROC_DIR)

  for (tx in c("UST", "IFX", "ADA")) {
    rows <- res[res$intervention == tx, ]
    rows <- rows[order(rows$draw), ]
    for (i in seq_len(5)) {
      expected <- run_comparator_arm_lifetime(
        tx, HORIZON_CYCLES_6YR, build_all_transition_matrices(RAW_DIR), raw_dir = RAW_DIR,
        proc_dir = PROC_DIR, utilities = utility_draws[i, ]
      )
      expect_equal(rows$qalys[i], expected$qalys, tolerance = 1e-9, info = paste(tx, i))
      expect_equal(rows$total_cost[i], expected$total_cost, tolerance = 1e-9, info = paste(tx, i))
    }
  }
})

test_that("run_psa's comparator arms have constant cost across draws (nothing sampled affects it yet) but varying QALYs (utility is sampled)", {
  res <- run_psa(n_draws = 15, raw_dir = RAW_DIR, proc_dir = PROC_DIR, seed = 2)
  ust <- res[res$intervention == "UST", ]
  expect_equal(length(unique(ust$total_cost)), 1)
  expect_true(length(unique(ust$qalys)) > 1)
  expect_true(all(is.na(ust$pi_sdr)) && all(is.na(ust$treg_price)))
})

test_that("run_psa's TREG rows carry sampled pi in [0,1], a positive price and a positive relapse hazard, all varying draw to draw", {
  res <- run_psa(n_draws = 15, raw_dir = RAW_DIR, proc_dir = PROC_DIR, seed = 3)
  treg <- res[res$intervention == "TREG", ]
  expect_true(all(treg$pi_sdr >= 0 & treg$pi_sdr <= 1))
  expect_true(all(treg$treg_price > 0))
  expect_true(all(treg$relapse_hazard_annual > 0))
  expect_true(length(unique(treg$pi_sdr)) > 1)
  expect_true(length(unique(treg$treg_price)) > 1)
  expect_true(length(unique(treg$relapse_hazard_annual)) > 1)

  # The comparator arms have no SDR state at all, so the hazard is meaningless for them and must
  # be recorded as NA rather than as a number that could be silently regressed against.
  comparators <- res[res$intervention != "TREG", ]
  expect_true(all(is.na(comparators$relapse_hazard_annual)))
})

test_that("run_psa's per-draw relapse hazard actually reaches the Treg simulation -- independent re-derivation, and it is NOT the pre-B3 h = 0 behaviour", {
  # The acceptance test for B3's sampling half. Two halves, both necessary:
  #
  #   (a) Re-run run_treg_arm_lifetime() OUTSIDE run_psa(), from that draw's own recorded
  #       (pi, price, utilities, h), and require it to reproduce run_psa()'s TREG row exactly. If
  #       the sampled hazard were being drawn but not passed (the precise shape of the original
  #       bug -- run_psa() sampled nothing and passed a hardcoded 0), this half fails.
  #   (b) Re-run the same draw with h = 0 and require the result to DIFFER materially. Without
  #       this half, (a) alone would still pass if run_psa() and the check both used 0, so (b) is
  #       what makes the test specific to the fix rather than merely self-consistent.
  res <- run_psa(n_draws = 5, raw_dir = RAW_DIR, proc_dir = PROC_DIR, seed = 77)
  # Reproduce the same RNG stream in the same order run_psa() consumes it (utility chain, price,
  # pi, then hazard last -- run_psa()'s own docstring on why hazard is drawn last).
  set.seed(77)
  utility_draws <- sample_utility_chain_draws(5, PROC_DIR)
  invisible(sample_treg_price_draws(5, PROC_DIR))
  invisible(stats::runif(5, 0, 1))
  hazard_draws <- sample_relapse_hazard_draws(5)

  treg <- res[res$intervention == "TREG", ]
  treg <- treg[order(treg$draw), ]
  expect_equal(treg$relapse_hazard_annual, hazard_draws, tolerance = 1e-12)

  matrices <- build_all_transition_matrices(RAW_DIR)
  for (i in seq_len(5)) {
    at_sampled_h <- run_treg_arm_lifetime(
      HORIZON_CYCLES_6YR, treg$pi_sdr[i], relapse_hazard_annual = treg$relapse_hazard_annual[i],
      price_usd = treg$treg_price[i], matrices = matrices, raw_dir = RAW_DIR, proc_dir = PROC_DIR,
      utilities = utility_draws[i, ]
    )
    expect_equal(treg$qalys[i], at_sampled_h$qalys, tolerance = 1e-9, info = paste("draw", i))
    expect_equal(treg$total_cost[i], at_sampled_h$total_cost, tolerance = 1e-9, info = paste("draw", i))

    at_zero_h <- run_treg_arm_lifetime(
      HORIZON_CYCLES_6YR, treg$pi_sdr[i], relapse_hazard_annual = 0,
      price_usd = treg$treg_price[i], matrices = matrices, raw_dir = RAW_DIR, proc_dir = PROC_DIR,
      utilities = utility_draws[i, ]
    )
    # Sampling h can only ever LOSE Treg QALYs relative to a permanent cure (SDR is never worse,
    # per QALY, than the non-cured track a relapse returns the patient to -- headroom_pi_star()'s
    # own monotonicity argument), so the direction is checkable, not just the magnitude. Draws
    # with pi ~ 0 have almost nobody in SDR and so barely move; the strict inequality is asserted
    # only where there is a real cured population to lose.
    expect_true(at_zero_h$qalys >= treg$qalys[i] - 1e-9, info = paste("draw", i))
    if (treg$pi_sdr[i] > 0.1) {
      expect_true(at_zero_h$qalys - treg$qalys[i] > 1e-4, info = paste("draw", i))
    }
  }
})

test_that("TREG's PSA QALY spread is not explained by the utility chain and pi alone -- the relapse hazard contributes real variance", {
  # (2) of B3's acceptance criteria: the sampled hazard must be a genuine per-draw driver of
  # Treg's outcome, not a recorded-but-inert column. Checked by re-simulating every draw at h = 0
  # while holding that draw's pi/price/utilities fixed: the difference between the two QALY
  # vectors is, by construction, the pure effect of the hazard draw. If h were inert that
  # difference would be identically zero.
  n <- 40
  res <- run_psa(n_draws = n, raw_dir = RAW_DIR, proc_dir = PROC_DIR, seed = 78)
  set.seed(78)
  utility_draws <- sample_utility_chain_draws(n, PROC_DIR)

  treg <- res[res$intervention == "TREG", ]
  treg <- treg[order(treg$draw), ]
  matrices <- build_all_transition_matrices(RAW_DIR)
  qalys_at_zero_h <- vapply(seq_len(n), function(i) {
    run_treg_arm_lifetime(
      HORIZON_CYCLES_6YR, treg$pi_sdr[i], relapse_hazard_annual = 0, price_usd = treg$treg_price[i],
      matrices = matrices, raw_dir = RAW_DIR, proc_dir = PROC_DIR, utilities = utility_draws[i, ]
    )$qalys
  }, numeric(1))

  hazard_effect <- qalys_at_zero_h - treg$qalys
  expect_true(all(hazard_effect >= -1e-9))          # never favours relapse over a permanent cure
  expect_true(sd(hazard_effect) > 1e-4)              # and it genuinely VARIES across draws
  expect_true(all(hazard_effect > 0))                # every draw loses something; none is inert

  # The effect must track the two things it structurally depends on. Its exact form is known, not
  # guessed: verify_pi_factorization() (R/05_deterministic_results.R) establishes that Treg's
  # incremental QALY over its own pi=0 track is EXACTLY pi * g(h) at any fixed h, so
  # hazard_effect_i = pi_i * (g_i(0) - g_i(h_i)). Dividing out pi therefore leaves a quantity that
  # is a monotone increasing function of h alone -- except that g_i also depends on that draw's own
  # sampled utility vector, which is why this is checked as a RANK correlation with a real
  # tolerance rather than as an identity. Observed 0.82-0.92 across four seeds; 0.75 is the floor.
  normalised_effect <- hazard_effect / treg$pi_sdr
  rank_cor_with_h <- stats::cor(normalised_effect, treg$relapse_hazard_annual, method = "spearman")
  expect_true(rank_cor_with_h > 0.75)

  # ...and the utility chain is NOT what is driving it, which is the specific alternative
  # explanation this test exists to rule out (before B3, utilities were the only thing that made
  # Treg's QALYs move at a given pi). The hazard's rank association must clearly beat the utility
  # chain's own association with the same effect.
  expect_true(rank_cor_with_h > abs(stats::cor(hazard_effect, treg$util_remission)) + 0.3)
})

test_that("psa_cost_effectiveness_plane hand-computed on a tiny toy PSA result", {
  toy <- data.frame(
    draw = c(1, 1, 2, 2),
    intervention = c("TREG", "IFX", "TREG", "IFX"),
    qalys = c(5, 4, 6, 4.5),
    total_cost = c(50000, 40000, 55000, 41000)
  )
  plane <- psa_cost_effectiveness_plane(toy, "TREG", "IFX")
  expect_equal(plane$delta_qalys, c(1, 1.5))
  expect_equal(plane$delta_cost, c(10000, 14000))
})

test_that("psa_ceac picks the arm with the highest NMB per draw and rows sum to 1 across arms", {
  toy <- data.frame(
    draw = rep(1:2, each = 2),
    intervention = rep(c("A", "B"), 2),
    qalys = c(4, 3, 4, 3),          # A always has more QALYs
    total_cost = c(100000, 0, 100000, 0)  # but B is free -- at low WTP, B wins; at high WTP, A wins
  )
  ceac <- psa_ceac(toy, wtp_grid = c(1000, 1000000))
  low_wtp <- ceac[ceac$wtp_usd == 1000, ]
  high_wtp <- ceac[ceac$wtp_usd == 1000000, ]
  expect_equal(low_wtp$B, 1)   # NMB_A = 4000-100000 <0 ; NMB_B = 3000 -> B wins
  expect_equal(high_wtp$A, 1)  # NMB_A = 4e6-1e5 ; NMB_B = 3e6 -> A wins
  expect_equal(low_wtp$A + low_wtp$B, 1)
  expect_equal(high_wtp$A + high_wtp$B, 1)
})
