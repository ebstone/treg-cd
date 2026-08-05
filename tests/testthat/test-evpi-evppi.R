# testthat::test_dir() runs with the working directory set to tests/testthat/, but this
# project's scripts source() each other with repo-root-relative paths (README.md convention).
# Same loading pattern as test-psa.R.
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
source(repo_root_relative("R", "06_psa.R"), local = TRUE)
source(repo_root_relative("R", "07_evpi_evppi.R"), local = TRUE)

RAW_DIR <- repo_root_relative("data", "raw")
PROC_DIR <- repo_root_relative("data", "processed")

test_that("evpi_from_nb is a simple, hand-computed formula", {
  # Draw 1: A wins (10 > 4); draw 2: B wins (8 > 2). Perfect foresight picks the winner each time
  # (mean 9); current information picks whichever has the higher AVERAGE (A: mean 6, B: mean 5).
  nb <- cbind(A = c(10, 2), B = c(4, 8))
  expect_equal(evpi_from_nb(nb), mean(c(10, 8)) - max(c(mean(c(10, 2)), mean(c(4, 8)))))
  expect_equal(evpi_from_nb(nb), 3)
})

test_that("evpi_from_nb is exactly 0 when one arm dominates in every draw", {
  nb <- cbind(A = c(10, 12, 9), B = c(5, 6, 4))
  expect_equal(evpi_from_nb(nb), 0)
})

test_that("net_benefit_matrix pivots correctly and computes NMB = wtp*qalys - cost", {
  toy <- data.frame(
    draw = c(1, 1, 2, 2), intervention = c("TREG", "IFX", "TREG", "IFX"),
    qalys = c(5, 4, 6, 4.5), total_cost = c(50000, 40000, 55000, 41000)
  )
  nb <- net_benefit_matrix(toy, wtp_usd = 100000)
  expect_equal(unname(nb[, "TREG"]), 100000 * c(5, 6) - c(50000, 55000))
  expect_equal(unname(nb[, "IFX"]), 100000 * c(4, 4.5) - c(40000, 41000))
})

test_that("net_benefit_matrix's treg_price_override exploits linearity exactly (algebraic swap, no re-simulation)", {
  toy <- data.frame(
    draw = c(1, 2), intervention = c("TREG", "TREG"),
    qalys = c(5, 6), total_cost = c(50000, 55000), treg_price = c(20000, 22000)
  )
  nb <- net_benefit_matrix(toy, wtp_usd = 100000, treg_price_override = 10000)
  # cost should drop by exactly (sampled_price - override) for each draw
  expected_cost <- toy$total_cost - toy$treg_price + 10000
  expect_equal(unname(nb[, "TREG"]), 100000 * toy$qalys - expected_cost)
})

test_that("pca_reduce recovers (almost) all variance in 1 component when the data is exactly rank-1", {
  set.seed(1)
  z <- rnorm(500)
  m <- cbind(a = z, b = 2 * z, c = -0.5 * z)  # every column is a scalar multiple of the same z
  reduced <- pca_reduce(m, n_components = 1)
  expect_equal(ncol(reduced), 1)
  pr <- stats::prcomp(m, center = TRUE, scale. = TRUE)
  expect_equal(pr$sdev[1]^2 / sum(pr$sdev^2), 1, tolerance = 1e-6)
})

test_that("reduce_for_gam passes small matrices through untouched and only reduces when needed", {
  small <- cbind(a = 1:5, b = 5:1)
  r <- reduce_for_gam(small, max_raw_dims = 2)
  expect_equal(r$matrix, small)
  expect_equal(r$n_raw, 2)
  expect_equal(r$n_used, 2)
  expect_equal(r$variance_explained, 1)

  set.seed(2)
  big <- matrix(rnorm(500 * 5), ncol = 5)
  r2 <- reduce_for_gam(big, max_raw_dims = 2, n_components = 2)
  expect_equal(r2$n_raw, 5)
  expect_equal(r2$n_used, 2)
  expect_true(r2$variance_explained > 0 && r2$variance_explained <= 1)
  expect_equal(ncol(r2$matrix), 2)
})

test_that("evppi_gam recovers close to the full EVPI when one arm's benefit is (near-)deterministic in the parameter", {
  set.seed(3)
  n <- 1500
  x <- stats::runif(n)
  noise <- stats::rnorm(n, 0, 0.02)
  nb <- cbind(A = 10 * x + noise, B = 5 + noise)  # A's benefit is essentially determined by x; B is flat
  full_evpi <- evpi_from_nb(nb)
  est <- evppi_gam(nb, cbind(x = x))
  expect_equal(est, full_evpi, tolerance = 0.05)
})

test_that("evppi_gam is close to 0 when net benefit doesn't actually depend on the parameter", {
  set.seed(4)
  n <- 1500
  x <- stats::runif(n)                    # a parameter that's pure noise w.r.t. the outcome
  nb <- cbind(A = stats::rnorm(n, 10, 1), B = stats::rnorm(n, 9, 1))  # A dominates regardless of x
  est <- evppi_gam(nb, cbind(x = x))
  expect_equal(est, 0, tolerance = 0.05)
})

test_that("population_evpi is a simple multiply and rejects a negative population", {
  expect_equal(population_evpi(10, 1000), 10000)
  expect_equal(population_evpi(0, 1000), 0)
  expect_error(population_evpi(10, -5))
})

test_that("cross_check_voi returns NULL with a message when voi isn't installed (true in this environment)", {
  skip_if(requireNamespace("voi", quietly = TRUE), "voi is actually installed here; cross-check path not exercised by this test")
  toy <- data.frame(draw = 1:5, intervention = "TREG", qalys = 1:5, total_cost = 1:5 * 100)
  expect_message(result <- cross_check_voi(toy, cbind(x = 1:5), 100000), "voi is not installed")
  expect_null(result)
})

# ---- Integration tests against a real (small) PSA run -------------------------------------------

test_that("evpi_surface returns one row per (price, wtp) grid point, all EVPI >= 0", {
  psa <- run_psa(n_draws = 200, raw_dir = RAW_DIR, proc_dir = PROC_DIR, seed = 11)
  price_grid <- c(5000, 15000, 25000)
  wtp_grid <- c(50000, 150000)
  surf <- evpi_surface(psa, price_grid, wtp_grid)
  expect_equal(nrow(surf), length(price_grid) * length(wtp_grid))
  expect_true(all(surf$evpi_per_patient >= -1e-6))  # allow GAM/float noise around exact 0
})

test_that("evpi_surface's price override matches a manual net_benefit_matrix computation at the same price", {
  psa <- run_psa(n_draws = 150, raw_dir = RAW_DIR, proc_dir = PROC_DIR, seed = 12)
  surf <- evpi_surface(psa, price_grid_usd = 12000, wtp_grid_usd = 100000)
  manual_nb <- net_benefit_matrix(psa, 100000, treg_price_override = 12000)
  expect_equal(surf$evpi_per_patient[1], evpi_from_nb(manual_nb))
})

test_that("evppi_by_subset returns all 7 expected subsets with sane, non-negative-ish values and correct dimension counts", {
  psa <- run_psa(n_draws = 250, raw_dir = RAW_DIR, proc_dir = PROC_DIR, seed = 13)
  res <- evppi_by_subset(psa, wtp_usd = 100000)

  expect_setequal(res$subset, c("A", "C", "E", "A u C", "A u E", "C u E", "A u C u E"))
  expect_true(all(res$evppi >= -1e-6))                        # allow GAM estimation noise near 0
  expect_true(all(res$evppi <= res$total_evpi + 1e-6))        # EVPPI can't exceed total EVPI
  expect_equal(res$n_raw_params[res$subset == "A"], 1)
  expect_equal(res$n_raw_params[res$subset == "E"], 5)
  expect_equal(res$n_params_used[res$subset == "E"], 2)       # PCA-reduced
  expect_equal(res$n_params_used[res$subset == "A u C"], 2)   # not reduced (<= max_raw_dims)
  expect_true(all(res$total_evpi == res$total_evpi[1]))       # same reference EVPI on every row
})

test_that("evppi_convergence returns one estimate per requested draw count, all finite", {
  psa <- run_psa(n_draws = 300, raw_dir = RAW_DIR, proc_dir = PROC_DIR, seed = 14)
  conv <- evppi_convergence(psa, "pi_sdr", wtp_usd = 100000, draw_counts = c(100, 200, 300))
  expect_equal(nrow(conv), 3)
  expect_equal(conv$n_draws, c(100, 200, 300))
  expect_true(all(is.finite(conv$evppi)))
})

test_that("evppi_convergence rejects a draw count larger than what's in psa_results", {
  psa <- run_psa(n_draws = 50, raw_dir = RAW_DIR, proc_dir = PROC_DIR, seed = 15)
  expect_error(evppi_convergence(psa, "pi_sdr", 100000, draw_counts = c(100)))
})

# ---- pi prior-sensitivity on EVPPI (resolutions memo sec 3.2) -----------------------------------

test_that("evpi_from_nb with weights = NULL is byte-identical to the original unweighted formula", {
  set.seed(20)
  nb <- cbind(A = rnorm(200, 10, 2), B = rnorm(200, 9, 2))
  expect_equal(evpi_from_nb(nb, weights = NULL), evpi_from_nb(nb))
  expect_equal(evpi_from_nb(nb, weights = rep(1, 200)), evpi_from_nb(nb))
})

test_that("evpi_from_nb's weights argument matches a hand-computed weighted expectation", {
  # Two draws, weight draw 2 three times as heavily as draw 1.
  nb <- cbind(A = c(10, 2), B = c(4, 8))
  w <- c(1, 3)
  wn <- w / sum(w)
  expected <- sum(c(10, 8) * wn) - max(sum(c(10, 2) * wn), sum(c(4, 8) * wn))
  expect_equal(evpi_from_nb(nb, weights = w), expected)
})

test_that("evpi_from_nb rejects a weights vector of the wrong length or with a negative entry", {
  nb <- cbind(A = c(10, 2), B = c(4, 8))
  expect_error(evpi_from_nb(nb, weights = c(1, 1, 1)))
  expect_error(evpi_from_nb(nb, weights = c(1, -1)))
})

test_that("prior_reweight under the reference U(0,1) density itself returns all-equal weights (no reweighting)", {
  set.seed(21)
  pi_values <- runif(500)
  w <- prior_reweight(pi_values, function(pi) stats::dunif(pi, 0, 1))
  expect_true(all(abs(w - 1) < 1e-10))
})

test_that("prior_reweight's weights have mean 1 and are strictly positive for a proper alternative density", {
  set.seed(22)
  pi_values <- runif(500)
  w_beta13 <- prior_reweight(pi_values, function(pi) stats::dbeta(pi, 1, 3))
  expect_equal(mean(w_beta13), 1, tolerance = 1e-10)
  expect_true(all(w_beta13 >= 0))
  # Beta(1,3) puts more mass near 0: low-pi draws should get upweighted relative to high-pi draws.
  low <- w_beta13[pi_values < 0.2]
  high <- w_beta13[pi_values > 0.8]
  expect_true(mean(low) > mean(high))
})

test_that("prior_reweight rejects a density that goes negative or non-finite anywhere it's evaluated", {
  pi_values <- c(0.1, 0.5, 0.9)
  expect_error(prior_reweight(pi_values, function(pi) -pi))
  expect_error(prior_reweight(pi_values, function(pi) rep(NA_real_, length(pi))))
})

test_that("evppi_gam with weights = NULL reproduces the unweighted estimate exactly", {
  set.seed(23)
  n <- 800
  x <- runif(n)
  nb <- cbind(A = 10 * x + rnorm(n, 0, 0.02), B = 5 + rnorm(n, 0, 0.02))
  est_default <- evppi_gam(nb, cbind(x = x))
  est_explicit_null <- evppi_gam(nb, cbind(x = x), weights = NULL)
  expect_equal(est_default, est_explicit_null)
})

test_that("evppi_gam with uniform weights reproduces the unweighted estimate", {
  set.seed(24)
  n <- 800
  x <- runif(n)
  nb <- cbind(A = 10 * x + rnorm(n, 0, 0.02), B = 5 + rnorm(n, 0, 0.02))
  est_unweighted <- evppi_gam(nb, cbind(x = x))
  est_uniform_weighted <- evppi_gam(nb, cbind(x = x), weights = rep(2, n))  # any constant, not just 1
  expect_equal(est_uniform_weighted, est_unweighted, tolerance = 1e-6)
})

test_that("evppi_by_subset's weights argument is threaded consistently, and shrinks subset A's EVPPI when pi's high-spread region is downweighted", {
  # A controlled synthetic PSA rather than a real (noisy, often near-zero-EVPI at small n) Markov
  # run: TREG's net benefit is constructed to depend STRONGLY and near-linearly on pi (mimicking
  # the headroom relationship -- higher cure fraction, better NMB), comparator flat aside from
  # independent noise, so essentially all decision-relevant uncertainty is attributable to subset
  # A (pi) by construction -- a regime where reweighting's effect is large enough to see clearly,
  # not lost in GAM estimation noise.
  set.seed(25)
  n <- 1500
  pi_vals <- stats::runif(n)
  treg_nmb <- 20000 * pi_vals + stats::rnorm(n, 0, 200)
  comp_nmb <- rep(10000, n) + stats::rnorm(n, 0, 200)
  wtp <- 100000
  qalys_common <- 5
  # util_* columns get tiny per-draw jitter (not a constant) purely so reduce_for_gam()'s PCA step
  # on subset E has non-zero variance to work with -- they don't otherwise enter this synthetic
  # nb's construction at all, so they contribute ~0 to any subset's EVPPI either way.
  jitter <- function(center) rep(center, n) + stats::rnorm(n, 0, 1e-6)
  # treg_price also gets jitter, not a literal constant: subsets that pool it with the util
  # columns (e.g. "C u E") PCA-reduce their raw matrix (reduce_for_gam()), and prcomp() itself
  # refuses to scale a genuinely zero-variance column -- a real constraint of that step, not
  # something this synthetic fixture should paper over by avoiding those subsets.
  psa_toy <- data.frame(
    draw = rep(1:n, 2), intervention = rep(c("TREG", "COMP"), each = n), qalys = qalys_common,
    total_cost = c(wtp * qalys_common - treg_nmb, wtp * qalys_common - comp_nmb),
    pi_sdr = c(pi_vals, rep(NA_real_, n)), treg_price = c(jitter(15000), rep(NA_real_, n)),
    util_modsev = rep(jitter(0.5), 2), util_resp = rep(jitter(0.6), 2),
    util_mild = rep(jitter(0.7), 2), util_remission = rep(jitter(0.8), 2),
    util_surgery = rep(jitter(0.4), 2),
    stringsAsFactors = FALSE
  )

  res_unweighted <- evppi_by_subset(psa_toy, wtp_usd = wtp)
  # Beta(1,3) puts most of its mass below pi ~ 0.3 -- downweighting the high-pi draws that drive
  # most of TREG's NMB spread should shrink subset A's estimated EVPPI (and the total EVPI it's
  # nearly all of) relative to the reference U(0,1) prior.
  w <- prior_reweight(pi_vals, function(pi) stats::dbeta(pi, 1, 3))
  res_weighted <- evppi_by_subset(psa_toy, wtp_usd = wtp, weights = w)

  expect_setequal(res_weighted$subset, res_unweighted$subset)
  expect_true(all(res_weighted$total_evpi == res_weighted$total_evpi[1]))  # internally consistent
  expect_true(res_weighted$evppi[res_weighted$subset == "A"] < res_unweighted$evppi[res_unweighted$subset == "A"])
  expect_true(res_weighted$total_evpi[1] < res_unweighted$total_evpi[1])
})

test_that("evppi_prior_sensitivity returns one full subset table per prior, with a stable subset set and a within-prior rank", {
  psa <- run_psa(n_draws = 300, raw_dir = RAW_DIR, proc_dir = PROC_DIR, seed = 17)
  res <- evppi_prior_sensitivity(psa, wtp_usd = 100000)

  expect_setequal(res$prior, names(PI_PRIOR_SENSITIVITY_SPECS))
  expect_equal(nrow(res), length(PI_PRIOR_SENSITIVITY_SPECS) * 7)  # 7 subsets, per evppi_by_subset()

  # Every prior's rows should carry the SAME 7 subsets.
  for (p in unique(res$prior)) {
    expect_setequal(res$subset[res$prior == p], c("A", "C", "E", "A u C", "A u E", "C u E", "A u C u E"))
  }

  # rank_within_prior is an integer rank within each prior group (ties -> shared minimum rank, so
  # small GAM-noise ties near 0 don't produce non-integer ranks), highest EVPPI = rank 1.
  for (p in unique(res$prior)) {
    sub <- res[res$prior == p, ]
    expect_true(all(sub$rank_within_prior == round(sub$rank_within_prior)))
    expect_true(min(sub$rank_within_prior) == 1)
    expect_true(max(sub$rank_within_prior) <= nrow(sub))
    expect_true(sub$rank_within_prior[which.max(sub$evppi)] == 1)
  }

  # The reference U(0,1) row should reproduce evppi_by_subset()'s own unweighted result exactly --
  # prior_reweight() under its own generating density returns all-1 weights (tested above), so
  # this prior's pass through evppi_by_subset() must be numerically identical to the plain call.
  reference <- res[res$prior == "U(0,1) [reference]", ]
  reference <- reference[order(reference$subset), ]
  plain <- evppi_by_subset(psa, wtp_usd = 100000)
  plain <- plain[order(plain$subset), ]
  expect_equal(reference$evppi, plain$evppi, tolerance = 1e-6)
})
