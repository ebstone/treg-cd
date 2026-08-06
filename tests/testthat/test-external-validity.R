# External face-validity tests (analysis_plan.md Aim 5; peer review 2026-08-05, M9): every other
# test in this suite checks INTERNAL consistency (rows sum to 1, cohort mass is conserved,
# round-trips agree with each other) -- none of them can catch a model that is internally perfect
# and externally wrong, which is exactly what happened with the induction-cycle bug (B1,
# docs/model_audit_v6.md A17). This file asserts reproduction of PUBLISHED quantities external to
# this codebase, not internal invariants, per the review's own explicit recommendation.
#
# testthat::test_dir() runs with the working directory set to tests/testthat/, but this project's
# scripts source() each other with repo-root-relative paths (README.md convention). Load in
# dependency order via repo-root-relative paths so each script's own `if (!exists(...))` guard
# skips its internal (otherwise-broken-from-here) source() call.
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

test_that("B1 acceptance test: end-of-induction Remission occupancy reproduces the published trial endpoint", {
  # Aliyev et al. 2019 Appendix S2's own worked example derives Table 3's UST Moderate-Severe ->
  # Remission entry (0.133) FROM the published UNITI week-6 remission proportion (0.349) via a
  # rate conversion -- so applying Table 3's 2-week matrix for the right number of cycles must
  # come back close to 0.349, not the raw 0.133 single-cycle figure. Similarly, IFX/ADA's 2-week
  # entry (0.203) is derived from a week-4 endpoint; CLASSIC-I/CHARM-style ADA induction data
  # (data/raw/aliyev2019_appendixS1_table2_parameters.csv: 0.50 week-4 response probability x 0.72
  # remitter:responder ratio) gives an external week-4 remission estimate of ~0.36. Tolerance
  # bands are the review's own explicit acceptance range (docs/model_audit_v6.md A17), not
  # machine-precision equality -- this is a face-validity check against an external trial
  # endpoint, not a round-trip of this codebase's own arithmetic against itself.
  out <- load_published_induction(RAW_DIR)

  ust <- out[out$therapy == "UST", ]
  expect_equal(ust$induction_cycles, 3)
  expect_equal(ust$to_remission, 0.349, tolerance = 0.03)  # published UNITI week-6 remission

  for (tx in c("IFX", "ADA")) {
    row <- out[out$therapy == tx, ]
    expect_equal(row$induction_cycles, 2)
    expect_equal(row$to_remission, 0.36, tolerance = 0.02, info = tx)  # ~week-4 remission estimate
  }

  # The bug this test exists to catch: reading Table 3's Moderate-Severe row as a 1-cycle terminal
  # split (the pre-B1-fix behaviour) understates every therapy's true end-of-induction Remission
  # occupancy by roughly half -- assert the corrected value is well above the raw, uncorrected
  # 2-week entry for every therapy, not just close to the external endpoint.
  raw_2wk_remission <- c(UST = 0.133, IFX = 0.203, ADA = 0.203)
  for (tx in out$therapy) {
    row <- out[out$therapy == tx, ]
    expect_true(row$to_remission > raw_2wk_remission[[tx]] * 1.5, info = tx)
  }
})

test_that("induction terminal occupancy conserves cohort mass (external check should not come at the cost of the internal one)", {
  out <- load_published_induction(RAW_DIR)
  totals <- out$to_moderate_severe + out$to_moderate_severe_responder + out$to_mild + out$to_remission
  expect_equal(totals, rep(1, nrow(out)), tolerance = 1e-9)
})

test_that("B3 acceptance test: the PSA's durability prior is centred inside PolTREG's own published observation window", {
  # Added 2026-08-06 with B3 (docs/model_audit_v6.md A19). The post-cure relapse hazard h has no
  # trial of its own to validate against -- no Treg product has ever been followed to relapse in
  # Crohn's disease -- so the strongest external anchor available is the one
  # R/05_deterministic_results.R's RELAPSE_DURATION_GRID_YEARS comment already cites: the PolTREG
  # T1D cohort, 54 patients from Phase I/II trials followed 7-12 years, a subset of whom remained
  # in clinical remission across that window (analysis_plan.md §7.2's analog table, recorded there
  # as company/conference communications rather than a peer-reviewed long-term publication -- so
  # this is a plausibility band, deliberately not a point estimate).
  #
  # The check is therefore a band check, not an equality: the CENTRE of the sampled durability
  # distribution, expressed as a median SDR duration in years, must land inside the 7-12 year
  # window that cohort was actually observed over. This is external in the sense that matters --
  # it fails if someone re-centres the prior on a durability nobody has ever observed (either
  # direction: a 3-year median would contradict PolTREG's plateau, a 40-year median would be a
  # claim about a period no Treg product has existed for), which no internal round-trip or
  # sum-to-one invariant in this suite could catch.
  set.seed(2026)
  draws <- sample_relapse_hazard_draws(20000)
  median_implied_duration <- hazard_to_duration(stats::median(draws))
  expect_true(median_implied_duration >= 7 && median_implied_duration <= 12,
              info = sprintf("median implied SDR duration = %.2f years", median_implied_duration))

  # The deterministic base case must sit in that same window, and for the same reason -- it is the
  # PSA prior's mean, not an independently chosen figure (sample_relapse_hazard_draws()'s own
  # docstring), so the two cannot drift apart without this failing.
  expect_true(DEFAULT_RELAPSE_DURATION_YEARS >= 7 && DEFAULT_RELAPSE_DURATION_YEARS <= 12)
  expect_equal(mean(draws), duration_to_hazard(DEFAULT_RELAPSE_DURATION_YEARS), tolerance = 0.03)

  # Permanent remission is NOT supported by any published source (that comment's own words:
  # "nothing published supports permanence"), so the prior must not put meaningful mass on
  # durations beyond anything ever observed for any cell therapy in any indication.
  expect_true(mean(vapply(draws, hazard_to_duration, numeric(1)) > 100) < 0.03)
})
