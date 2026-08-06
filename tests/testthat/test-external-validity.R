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
source(repo_root_relative("R", "00_derive_transition_probs.R"), local = TRUE)

RAW_DIR <- repo_root_relative("data", "raw")

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
