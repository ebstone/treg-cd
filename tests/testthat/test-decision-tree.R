# testthat::test_dir() runs with the working directory set to tests/testthat/, but this
# project's scripts source() each other with repo-root-relative paths (README.md convention).
# Load in dependency order via repo-root-relative paths so each script's own
# `if (!exists(...))` guard skips its internal (otherwise-broken-from-here) source() call.
repo_root_relative <- function(...) file.path("..", "..", ...)
source(repo_root_relative("R", "utils", "transition_matrix.R"), local = TRUE)
source(repo_root_relative("R", "utils", "life_table.R"), local = TRUE)
source(repo_root_relative("R", "00_derive_transition_probs.R"), local = TRUE)
source(repo_root_relative("R", "01_decision_tree.R"), local = TRUE)
source(repo_root_relative("R", "02_markov_engine.R"), local = TRUE)

test_that("run_decision_tree splits UST's real induction row correctly and sums to exactly 1", {
  induction <- load_published_induction(repo_root_relative("data", "raw"))
  ust_row <- induction[induction$therapy == "UST", ]
  split <- run_decision_tree(ust_row)

  # Non-responders never start on the biologic.
  expect_equal(unname(split$initial_on_biologic["Moderate-Severe"]), 0)
  # ...they start directly on CT instead, at (normalised) to_moderate_severe.
  expect_equal(
    unname(split$initial_on_ct["Moderate-Severe"]),
    ust_row$to_moderate_severe / (ust_row$to_moderate_severe + ust_row$to_moderate_severe_responder +
      ust_row$to_mild + ust_row$to_remission),
    tolerance = 1e-9
  )
  # Everything else in initial_on_ct is 0; everything else non-zero in initial_on_biologic is
  # M-SR/Mild/Remission only.
  expect_equal(unname(split$initial_on_ct[c("Moderate-Severe Responder", "Mild", "Remission", "Surgery", "Death")]), rep(0, 5))
  expect_equal(unname(split$initial_on_biologic[c("Surgery", "Death")]), c(0, 0))

  expect_equal(sum(split$initial_on_biologic) + sum(split$initial_on_ct), 1, tolerance = 1e-12)
})

test_that("run_decision_tree hand-computed: exact split and normalisation arithmetic", {
  # Already sums to exactly 1 -- confirms the split itself, independent of normalisation.
  row_exact <- data.frame(
    therapy = "TEST", to_moderate_severe = 0.5, to_moderate_severe_responder = 0.2,
    to_mild = 0.2, to_remission = 0.1
  )
  split_exact <- run_decision_tree(row_exact)
  expect_equal(unname(split_exact$initial_on_ct["Moderate-Severe"]), 0.5)
  expect_equal(unname(split_exact$initial_on_biologic["Moderate-Severe Responder"]), 0.2)
  expect_equal(unname(split_exact$initial_on_biologic["Mild"]), 0.2)
  expect_equal(unname(split_exact$initial_on_biologic["Remission"]), 0.1)

  # Sums to 1.02 (publication-rounding-style overshoot) -- every component should be divided by
  # 1.02, not left as raw published proportions.
  row_over <- data.frame(
    therapy = "TEST", to_moderate_severe = 0.51, to_moderate_severe_responder = 0.204,
    to_mild = 0.204, to_remission = 0.102
  )
  split_over <- run_decision_tree(row_over)
  expect_equal(unname(split_over$initial_on_ct["Moderate-Severe"]), 0.5, tolerance = 1e-9)
  expect_equal(unname(split_over$initial_on_biologic["Remission"]), 0.1, tolerance = 1e-9)
  expect_equal(sum(split_over$initial_on_biologic) + sum(split_over$initial_on_ct), 1, tolerance = 1e-12)
})

test_that("run_decision_tree_all covers UST/IFX/ADA, each summing to exactly 1", {
  all_splits <- run_decision_tree_all(repo_root_relative("data", "raw"))
  expect_setequal(names(all_splits), c("UST", "IFX", "ADA"))
  for (tx in names(all_splits)) {
    total <- sum(all_splits[[tx]]$initial_on_biologic) + sum(all_splits[[tx]]$initial_on_ct)
    expect_equal(total, 1, tolerance = 1e-12, info = tx)
  }
})

test_that("end-to-end: R/00 induction -> R/01 decision tree -> R/02 maintenance engine conserves cohort for UST", {
  maintenance <- load_published_maintenance(repo_root_relative("data", "raw"))
  ust_m <- build_transition_matrix(maintenance[maintenance$therapy == "UST", ], MAINTENANCE_STATES)
  ct_m <- build_transition_matrix(maintenance[maintenance$therapy == "CT", ], MAINTENANCE_STATES)

  induction <- load_published_induction(repo_root_relative("data", "raw"))
  split <- run_decision_tree(induction[induction$therapy == "UST", ])

  res <- run_maintenance_arm(
    ust_m, ct_m, split$initial_on_biologic, n_cycles = 300,
    cap_cycle = 52, apply_cap = TRUE, initial_on_ct = split$initial_on_ct
  )

  # Full-pipeline cohort conservation, not just R/02 in isolation.
  expect_equal(rowSums(res$total), rep(1, 301), tolerance = 1e-9)

  # Non-responders land in on_ct immediately at cycle 0 (row 1) -- not just accumulating from
  # cycle 1 onward via the Moderate-Severe switch.
  expect_equal(unname(res$on_ct[1, "Moderate-Severe"]), unname(split$initial_on_ct["Moderate-Severe"]))
  expect_true(res$on_ct[1, "Moderate-Severe"] > 0)
  # Correspondingly, on_biologic starts with zero Moderate-Severe mass.
  expect_equal(unname(res$on_biologic[1, "Moderate-Severe"]), 0)
})
