# Tests for R/utils/refractory_multipliers.R and its wiring into R/05's run_comparator_arm_lifetime()/
# run_refractory_scenario() -- scenario S3 (analysis_plan.md §10.3, Decision 3), closed 2026-08-05.

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

RAW_DIR <- repo_root_relative("data", "raw")
PROC_DIR <- repo_root_relative("data", "processed")

test_that("load_refractory_multipliers reproduces the sourced UNITI-1/UNITI-2 rates exactly", {
  m <- load_refractory_multipliers(RAW_DIR)
  # Figure 1A/1B, week 8, 6 mg/kg (data/raw/uniti_refractory_population_multipliers.csv)
  expect_equal(m$induction_response_multiplier, 0.378 / 0.579)
  expect_equal(m$induction_remission_multiplier, 0.209 / 0.402)
  # Figure 3B, week 44, q8w
  expect_equal(m$maintenance_remission_multiplier_cumulative, 0.411 / 0.625)
  expect_equal(m$maintenance_cumulative_n_cycles, 18)
  # All four sourced rates are, correctly, refractory < naive (UNITI-1 < UNITI-2 in every pairing)
  expect_true(m$induction_response_multiplier < 1)
  expect_true(m$induction_remission_multiplier < 1)
  expect_true(m$maintenance_remission_multiplier_cumulative < 1)
})

test_that("load_refractory_multipliers reproduces the sourced surgery-hazard rates exactly (upper-bound proxy)", {
  m <- load_refractory_multipliers(RAW_DIR)
  # data/raw/refractory_surgery_hazard_multiplier.csv: SOJOURN 11.6% (naive) vs. Kassouri et al.
  # 2020 23.5% (double-biologic-failure, a MORE severe population than this project's refractory
  # scenario) -- a deliberate upper-bound proxy, documented as such in that file's own header.
  expect_equal(m$surgery_hazard_multiplier_cumulative, 0.235 / 0.116)
  expect_equal(m$surgery_cumulative_n_cycles, 26)
  # Surgery hazard is elevated (>1), the opposite direction from the response/remission
  # multipliers (<1) -- refractory disease makes surgery MORE likely, not less.
  expect_true(m$surgery_hazard_multiplier_cumulative > 1)
})

test_that("apply_refractory_multiplier_induction returns a valid probability row that sums to 1", {
  induction <- load_published_induction(RAW_DIR)
  mult <- load_refractory_multipliers(RAW_DIR)
  for (tx in COMPARATOR_THERAPIES) {
    row <- induction[induction$therapy == tx, ]
    adj <- apply_refractory_multiplier_induction(row, mult$induction_response_multiplier,
                                                  mult$induction_remission_multiplier)
    total <- adj$to_moderate_severe + adj$to_moderate_severe_responder + adj$to_mild + adj$to_remission
    expect_equal(total, 1, tolerance = 1e-9)
    expect_true(all(c(adj$to_moderate_severe, adj$to_moderate_severe_responder,
                       adj$to_mild, adj$to_remission) >= 0))
    # The refractory adjustment can only ever reduce remission and total response relative to the
    # naive row -- both multipliers are < 1 and remission is floor-clamped, never inflated.
    expect_true(adj$to_remission < row$to_remission)
    expect_true(adj$to_remission + adj$to_mild + adj$to_moderate_severe_responder <
                  row$to_remission + row$to_mild + row$to_moderate_severe_responder)
    expect_true(adj$to_moderate_severe > row$to_moderate_severe)
  }
})

test_that("apply_refractory_multiplier_induction's response floor-clamp never lets remission exceed total response", {
  # A deliberately adversarial row: remission_multiplier much closer to 1 than response_multiplier,
  # so a naive (unclamped) implementation would produce a negative "response-only" remainder.
  row <- data.frame(to_moderate_severe = 0.5, to_moderate_severe_responder = 0.05,
                     to_mild = 0.05, to_remission = 0.4)
  adj <- apply_refractory_multiplier_induction(row, response_multiplier = 0.1, remission_multiplier = 0.99)
  total <- adj$to_moderate_severe + adj$to_moderate_severe_responder + adj$to_mild + adj$to_remission
  expect_equal(total, 1, tolerance = 1e-9)
  expect_true(adj$to_mild >= 0 && adj$to_moderate_severe_responder >= 0)
  expect_true(adj$to_remission <= adj$to_remission + adj$to_mild + adj$to_moderate_severe_responder)
})

test_that("apply_refractory_multiplier_maintenance preserves row-stochasticity for every therapy's matrix", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  mult <- load_refractory_multipliers(RAW_DIR)
  for (tx in COMPARATOR_THERAPIES) {
    m <- matrices[[tx]]
    adj <- apply_refractory_multiplier_maintenance(m, mult$maintenance_remission_multiplier_cumulative,
                                                    mult$maintenance_cumulative_n_cycles)
    expect_equal(unname(rowSums(adj)), rep(1, nrow(adj)), tolerance = 1e-9)
    expect_true(all(adj >= -1e-12 & adj <= 1 + 1e-12))
    # Every row's Remission entry can only shrink (or stay at 0); Mild only grows (or stays put)
    # to absorb exactly the removed mass -- the redistribution moves mass one notch down the
    # severity axis, never anywhere else.
    expect_true(all(adj[, "Remission"] <= m[, "Remission"] + 1e-12))
    expect_true(all(adj[, "Mild"] >= m[, "Mild"] - 1e-12))
    other_cols <- setdiff(colnames(m), c("Remission", "Mild"))
    expect_equal(adj[, other_cols], m[, other_cols])
  }
})

test_that("apply_refractory_multiplier_maintenance's surgery_hazard_multiplier_cumulative = NULL default is byte-identical to omitting it", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  mult <- load_refractory_multipliers(RAW_DIR)
  for (tx in COMPARATOR_THERAPIES) {
    m <- matrices[[tx]]
    without_arg <- apply_refractory_multiplier_maintenance(m, mult$maintenance_remission_multiplier_cumulative,
                                                            mult$maintenance_cumulative_n_cycles)
    explicit_null <- apply_refractory_multiplier_maintenance(m, mult$maintenance_remission_multiplier_cumulative,
                                                              mult$maintenance_cumulative_n_cycles,
                                                              surgery_hazard_multiplier_cumulative = NULL,
                                                              surgery_cumulative_n_cycles = NULL)
    expect_equal(explicit_null, without_arg)
  }
})

test_that("apply_refractory_multiplier_maintenance's surgery-hazard step preserves row-stochasticity and only ever raises Surgery", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  mult <- load_refractory_multipliers(RAW_DIR)
  for (tx in COMPARATOR_THERAPIES) {
    m <- matrices[[tx]]
    adj <- apply_refractory_multiplier_maintenance(
      m, mult$maintenance_remission_multiplier_cumulative, mult$maintenance_cumulative_n_cycles,
      surgery_hazard_multiplier_cumulative = mult$surgery_hazard_multiplier_cumulative,
      surgery_cumulative_n_cycles = mult$surgery_cumulative_n_cycles
    )
    expect_equal(unname(rowSums(adj)), rep(1, nrow(adj)), tolerance = 1e-9)
    expect_true(all(adj >= -1e-12 & adj <= 1 + 1e-12))
    # Surgery can only rise (or stay at 0, for rows with no baseline surgery probability at all) --
    # the multiplier is > 1 by construction (surgery_hazard_multiplier_cumulative > 1).
    expect_true(all(adj[, "Surgery"] >= m[, "Surgery"] - 1e-12))
    # Rows with zero baseline Surgery probability are untouched by the surgery step specifically
    # (0 * multiplier = 0 -- nothing to elevate, no rescale triggered).
    zero_surgery_rows <- which(m[, "Surgery"] <= 0)
    if (length(zero_surgery_rows) > 0) {
      expect_equal(adj[zero_surgery_rows, "Surgery"], m[zero_surgery_rows, "Surgery"])
    }
  }
})

test_that("apply_refractory_multiplier_maintenance requires surgery_cumulative_n_cycles whenever a surgery multiplier is passed", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  mult <- load_refractory_multipliers(RAW_DIR)
  m <- matrices[["UST"]]
  expect_error(
    apply_refractory_multiplier_maintenance(
      m, mult$maintenance_remission_multiplier_cumulative, mult$maintenance_cumulative_n_cycles,
      surgery_hazard_multiplier_cumulative = mult$surgery_hazard_multiplier_cumulative
      # surgery_cumulative_n_cycles omitted -- should error, not silently divide by NULL
    )
  )
})

test_that("run_comparator_arm_lifetime(refractory = FALSE) is byte-identical to the pre-existing default", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  for (tx in COMPARATOR_THERAPIES) {
    plain <- run_comparator_arm_lifetime(tx, HORIZON_CYCLES_6YR, matrices, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
    explicit_false <- run_comparator_arm_lifetime(tx, HORIZON_CYCLES_6YR, matrices, raw_dir = RAW_DIR,
                                                   proc_dir = PROC_DIR, refractory = FALSE)
    expect_equal(explicit_false, plain)
  }
})

test_that("run_comparator_arm_lifetime(refractory = TRUE) runs cleanly and reduces QALYs relative to naive", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  for (tx in COMPARATOR_THERAPIES) {
    naive <- run_comparator_arm_lifetime(tx, HORIZON_CYCLES_6YR, matrices, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
    refractory <- run_comparator_arm_lifetime(tx, HORIZON_CYCLES_6YR, matrices, raw_dir = RAW_DIR,
                                               proc_dir = PROC_DIR, refractory = TRUE)
    expect_false(anyNA(unlist(refractory)))
    expect_true(refractory$qalys > 0)
    expect_true(refractory$total_cost > 0)
    # A refractory population responds worse to every biologic -- QALYs must fall, strictly, for
    # all three comparators (each has a nonzero remission/response probability to begin with).
    expect_true(refractory$qalys < naive$qalys)
  }
})

test_that("run_refractory_scenario's biologic_naive rows exactly reproduce run_base_case()", {
  scenario <- run_refractory_scenario(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  base <- run_base_case(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  naive_rows <- scenario[scenario$population == "biologic_naive", setdiff(names(scenario), "population")]
  rownames(naive_rows) <- NULL
  rownames(base) <- NULL
  expect_equal(naive_rows, base)
})

test_that("run_refractory_scenario leaves TREG identical across both population rows", {
  scenario <- run_refractory_scenario(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  treg_rows <- scenario[scenario$intervention == "TREG", ]
  expect_equal(nrow(treg_rows), 2)
  expect_equal(treg_rows$qalys[1], treg_rows$qalys[2])
  expect_equal(treg_rows$total_cost[1], treg_rows$total_cost[2])
})

test_that("run_refractory_scenario's refractory rows have strictly lower QALYs than biologic_naive for every comparator", {
  scenario <- run_refractory_scenario(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  for (tx in COMPARATOR_THERAPIES) {
    naive_qalys <- scenario$qalys[scenario$population == "biologic_naive" & scenario$intervention == tx]
    refractory_qalys <- scenario$qalys[scenario$population == "refractory" & scenario$intervention == tx]
    expect_true(refractory_qalys < naive_qalys)
  }
})
