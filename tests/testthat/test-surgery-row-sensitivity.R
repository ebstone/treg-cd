# Tests for R/utils/surgery_row_sensitivity.R and its wiring into R/09_scenarios.R's
# run_scenario_r2_surgery_sensitivity() -- the R2 sensitivity check (peer review 2026-08-05,
# docs/model_audit_v6.md A17's own follow-up note), added 2026-08-06.

repo_root_relative <- function(...) file.path("..", "..", ...)
source(repo_root_relative("R", "utils", "transition_matrix.R"), local = TRUE)
source(repo_root_relative("R", "utils", "life_table.R"), local = TRUE)
source(repo_root_relative("R", "utils", "population_utility.R"), local = TRUE)
source(repo_root_relative("R", "utils", "refractory_multipliers.R"), local = TRUE)
source(repo_root_relative("R", "utils", "surgery_row_sensitivity.R"), local = TRUE)
source(repo_root_relative("R", "00_derive_transition_probs.R"), local = TRUE)
source(repo_root_relative("R", "01_decision_tree.R"), local = TRUE)
source(repo_root_relative("R", "02_markov_engine.R"), local = TRUE)
source(repo_root_relative("R", "03_cure_fraction_module.R"), local = TRUE)
source(repo_root_relative("R", "04_costs_utilities.R"), local = TRUE)
source(repo_root_relative("R", "05_deterministic_results.R"), local = TRUE)
source(repo_root_relative("R", "09_scenarios.R"), local = TRUE)

RAW_DIR <- repo_root_relative("data", "raw")
PROC_DIR <- repo_root_relative("data", "processed")

test_that("build_table2_surgery_row sums to exactly 1 and independently reproduces the hand-derived figures", {
  row <- build_table2_surgery_row()
  expect_equal(sum(row), 1, tolerance = 1e-9)
  expect_setequal(names(row), MAINTENANCE_STATES)

  # Independent re-derivation, NOT calling build_table2_surgery_row()'s own machinery, so this is
  # a real cross-check rather than a tautology. Table 2's raw 8-week entries: to M-S 0.060, to
  # Mild/M-SR combined 0.075, self-loop 0.34; residual (unreported) 1 - 0.060 - 0.075 - 0.34 =
  # 0.525 assumed to Remission (this file's own header, assumption 1).
  leave_8wk <- 1 - 0.34
  rate <- -log(1 - leave_8wk) / (8 / 52)
  leave_2wk <- 1 - exp(-rate * 2 / 52)
  expect_equal(unname(row[["Surgery"]]), 1 - leave_2wk, tolerance = 1e-9)
  expect_equal(unname(row[["Moderate-Severe"]]), leave_2wk * (0.060 / leave_8wk), tolerance = 1e-9)
  expect_equal(unname(row[["Remission"]]), leave_2wk * (0.525 / leave_8wk), tolerance = 1e-9)
  # Assumption 2: Mild/M-SR combined (0.075) split evenly
  expect_equal(unname(row[["Mild"]]), unname(row[["Moderate-Severe Responder"]]))
  expect_equal(unname(row[["Mild"]]), leave_2wk * (0.075 / leave_8wk) / 2, tolerance = 1e-9)
  # Assumption 3: no sourced Surgery-to-Death figure at any cycle length
  expect_equal(unname(row[["Death"]]), 0)
})

test_that("build_table2_surgery_row implies a MUCH stickier Surgery state than Table 4's own sourced row -- the whole point of this sensitivity check", {
  row <- build_table2_surgery_row()
  matrices <- build_all_transition_matrices(RAW_DIR)
  table4_self_loop <- matrices[["UST"]]["Surgery", "Surgery"]
  # Table 4's own 2-week self-loop is ~0.0976; Table 2's own 8-week self-loop (0.34) converts to a
  # 2-week-equivalent self-loop around 0.76 -- an order of magnitude apart (R2, project memory
  # verification note). Asserting "much larger", not a specific value, since this is the
  # qualitative claim the whole sensitivity check rests on.
  expect_true(row[["Surgery"]] > 5 * table4_self_loop)
})

test_that("apply_table2_surgery_row replaces ONLY the Surgery row, leaves every other row byte-identical, and never mutates its input", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  before <- matrices[["UST"]]
  altered <- apply_table2_surgery_row(matrices[["UST"]])

  expect_equal(altered["Surgery", ], build_table2_surgery_row())
  for (state in setdiff(MAINTENANCE_STATES, "Surgery")) {
    expect_identical(altered[state, ], before[state, ])
  }
  # The real UST matrix (as build_all_transition_matrices() returns it) is untouched -- this
  # function must not mutate its argument in place.
  expect_identical(matrices[["UST"]], before)
})

test_that("run_scenario_r2_surgery_sensitivity's table4 rows exactly reproduce run_scenario_s12_non_cured_hr's own base-case finding", {
  s12 <- run_scenario_s12_non_cured_hr(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  r2 <- run_scenario_r2_surgery_sensitivity(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  table4_rows <- r2[r2$surgery_source == "table4_2wk_sourced", ]

  for (hr in NON_CURED_HR_GRID) {
    expect_equal(
      table4_rows$qalys[table4_rows$non_cured_hazard_ratio == hr],
      s12$qalys[s12$non_cured_hazard_ratio == hr],
      info = paste("HR =", hr)
    )
  }
  # This is the finding R/09_scenarios.R's own module header documents for S12: QALYs fall (very
  # slightly) as the advantage strengthens, under the real Table-4-sourced Surgery row.
  expect_equal(unique(table4_rows$qaly_direction), "falls")
})

test_that("run_scenario_r2_surgery_sensitivity's table2 rows show a DIFFERENT sign than table4's -- the sensitivity check's actual finding, verified 2026-08-06", {
  r2 <- run_scenario_r2_surgery_sensitivity(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  table2_rows <- r2[r2$surgery_source == "table2_8wk_converted", ]
  table4_rows <- r2[r2$surgery_source == "table4_2wk_sourced", ]

  # The sign genuinely flips under Table 2's stickier Surgery dynamics: QALYs RISE as the
  # advantage strengthens (removing entries into a state that's now a SLOW lane, not a fast one),
  # the "intuitively expected" direction S12's own base case does NOT show. This is the concrete,
  # machine-checkable answer to R2's own question ("don't report the counterintuitive finding as
  # settled if the sign flips") -- it does flip.
  expect_equal(unique(table2_rows$qaly_direction), "rises")
  expect_false(identical(unique(table2_rows$qaly_direction), unique(table4_rows$qaly_direction)))

  # Cost still falls monotonically as the advantage strengthens under EITHER Surgery source --
  # Surgery remains an expensive state to enter regardless of how sticky it is once entered, so
  # this part of S12's own finding (NMB rises regardless of the QALY sign) should be robust.
  for (src in unique(r2$surgery_source)) {
    sub <- r2[r2$surgery_source == src, ]
    ordered <- sub[order(-sub$non_cured_hazard_ratio), ]
    expect_true(all(diff(ordered$total_cost) <= 1e-6), info = src)
  }
})

test_that("run_scenario_r2_surgery_sensitivity never mutates the real UST comparator matrix", {
  matrices_before <- build_all_transition_matrices(RAW_DIR)
  invisible(run_scenario_r2_surgery_sensitivity(raw_dir = RAW_DIR, proc_dir = PROC_DIR))
  matrices_after <- build_all_transition_matrices(RAW_DIR)
  expect_identical(matrices_before[["UST"]], matrices_after[["UST"]])
})
