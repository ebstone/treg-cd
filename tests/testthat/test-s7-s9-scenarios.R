# Tests for S7 (SDR utility = Remission vs. general-population) and S9 (healthcare-sector vs.
# societal perspective) -- both needed real sourcing (Jiang et al. 2021; Manceur et al. 2020),
# so tests live in their own file, matching the precedent set by S3
# (test-refractory-multipliers.R), not the no-new-sourcing scenarios in test-scenarios.R.

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

# ---- S7: population-utility data and lookup ----------------------------------------------------

test_that("load_population_utility_norms loads 14 rows (7 age bands x 2 sexes), all valid utilities", {
  norms <- load_population_utility_norms(RAW_DIR)
  expect_equal(nrow(norms), 14)
  expect_setequal(norms$sex, c("male", "female"))
  expect_true(all(norms$mean_utility >= -1 & norms$mean_utility <= 1))
  expect_equal(table(norms$age_band)[[1]], 2)  # exactly one male + one female row per band
})

test_that("general_population_utility_at_age handles exact integer band boundaries with no gaps", {
  norms <- load_population_utility_norms(RAW_DIR)
  # <25 (0.906, 0.931) vs 25-34 (0.907, 0.916): boundary at age 25
  expect_equal(general_population_utility_at_age(24, norms), mean(c(0.906, 0.931)))
  expect_equal(general_population_utility_at_age(24.9, norms), mean(c(0.906, 0.931)))
  expect_equal(general_population_utility_at_age(25, norms), mean(c(0.907, 0.916)))
  # 35-44 (0.841, 0.845): this project's cohort baseline age (Aliyev, 35)
  expect_equal(general_population_utility_at_age(35, norms), mean(c(0.841, 0.845)))
  expect_equal(general_population_utility_at_age(44.99, norms), mean(c(0.841, 0.845)))
  # Terminal open-ended band (75+): ages far beyond the table still resolve, not an error
  expect_equal(general_population_utility_at_age(75, norms), mean(c(0.786, 0.825)))
  expect_equal(general_population_utility_at_age(150, norms), mean(c(0.786, 0.825)))
})

test_that("general_population_utility_schedule(baseline_age = NULL) returns a length-1 vector at the reference age", {
  sched <- general_population_utility_schedule(10, baseline_age = NULL, raw_dir = RAW_DIR)
  expect_equal(length(sched), 1)
  expect_equal(sched, general_population_utility_at_age(ASSUMED_PATIENT_AGE_YEARS,
                                                          load_population_utility_norms(RAW_DIR)))
})

test_that("general_population_utility_schedule(baseline_age = X) is age-varying and the right length over a long horizon", {
  sched <- general_population_utility_schedule(HORIZON_CYCLES_LIFETIME, cycle_weeks = 2,
                                                baseline_age = 35, raw_dir = RAW_DIR)
  expect_equal(length(sched), HORIZON_CYCLES_LIFETIME + 1)
  expect_false(anyNA(sched))
  expect_true(length(unique(sched)) > 1)  # genuinely varies across a 65-year span
  norms <- load_population_utility_norms(RAW_DIR)
  expect_equal(sched[1], general_population_utility_at_age(35, norms))
})

# ---- S7: attach_sdr_costs_utilities()/attach_treg_costs_utilities() sdr_utility passthrough ----

test_that("attach_sdr_costs_utilities(sdr_utility = NULL) is byte-identical to omitting it, and uses Remission's own utility", {
  on_sdr <- c(0, 100, 95, 90)
  utilities <- c("Moderate-Severe" = 0.5, "Moderate-Severe Responder" = 0.6, "Mild" = 0.7,
                 "Remission" = 0.82, "Surgery" = 0.5, "Death" = 0)
  monitoring_costs <- c("Moderate-Severe" = 200, "Moderate-Severe Responder" = 200, "Mild" = 90,
                        "Remission" = 10, "Surgery" = 800, "Death" = 0)
  without_arg <- attach_sdr_costs_utilities(on_sdr, utilities, monitoring_costs)
  explicit_null <- attach_sdr_costs_utilities(on_sdr, utilities, monitoring_costs, sdr_utility = NULL)
  expect_equal(without_arg, explicit_null)

  custom <- attach_sdr_costs_utilities(on_sdr, utilities, monitoring_costs, sdr_utility = 0.9)
  expect_false(isTRUE(all.equal(custom$qalys_by_cycle, without_arg$qalys_by_cycle)))
  expect_equal(custom$qalys_by_cycle, without_arg$qalys_by_cycle / utilities[["Remission"]] * 0.9)
})

test_that("attach_treg_costs_utilities's sdr_utility = NULL default is byte-identical to omitting it", {
  # Build a minimal but real Treg arm result to exercise the full attach_treg_costs_utilities() path.
  matrices <- build_all_transition_matrices(RAW_DIR)
  induction <- load_published_induction(RAW_DIR)
  ust_row <- induction[induction$therapy == "UST", ]
  split <- run_decision_tree(ust_row)
  arm <- run_treg_arm(matrices[["UST"]], matrices[["CT"]], split$initial_on_biologic,
                       split$initial_on_ct, HORIZON_CYCLES_6YR, pi_sdr = 0.5, relapse_hazard_annual = 0.07)
  utilities <- load_health_state_utilities(PROC_DIR)
  monitoring_costs <- health_state_monitoring_costs()

  without_arg <- attach_treg_costs_utilities(arm, utilities, monitoring_costs)
  explicit_null <- attach_treg_costs_utilities(arm, utilities, monitoring_costs, sdr_utility = NULL)
  expect_equal(without_arg, explicit_null)
})

test_that("run_treg_arm_lifetime's sdr_utility_source = 'remission' default is byte-identical to omitting it", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  treg_price <- load_treg_dose_acquisition_cost(PROC_DIR)
  without_arg <- run_treg_arm_lifetime(HORIZON_CYCLES_6YR, pi_sdr = 0.5, relapse_hazard_annual = 0.07,
                                        price_usd = treg_price, matrices = matrices,
                                        raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  explicit_default <- run_treg_arm_lifetime(HORIZON_CYCLES_6YR, pi_sdr = 0.5, relapse_hazard_annual = 0.07,
                                             price_usd = treg_price, matrices = matrices,
                                             raw_dir = RAW_DIR, proc_dir = PROC_DIR,
                                             sdr_utility_source = "remission")
  expect_equal(without_arg, explicit_default)
})

test_that("run_scenario_s7_sdr_utility: identical at pi=0, general-population QALYs exceed remission-based QALYs at pi>0", {
  res <- run_scenario_s7_sdr_utility(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  expect_equal(nrow(res), length(PI_ANCHOR_GRID) * 2)

  pi0 <- res[res$pi_sdr == 0, ]
  expect_equal(pi0$qalys[pi0$sdr_utility_source == "remission"],
               pi0$qalys[pi0$sdr_utility_source == "general_population"])

  for (pi_val in setdiff(PI_ANCHOR_GRID, 0)) {
    rem <- res$qalys[res$pi_sdr == pi_val & res$sdr_utility_source == "remission"]
    gen <- res$qalys[res$pi_sdr == pi_val & res$sdr_utility_source == "general_population"]
    # General-population utility at age 35 (~0.843) exceeds Remission's utility (0.82), so the
    # general-population scenario gives strictly higher QALYs once any patient is actually cured.
    expect_true(gen > rem)
  }
})

# ---- S9: societal-perspective productivity cost --------------------------------------------------

test_that("productivity_cost_per_cycle scales the annual Manceur figure by cycle length", {
  expect_equal(productivity_cost_per_cycle(2), ANNUAL_INDIRECT_COST_USD * 2 / 52)
  expect_equal(productivity_cost_per_cycle(8), ANNUAL_INDIRECT_COST_USD * 8 / 52)
})

test_that("societal_monitoring_costs adds the productivity cost to every living state and leaves Death at 0", {
  base <- health_state_monitoring_costs()
  societal <- societal_monitoring_costs()
  living <- setdiff(names(base), "Death")
  expect_equal(societal[living], base[living] + productivity_cost_per_cycle())
  expect_equal(unname(societal["Death"]), 0)
})

test_that("run_comparator_arm_lifetime's perspective = 'healthcare_sector' default is byte-identical to omitting it", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  for (tx in COMPARATOR_THERAPIES) {
    without_arg <- run_comparator_arm_lifetime(tx, HORIZON_CYCLES_6YR, matrices, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
    explicit_default <- run_comparator_arm_lifetime(tx, HORIZON_CYCLES_6YR, matrices, raw_dir = RAW_DIR,
                                                     proc_dir = PROC_DIR, perspective = "healthcare_sector")
    expect_equal(without_arg, explicit_default)
  }
})

test_that("run_comparator_arm_lifetime(perspective = 'societal') raises cost but leaves QALYs unchanged", {
  matrices <- build_all_transition_matrices(RAW_DIR)
  for (tx in COMPARATOR_THERAPIES) {
    healthcare <- run_comparator_arm_lifetime(tx, HORIZON_CYCLES_6YR, matrices, raw_dir = RAW_DIR, proc_dir = PROC_DIR)
    societal <- run_comparator_arm_lifetime(tx, HORIZON_CYCLES_6YR, matrices, raw_dir = RAW_DIR,
                                             proc_dir = PROC_DIR, perspective = "societal")
    expect_equal(healthcare$qalys, societal$qalys)
    expect_true(societal$total_cost > healthcare$total_cost)
  }
})

test_that("run_scenario_s9_perspective's healthcare_sector rows reproduce run_base_case(), and societal costs exceed healthcare_sector costs for every arm", {
  scenario <- run_scenario_s9_perspective(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  base <- run_base_case(raw_dir = RAW_DIR, proc_dir = PROC_DIR)
  hc_rows <- scenario[scenario$perspective == "healthcare_sector", setdiff(names(scenario), "perspective")]
  rownames(hc_rows) <- NULL; rownames(base) <- NULL
  expect_equal(hc_rows, base)

  for (arm in c(COMPARATOR_THERAPIES, "TREG")) {
    hc_cost <- scenario$total_cost[scenario$perspective == "healthcare_sector" & scenario$intervention == arm]
    soc_cost <- scenario$total_cost[scenario$perspective == "societal" & scenario$intervention == arm]
    hc_qalys <- scenario$qalys[scenario$perspective == "healthcare_sector" & scenario$intervention == arm]
    soc_qalys <- scenario$qalys[scenario$perspective == "societal" & scenario$intervention == arm]
    expect_true(soc_cost > hc_cost)
    expect_equal(hc_qalys, soc_qalys)
  }
})
