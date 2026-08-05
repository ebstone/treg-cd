# Runs the structural scenarios S1-S12 (analysis_plan.md §10.3) that don't already have their own
# dedicated entry point. S3 (refractory population) is produced by
# analysis/run_full_analysis.R's own step 9/9 (run_refractory_scenario(), R/05) -- not repeated
# here. S6/S7/S9/S12 are NOT run here: they need new sourcing or new mechanics not yet built
# (R/09_scenarios.R's own module header explains each). This script covers S1, S2, S4, S5, S10,
# S11 -- the ones R/09_scenarios.R implements, first pass 2026-08-05.

source("R/utils/transition_matrix.R")
source("R/utils/life_table.R")
source("R/utils/refractory_multipliers.R")
source("R/00_derive_transition_probs.R")
source("R/01_decision_tree.R")
source("R/02_markov_engine.R")
source("R/03_cure_fraction_module.R")
source("R/04_costs_utilities.R")
source("R/05_deterministic_results.R")
source("R/09_scenarios.R")

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)

write_table <- function(df, name) {
  path <- file.path("output/tables", name)
  utils::write.csv(df, path, row.names = FALSE)
  cat("  wrote", path, sprintf("(%d rows)\n", nrow(df)))
}

# ---- S1: 2-year maintenance cap on vs. off -----------------------------------------------------
cat("=== S1: 2-year maintenance cap on vs. off ===\n")
s1 <- run_scenario_s1_cap()
print(s1, row.names = FALSE)
write_table(s1, "scenario_s1_cap.csv")

# ---- S2: ADA included vs. excluded -------------------------------------------------------------
cat("\n=== S2: ADA included vs. excluded (NMB table) ===\n")
s2_nmb <- run_scenario_s2_comparator_set()
print(s2_nmb, row.names = FALSE)
write_table(s2_nmb, "scenario_s2_comparator_set.csv")

cat("\n=== S2: required pi* at Treg's sourced price, lifetime horizon ===\n")
s2_headroom <- run_scenario_s2_headroom_at_sourced_price(n_cycles = HORIZON_CYCLES_LIFETIME,
                                                          baseline_age = ASSUMED_PATIENT_AGE_YEARS)
print(s2_headroom, row.names = FALSE)
write_table(s2_headroom, "scenario_s2_headroom_at_sourced_price.csv")

# ---- S4: cure-fraction anchor (named pi points x relapse-duration grid) ------------------------
cat("\n=== S4: cure-fraction anchor table (pi x relapse duration, at Treg's sourced price) ===\n")
s4 <- run_cure_fraction_anchor_table()
print(s4, row.names = FALSE)
write_table(s4, "scenario_s4_cure_fraction_anchor.csv")

# ---- S5: time horizon (6.15-year / 10-year / lifetime) -----------------------------------------
cat("\n=== S5: time horizon (6.15-year / 10-year / lifetime) ===\n")
s5 <- run_scenario_s5_horizon()
print(s5, row.names = FALSE)
write_table(s5, "scenario_s5_horizon.csv")

# ---- S10: discount rate (0% / 1.5% / 3% / 5%) --------------------------------------------------
cat("\n=== S10: discount rate (0% / 1.5% / 3% / 5%) ===\n")
s10 <- run_scenario_s10_discount_rate()
print(s10, row.names = FALSE)
write_table(s10, "scenario_s10_discount_rate.csv")

# ---- S11: SDR relapse re-enters Mild vs. Moderate-Severe Responder -----------------------------
cat("\n=== S11: SDR relapse destination (Mild vs. Moderate-Severe Responder) ===\n")
s11 <- run_scenario_s11_relapse_destination()
print(s11, row.names = FALSE)
write_table(s11, "scenario_s11_relapse_destination.csv")

# ---- S12: non-cured Treg = UST-equivalent vs. HR-advantaged --------------------------------------
cat("\n=== S12: non-cured Treg HR-advantage (illustrative grid, pi=0) ===\n")
s12 <- run_scenario_s12_non_cured_hr()
print(s12, row.names = FALSE)
write_table(s12, "scenario_s12_non_cured_hr.csv")

cat("\nDone. S6/S7/S9 NOT run here -- see R/09_scenarios.R's own module header for what each",
    "still needs (new sourcing or new mechanics, not built as of this pass).\n")
