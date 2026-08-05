# Runs the base-case deterministic analysis end-to-end (source()s the R/ pipeline in order).
#
# FIRST PASS, 2026-08-05, alongside R/05_deterministic_results.R's first pass. Produces
# output/tables/base_case_results.csv: UST/IFX/ADA at their real base case, plus TREG at pi=0 (the
# Null floor scenario -- see R/05_deterministic_results.R's module header for why there is no
# other numeric "Treg base case" to report yet, per Decision 4, analysis_plan.md sec 15). Does NOT
# yet produce the headroom frontier or any figures -- those need a plotting decision (ggplot2 vs.
# base graphics, not yet made) and are left for a follow-up pass; headroom_frontier() itself is
# already usable interactively from R/05_deterministic_results.R in the meantime.

source("R/utils/transition_matrix.R")
source("R/00_derive_transition_probs.R")
source("R/01_decision_tree.R")
source("R/02_markov_engine.R")
source("R/03_cure_fraction_module.R")
source("R/04_costs_utilities.R")
source("R/05_deterministic_results.R")

base_case <- run_base_case()

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(base_case, "output/tables/base_case_results.csv", row.names = FALSE)

cat("Base case (6.15-year horizon; TREG row is the pi=0 Null scenario, not a Treg efficacy claim):\n\n")
print(base_case, row.names = FALSE)
cat("\nWrote output/tables/base_case_results.csv\n")
