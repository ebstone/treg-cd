# Single entry point: reproduces every number, figure and table in the manuscript from a clean
# checkout (analysis_plan.md §12.2 convention). If a manuscript number can't be produced by this
# script, it does not go in the paper.
#
# FIRST FULL PASS, 2026-08-05, run on commit 77bcda3 (post-A16-fix, post-biosimilar-re-pricing --
# see docs/model_audit_v6.md A15/A16 and README.md's Status section before citing any number this
# script produces; everything from before that commit is superseded). Wires together every module
# built so far (R/00-R/08) into one reproducible run:
#   1. Base case (R/05 run_base_case()) -- UST/IFX/ADA at their real base case, TREG at the pi=0
#      Null floor (Decision 4: no other numeric Treg base case exists yet).
#   2. Headroom frontier (Aim 4, R/05 headroom_frontier()) -- (pi, price) indifference curve,
#      swept over price at each WTP threshold analysis_plan.md sec 7.1 asks for ($50k/$100k/
#      $150k, primary $150k).
#   3. 10,000-draw PSA (R/06 run_psa()) -- analysis_plan.md sec 9.3's own stated minimum, plus the
#      cost-effectiveness plane (TREG vs. each comparator) and CEAC.
#   4. EVPI surface (R/07 evpi_surface(), sec 9.2) over the same price x WTP grid, and EVPPI by
#      parameter subset (R/07 evppi_by_subset(), sec 9.3) at the primary WTP, plus a convergence
#      check on subset A (sec 9.3's explicit requirement) -- NOT the voi/BCEA cross-check
#      (cross_check_voi(), inert in this sandbox -- see R/07's own header and
#      docs/treg-cd_decision_resolutions_2026-08-05.md item 8 for why that's deprioritised, not
#      skipped by oversight). population_evpi() is also not called here: effective_population has
#      no default because the fractions Decision 6 needs to compute one aren't sourced yet
#      (R/07's own header) -- inventing one here would be exactly the fabrication this project's
#      methodology exists to avoid.
#   5. Deterministic + probabilistic EJP (Aim 1, R/08) -- ejp_deterministic() across a pi grid at
#      each WTP threshold, ejp_frontier() as an independent cross-check of headroom_frontier()
#      from the opposite direction (module header), ejp_probabilistic()/ejp_p50() from the PSA
#      draws already in hand, and gross_margin_over_cogs() alongside every P* this script produces.
#
# Explicitly OUT of scope for this pass, same as prior module-level passes flagged in their own
# headers: any figure (README: ggplot2-vs-base-graphics plotting decision not yet made -- tables
# only here), the structural scenarios S1-S12 (analysis/run_scenario_analyses.R, still a stub),
# and everything in the handoff doc's "still open" list (h sweep, pi prior-sensitivity reweighting,
# lifetime horizon, half-cycle correction, refractory multipliers) -- this script runs the pipeline
# as it exists today, it does not extend it.

source("R/utils/transition_matrix.R")
source("R/00_derive_transition_probs.R")
source("R/01_decision_tree.R")
source("R/02_markov_engine.R")
source("R/03_cure_fraction_module.R")
source("R/04_costs_utilities.R")
source("R/05_deterministic_results.R")
source("R/06_psa.R")
source("R/07_evpi_evppi.R")
source("R/08_ejp.R")

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)

WTP_PRIMARY <- 150000
WTP_GRID <- c(50000, 100000, 150000)       # analysis_plan.md sec 7.1: report all three, $150k primary
PRICE_GRID <- seq(0, 100000, by = 5000)    # spans COGS ($4,979) to well above acquisition ($19,917)
PI_GRID <- seq(0, 1, by = 0.1)

write_table <- function(df, name) {
  path <- file.path("output/tables", name)
  utils::write.csv(df, path, row.names = FALSE)
  cat("  wrote", path, sprintf("(%d rows)\n", nrow(df)))
}

# ---- 1. Base case ---------------------------------------------------------------------------
cat("=== 1/5: Base case ===\n")
base_case <- run_base_case()
print(base_case, row.names = FALSE)
write_table(base_case, "base_case_results.csv")

# ---- 2. Headroom frontier (Aim 4) ------------------------------------------------------------
cat("\n=== 2/5: Headroom frontier (Aim 4) ===\n")
frontier <- do.call(rbind, lapply(WTP_GRID, function(w) headroom_frontier(PRICE_GRID, w)))
write_table(frontier, "headroom_frontier.csv")

# ---- 3. PSA: draws, CE plane, CEAC ------------------------------------------------------------
cat("\n=== 3/5: PSA (10,000 draws; this is the slow step, ~3-4 min) ===\n")
psa_results <- run_psa(n_draws = 10000)
write_table(psa_results, "psa_draws.csv")

comparators <- setdiff(unique(psa_results$intervention), "TREG")
ce_plane <- do.call(rbind, lapply(comparators, function(comp) {
  d <- psa_cost_effectiveness_plane(psa_results, "TREG", comp)
  d$comparator <- comp
  d
}))
write_table(ce_plane, "psa_ce_plane.csv")

ceac_grid <- seq(0, 300000, by = 10000)
ceac <- psa_ceac(psa_results, ceac_grid)
write_table(ceac, "psa_ceac.csv")

# ---- 4. EVPI / EVPPI --------------------------------------------------------------------------
cat("\n=== 4/5: EVPI surface + EVPPI by subset ===\n")
evpi_wtp_grid <- seq(0, 300000, by = 25000)
evpi_surf <- evpi_surface(psa_results, PRICE_GRID, evpi_wtp_grid)
write_table(evpi_surf, "evpi_surface.csv")

evppi_subset <- evppi_by_subset(psa_results, WTP_PRIMARY)
print(evppi_subset, row.names = FALSE)
write_table(evppi_subset, "evppi_by_subset.csv")

convergence_counts <- c(500, 1000, 2000, 5000, 10000)
evppi_conv <- evppi_convergence(psa_results, "pi_sdr", WTP_PRIMARY, convergence_counts)
write_table(evppi_conv, "evppi_convergence_subset_A.csv")

# ---- 5. EJP (deterministic + probabilistic) + gross margin over COGS --------------------------
cat("\n=== 5/5: EJP ===\n")
ejp_det <- do.call(rbind, lapply(WTP_GRID, function(w) {
  do.call(rbind, lapply(PI_GRID, function(pi) {
    res <- ejp_deterministic(pi, w)
    data.frame(pi_sdr = pi, wtp_usd = w, p_star = res$p_star, comparator = res$comparator,
               gross_margin_over_cogs = gross_margin_over_cogs(res$p_star))
  }))
}))
write_table(ejp_det, "ejp_deterministic.csv")

ejp_front <- ejp_frontier(seq(0, 1, by = 0.05), WTP_PRIMARY)
write_table(ejp_front, "ejp_frontier.csv")

ejp_prob <- do.call(rbind, lapply(WTP_GRID, function(w) {
  res <- ejp_probabilistic(psa_results, w)
  data.frame(wtp_usd = w, median = res$median, ci_low = res$ci_low, ci_high = res$ci_high,
             p50 = ejp_p50(res$draws$p_star),
             gross_margin_over_cogs_at_median = gross_margin_over_cogs(res$median))
}))
print(ejp_prob, row.names = FALSE)
write_table(ejp_prob, "ejp_probabilistic.csv")

cat("\nDone. All tables in output/tables/.\n")
