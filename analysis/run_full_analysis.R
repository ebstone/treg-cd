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
#   6. Lifetime-horizon base case + headroom frontier (added 2026-08-05, same day as the
#      lifetime-horizon build itself -- R/05's own module header, R/utils/life_table.R) --
#      deterministic only in this pass (PSA/EVPI/EVPPI/probabilistic EJP still run at
#      HORIZON_CYCLES_6YR; see R/05's module header for why extending age_adjust_matrix() to
#      ~40,000 PSA arm-runner calls is a real, unbenchmarked performance question, not done here).
#   7. h sweep (added 2026-08-05, docs/treg-cd_decision_resolutions_2026-08-05.md sec 3.3): the
#      (pi, T, price) headroom surface at the LIFETIME horizon (R/05 headroom_frontier_by_duration()
#      -- "h interacts with the horizon... must be run jointly"), plus verify_pi_factorization()'s
#      numeric check that Treg's incremental QALY is exactly pi * g(h) (found to hold to
#      floating-point precision, not merely approximately -- R/05's own docstring on that
#      function).
#   8. pi prior-sensitivity on EVPPI (added 2026-08-05, same memo sec 3.2): evppi_by_subset()
#      re-estimated under Beta(1,3) and Beta(2,2) alternative priors for pi via importance
#      reweighting of the SAME PSA draws already produced in step 3 (R/07 evppi_prior_sensitivity())
#      -- no new simulation. Reports the ranking-is-stable-but-level-is-prior-dependent finding the
#      memo asks for.
#
# Explicitly OUT of scope for this pass, same as prior module-level passes flagged in their own
# headers: any figure (README: ggplot2-vs-base-graphics plotting decision not yet made -- tables
# only here), the structural scenarios S1-S12 (analysis/run_scenario_analyses.R, still a stub),
# half-cycle correction (A12), and refractory-population multipliers -- still open per the
# resolutions memo's own priority order (README.md's Status section).

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
cat("=== 1/8: Base case ===\n")
base_case <- run_base_case()
print(base_case, row.names = FALSE)
write_table(base_case, "base_case_results.csv")

# ---- 2. Headroom frontier (Aim 4) ------------------------------------------------------------
cat("\n=== 2/8: Headroom frontier (Aim 4) ===\n")
frontier <- do.call(rbind, lapply(WTP_GRID, function(w) headroom_frontier(PRICE_GRID, w)))
write_table(frontier, "headroom_frontier.csv")

# ---- 3. PSA: draws, CE plane, CEAC ------------------------------------------------------------
cat("\n=== 3/8: PSA (10,000 draws; this is the slow step, ~3-4 min) ===\n")
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
cat("\n=== 4/8: EVPI surface + EVPPI by subset ===\n")
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
cat("\n=== 5/8: EJP ===\n")
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

# ---- 6. Lifetime horizon (deterministic only this pass -- see module header) -------------------
cat("\n=== 6/8: Lifetime horizon (deterministic base case + headroom frontier) ===\n")
base_case_lifetime <- run_base_case(n_cycles = HORIZON_CYCLES_LIFETIME, baseline_age = ASSUMED_PATIENT_AGE_YEARS)
print(base_case_lifetime, row.names = FALSE)
write_table(base_case_lifetime, "base_case_results_lifetime.csv")

frontier_lifetime <- do.call(rbind, lapply(WTP_GRID, function(w) {
  headroom_frontier(PRICE_GRID, w, n_cycles = HORIZON_CYCLES_LIFETIME, baseline_age = ASSUMED_PATIENT_AGE_YEARS)
}))
write_table(frontier_lifetime, "headroom_frontier_lifetime.csv")

# The headline check this whole addition exists to make: at Treg's actual sourced acquisition
# price, what durable cure fraction does the lifetime horizon require, at each WTP threshold?
treg_price <- load_treg_dose_acquisition_cost()
headroom_at_sourced_price <- do.call(rbind, lapply(WTP_GRID, function(w) {
  res <- headroom_pi_star(treg_price, wtp_usd = w, n_cycles = HORIZON_CYCLES_LIFETIME,
                           baseline_age = ASSUMED_PATIENT_AGE_YEARS)
  data.frame(wtp_usd = w, price_usd = treg_price, pi_star = res$pi_star, feasible = res$feasible)
}))
cat("\nRequired durable cure fraction (pi*) at Treg's sourced acquisition price ($", format(treg_price, big.mark = ","), "), lifetime horizon:\n", sep = "")
print(headroom_at_sourced_price, row.names = FALSE)
write_table(headroom_at_sourced_price, "headroom_at_sourced_price_lifetime.csv")

# ---- 7. h sweep: (pi, T, price) headroom surface + the pi*g(h) factorisation check ------------
cat("\n=== 7/8: h sweep (relapse duration T) at the lifetime horizon ===\n")
# A coarser price grid than PRICE_GRID (11 points, not 21): this section multiplies price x WTP x
# duration (11 x 3 x 5 = 165 headroom_pi_star() evaluations, each a full lifetime-horizon Markov
# trace) rather than PRICE_GRID's own price x WTP (63) -- kept coarser so this section's added
# runtime stays a fraction of, not a multiple of, the rest of this script's.
PRICE_GRID_COARSE <- seq(0, 100000, by = 10000)
duration_surface <- do.call(rbind, lapply(WTP_GRID, function(w) {
  headroom_frontier_by_duration(PRICE_GRID_COARSE, w, n_cycles = HORIZON_CYCLES_LIFETIME,
                                 baseline_age = ASSUMED_PATIENT_AGE_YEARS)
}))
write_table(duration_surface, "headroom_frontier_by_duration_lifetime.csv")

# The factorisation check itself, at Treg's own sourced price and the base-case T = 10yr duration
# -- documents, reproducibly, that pi * g(h) is exact to floating-point precision at this horizon
# (R/05's verify_pi_factorization() docstring has the full argument for why).
pi_factorization_check <- verify_pi_factorization(
  price_usd = treg_price, relapse_hazard_annual = duration_to_hazard(DEFAULT_RELAPSE_DURATION_YEARS),
  n_cycles = HORIZON_CYCLES_LIFETIME, baseline_age = ASSUMED_PATIENT_AGE_YEARS
)
cat("pi * g(h) factorisation check (max relative error across the pi grid): ",
    format(max(pi_factorization_check$rel_error), scientific = TRUE), "\n", sep = "")
write_table(pi_factorization_check, "pi_factorization_check.csv")

# ---- 8. pi prior-sensitivity on EVPPI --------------------------------------------------------
cat("\n=== 8/8: pi prior-sensitivity on EVPPI ===\n")
# Re-estimates evppi_by_subset() (already computed once, unweighted, in step 4) under Beta(1,3)
# and Beta(2,2) alternative priors for pi via importance reweighting of the SAME psa_results
# drawn in step 3 -- no new simulation (R/07's own module header, sec 3.2).
evppi_priors <- evppi_prior_sensitivity(psa_results, WTP_PRIMARY)
print(evppi_priors[, c("prior", "subset", "evppi", "rank_within_prior")], row.names = FALSE)
write_table(evppi_priors, "evppi_prior_sensitivity.csv")

cat("\nDone. All tables in output/tables/.\n")
