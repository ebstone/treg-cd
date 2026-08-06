# Single entry point: reproduces every number, figure and table in the manuscript from a clean
# checkout (analysis_plan.md §12.2 convention). If a manuscript number can't be produced by this
# script, it does not go in the paper.
#
# FIRST FULL PASS, 2026-08-05, run on commit 77bcda3 (post-A16-fix, post-biosimilar-re-pricing --
# see docs/model_audit_v6.md A15/A16 and README.md's Status section before citing any number this
# script produces; everything from before that commit is superseded).
#
# RESTRUCTURED 2026-08-06 (peer review 2026-08-05, B1/B2, docs/model_audit_v6.md A17,
# README.md's Status section). Two changes, both superseding every number this script produced
# before this commit:
#
# 1. **B1 (induction cycles) is fixed upstream** (R/00_derive_transition_probs.R) -- nothing in
#    THIS script changed to pick that up; every arm-runner call below inherits the corrected
#    induction split automatically.
# 2. **The lifetime horizon, not the 6.15-year comparator, is now what this script treats as
#    primary.** analysis_plan.md §4.3 adopts the lifetime horizon as the base case; the review
#    found the PSA, EVPI surface, EVPPI, and deterministic/probabilistic EJP were still running
#    at the 6.15-year horizon by default -- the horizon §4.3 itself says is structurally
#    incapable of representing a cure. Every primary output below (base case, headroom frontier,
#    PSA, EVPI/EVPPI, EJP) now runs at HORIZON_CYCLES_LIFETIME with age/sex-specific mortality
#    (ASSUMED_PATIENT_AGE_YEARS) FIRST, written to the unsuffixed table names. The 6.15-year and
#    10-year horizons are retained, exactly as analysis_plan.md §10.3 already frames them, as the
#    named S5 comparability scenario -- run_scenario_s5_horizon() (R/09_scenarios.R) reports the
#    deterministic base case and headroom frontier at all three horizons side by side (with an
#    explicit `horizon` column, closing the review's own B2 acceptance test), written to
#    `..._by_horizon.csv`. The full 10,000-draw PSA/EVPI/EVPPI/probabilistic-EJP pipeline is
#    reported ONLY at the primary (lifetime) horizon, not tripled across all three horizons --
#    the original script never ran that pipeline at more than one horizon either (it ran only at
#    6.15yr, deterministic-only at lifetime); this pass flips which single horizon that
#    10,000-draw pipeline targets, it does not multiply the pipeline's own cost by 3x.
#
# Performance: a lifetime-horizon PSA (1,691 cycles, with age/sex-specific mortality) is
# substantially more expensive than the previous 6.15-year default (160 cycles, no mortality
# adjustment). R/06_psa.R's run_psa() now hoists the three comparator arms' occupancy-trace
# simulation out of the per-draw loop (simulate_comparator_arm_lifetime(), R/05 -- see run_psa()'s
# own "Performance note") since R/06 doesn't sample any transition probability, which is most of
# what makes this tractable at all; TREG still simulates once per draw (pi genuinely varies).
# **No R toolchain was available in the environment this restructuring was done in, so this
# script has not actually been re-run end to end since the change** -- the next session with R
# available should run it, watch actual wall-clock time on the PSA step, and reduce
# `PSA_N_DRAWS` (with that reduction stated explicitly, not silently) if a full 10,000-draw
# lifetime PSA proves impractical, per the review's own stated fallback.
#
# Everything else about what this script computes (Aims 1-4, the h sweep, prior-sensitivity,
# refractory scenario) is unchanged from the first pass -- see the per-step comments below.
#
# RELAPSE HAZARD, 2026-08-06 (peer review's B3, docs/model_audit_v6.md A19). **Every Treg-involving
# output below is superseded again by this change, for the third time in this file's short history
# (A16, then B1, now B3).** Nothing in THIS script changed to pick it up, exactly as with B1: the
# post-cure relapse hazard h now defaults to duration_to_hazard(DEFAULT_RELAPSE_DURATION_YEARS)
# (a median 10 years of drug-free remission) in headroom_pi_star()/headroom_frontier() (R/05) and
# ejp_deterministic()/ejp_frontier() (R/08) instead of 0 (a permanent cure), and R/06_psa.R samples
# it per draw around that same mean, so every call below inherits the base-case durability
# assumption automatically. Two things follow, both worth stating before anyone diffs the tables:
#
# - **What moves:** headroom_frontier.csv, headroom_at_sourced_price.csv, ejp_deterministic.csv,
#   ejp_frontier.csv, psa_draws.csv (plus a new relapse_hazard_annual column), psa_ce_plane.csv,
#   psa_ceac.csv, evpi_surface.csv, evppi_by_subset.csv (plus new B / A u B rows),
#   evppi_prior_sensitivity.csv, ejp_probabilistic.csv, headroom_frontier_by_horizon.csv. pi* rises
#   and P* falls, both substantially -- a cure that half-decays in a decade is worth much less over
#   a lifetime horizon than one that never decays.
# - **What does NOT move:** base_case_results.csv, base_case_results_by_horizon.csv and
#   refractory_scenario_results.csv are bit-identical, because they report Treg only at pi = 0,
#   where no patient ever enters the SDR state and h multiplies nothing (run_base_case()'s own
#   comment in R/05 has the argument, and a test asserts it). Every comparator (UST/IFX/ADA) figure
#   everywhere in this script is likewise untouched: h reaches the model only through Treg's SDR
#   track. headroom_frontier_by_duration.csv is also unchanged -- it always passed T explicitly.

source("R/utils/transition_matrix.R")
source("R/utils/life_table.R")
source("R/00_derive_transition_probs.R")
source("R/01_decision_tree.R")
source("R/02_markov_engine.R")
source("R/03_cure_fraction_module.R")
source("R/04_costs_utilities.R")
source("R/05_deterministic_results.R")
source("R/06_psa.R")
source("R/07_evpi_evppi.R")
source("R/08_ejp.R")
source("R/09_scenarios.R")

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)

WTP_PRIMARY <- 150000
WTP_GRID <- c(50000, 100000, 150000)       # analysis_plan.md sec 7.1: report all three, $150k primary
PRICE_GRID <- seq(0, 100000, by = 5000)    # spans COGS ($4,979) to well above acquisition ($19,917)
PI_GRID <- seq(0, 1, by = 0.1)
PSA_N_DRAWS <- 10000                       # analysis_plan.md sec 9.3's own stated minimum

# Primary-horizon PSA arguments, threaded through every step that needs them (PSA itself, and
# every downstream EVPI/EVPPI/EJP-probabilistic step that consumes its draws).
PRIMARY_HORIZON_ARGS <- list(n_cycles = HORIZON_CYCLES_LIFETIME, baseline_age = ASSUMED_PATIENT_AGE_YEARS)

write_table <- function(df, name) {
  path <- file.path("output/tables", name)
  utils::write.csv(df, path, row.names = FALSE)
  cat("  wrote", path, sprintf("(%d rows)\n", nrow(df)))
}

# =============================================================================================
# PRIMARY OUTPUTS -- lifetime horizon (analysis_plan.md §4.3's adopted base case)
# =============================================================================================

# ---- 1. Base case -----------------------------------------------------------------------------
cat("=== 1/10: Base case (lifetime horizon) ===\n")
base_case <- run_base_case(n_cycles = HORIZON_CYCLES_LIFETIME, baseline_age = ASSUMED_PATIENT_AGE_YEARS)
print(base_case, row.names = FALSE)
write_table(base_case, "base_case_results.csv")

# ---- 2. Headroom frontier (Aim 4) --------------------------------------------------------------
cat("\n=== 2/10: Headroom frontier (Aim 4, lifetime horizon) ===\n")
frontier <- do.call(rbind, lapply(WTP_GRID, function(w) {
  headroom_frontier(PRICE_GRID, w, n_cycles = HORIZON_CYCLES_LIFETIME, baseline_age = ASSUMED_PATIENT_AGE_YEARS)
}))
write_table(frontier, "headroom_frontier.csv")

# The headline check the lifetime horizon exists to make: at Treg's actual sourced acquisition
# price, what durable cure fraction does the lifetime horizon require, at each WTP threshold?
# `relapse_hazard_annual` is left implicit here (and in the frontier above, and in step 5's EJP
# calls) so it picks up R/05's own base-case default -- the 10-year median SDR duration, since
# 2026-08-06. That is deliberate, not an oversight: the base case belongs in one place, not
# restated at every call site where it could drift. The h SWEEP in step 6 is where h is passed
# explicitly, precisely because that step's job is to vary it.
treg_price <- load_treg_dose_acquisition_cost()
headroom_at_sourced_price <- do.call(rbind, lapply(WTP_GRID, function(w) {
  res <- headroom_pi_star(treg_price, wtp_usd = w, n_cycles = HORIZON_CYCLES_LIFETIME,
                           baseline_age = ASSUMED_PATIENT_AGE_YEARS)
  data.frame(wtp_usd = w, price_usd = treg_price, pi_star = res$pi_star, feasible = res$feasible)
}))
cat("\nRequired durable cure fraction (pi*) at Treg's sourced acquisition price ($", format(treg_price, big.mark = ","), "), lifetime horizon:\n", sep = "")
print(headroom_at_sourced_price, row.names = FALSE)
write_table(headroom_at_sourced_price, "headroom_at_sourced_price.csv")

# ---- 3. PSA: draws, CE plane, CEAC (lifetime horizon) ------------------------------------------
cat("\n=== 3/10: PSA (", PSA_N_DRAWS, " draws, lifetime horizon -- the slow step; see this",
    " script's own header on expected runtime) ===\n", sep = "")
psa_results <- do.call(run_psa, c(list(n_draws = PSA_N_DRAWS), PRIMARY_HORIZON_ARGS))
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

# ---- 4. EVPI / EVPPI (lifetime horizon) ---------------------------------------------------------
cat("\n=== 4/10: EVPI surface + EVPPI by subset (lifetime horizon) ===\n")
evpi_wtp_grid <- seq(0, 300000, by = 25000)
evpi_surf <- evpi_surface(psa_results, PRICE_GRID, evpi_wtp_grid)
write_table(evpi_surf, "evpi_surface.csv")

evppi_subset <- evppi_by_subset(psa_results, WTP_PRIMARY)
print(evppi_subset, row.names = FALSE)
write_table(evppi_subset, "evppi_by_subset.csv")

convergence_counts <- unique(pmin(c(500, 1000, 2000, 5000, PSA_N_DRAWS), PSA_N_DRAWS))
evppi_conv <- evppi_convergence(psa_results, "pi_sdr", WTP_PRIMARY, convergence_counts)
write_table(evppi_conv, "evppi_convergence_subset_A.csv")

# Subset B gets its own convergence check (added 2026-08-06 alongside B3), not just subset A's:
# analysis_plan.md §9.3 asks for a convergence check because regression-based EVPPI estimates are
# unstable at small draw counts, and that instability is worst for a SMALL EVPPI -- which subset B
# is expected to be relative to A. A converged A tells you nothing about whether B's number is
# real, and Aim 3's comparison rests on B being a real number, so it is checked directly.
evppi_conv_b <- evppi_convergence(psa_results, "relapse_hazard_annual", WTP_PRIMARY, convergence_counts)
write_table(evppi_conv_b, "evppi_convergence_subset_B.csv")

# ---- 5. EJP (deterministic + probabilistic) + gross margin over COGS (lifetime horizon) --------
cat("\n=== 5/10: EJP (lifetime horizon) ===\n")
ejp_det <- do.call(rbind, lapply(WTP_GRID, function(w) {
  do.call(rbind, lapply(PI_GRID, function(pi) {
    res <- ejp_deterministic(pi, w, n_cycles = HORIZON_CYCLES_LIFETIME, baseline_age = ASSUMED_PATIENT_AGE_YEARS)
    data.frame(pi_sdr = pi, wtp_usd = w, p_star = res$p_star, comparator = res$comparator,
               gross_margin_over_cogs = gross_margin_over_cogs(res$p_star))
  }))
}))
write_table(ejp_det, "ejp_deterministic.csv")

ejp_front <- ejp_frontier(seq(0, 1, by = 0.05), WTP_PRIMARY, n_cycles = HORIZON_CYCLES_LIFETIME,
                           baseline_age = ASSUMED_PATIENT_AGE_YEARS)
write_table(ejp_front, "ejp_frontier.csv")

# ejp_deterministic(pi, wtp)/headroom_pi_star(P*, wtp) round-trip against each other regardless of
# horizon (test-ejp.R) -- this is Aim 1/Aim 4's own internal cross-check, unaffected by which
# horizon is primary.
ejp_prob <- do.call(rbind, lapply(WTP_GRID, function(w) {
  res <- ejp_probabilistic(psa_results, w)
  data.frame(wtp_usd = w, median = res$median, ci_low = res$ci_low, ci_high = res$ci_high,
             p50 = ejp_p50(res$draws$p_star),
             gross_margin_over_cogs_at_median = gross_margin_over_cogs(res$median))
}))
print(ejp_prob, row.names = FALSE)
write_table(ejp_prob, "ejp_probabilistic.csv")

# ---- 6. h sweep: (pi, T, price) headroom surface + the pi*g(h) factorisation check -------------
cat("\n=== 6/10: h sweep (relapse duration T, lifetime horizon) ===\n")
# A coarser price grid than PRICE_GRID (11 points, not 21): this section multiplies price x WTP x
# duration (11 x 3 x 5 = 165 headroom_pi_star() evaluations, each a full lifetime-horizon Markov
# trace) rather than PRICE_GRID's own price x WTP (63) -- kept coarser so this section's added
# runtime stays a fraction of, not a multiple of, the rest of this script's.
PRICE_GRID_COARSE <- seq(0, 100000, by = 10000)
duration_surface <- do.call(rbind, lapply(WTP_GRID, function(w) {
  headroom_frontier_by_duration(PRICE_GRID_COARSE, w, n_cycles = HORIZON_CYCLES_LIFETIME,
                                 baseline_age = ASSUMED_PATIENT_AGE_YEARS)
}))
write_table(duration_surface, "headroom_frontier_by_duration.csv")

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

# ---- 7. pi prior-sensitivity on EVPPI (lifetime horizon) ---------------------------------------
cat("\n=== 7/10: pi prior-sensitivity on EVPPI (lifetime horizon) ===\n")
# Re-estimates evppi_by_subset() (already computed once, unweighted, in step 4) under Beta(1,3)
# and Beta(2,2) alternative priors for pi via importance reweighting of the SAME psa_results
# drawn in step 3 -- no new simulation (R/07's own module header, sec 3.2).
evppi_priors <- evppi_prior_sensitivity(psa_results, WTP_PRIMARY)
print(evppi_priors[, c("prior", "subset", "evppi", "rank_within_prior")], row.names = FALSE)
write_table(evppi_priors, "evppi_prior_sensitivity.csv")

# =============================================================================================
# COMPARABILITY SCENARIOS -- S5 (horizon), S3 (refractory population)
# =============================================================================================
# analysis_plan.md §10.3: the 6.15-year and 10-year horizons are retained as NAMED comparability
# scenarios against Aliyev's own original design and the paper's earlier drafts, not as
# alternative "base cases" competing with the lifetime horizon above. Deterministic only -- the
# full PSA/EVPI/EVPPI/probabilistic-EJP pipeline is reported once, at the primary horizon, above
# (this script's own header explains why tripling it wasn't done).

# ---- 8. S5: base case + headroom frontier at all three horizons, side by side ------------------
cat("\n=== 8/10: S5 horizon comparability (6.15yr / 10yr / lifetime, deterministic only) ===\n")
base_case_by_horizon <- run_scenario_s5_horizon()
print(base_case_by_horizon, row.names = FALSE)
write_table(base_case_by_horizon, "base_case_results_by_horizon.csv")

frontier_by_horizon <- do.call(rbind, lapply(WTP_GRID, function(w) {
  rbind(
    cbind(horizon = "6.15_year", headroom_frontier(PRICE_GRID, w, n_cycles = HORIZON_CYCLES_6YR)),
    cbind(horizon = "10_year", headroom_frontier(PRICE_GRID, w, n_cycles = HORIZON_CYCLES_10YR)),
    cbind(horizon = "lifetime", headroom_frontier(PRICE_GRID, w, n_cycles = HORIZON_CYCLES_LIFETIME,
                                                   baseline_age = ASSUMED_PATIENT_AGE_YEARS))
  )
}))
write_table(frontier_by_horizon, "headroom_frontier_by_horizon.csv")

# ---- 9. Refractory-population scenario (S3) -----------------------------------------------------
cat("\n=== 9/10: Refractory-population scenario (S3, 6.15-year comparability horizon) ===\n")
# run_refractory_scenario() (R/05, module header) -- the biologic_naive rows reproduce step 1's
# base_case exactly at the SAME horizon it's run at; the refractory rows apply the
# UNITI-1-vs-UNITI-2-sourced multipliers (R/utils/refractory_multipliers.R) to UST/IFX/ADA's
# induction and maintenance matrices only (CT and TREG are unaffected -- see that file's own
# header for why). Run at the 6.15-year horizon, matching analysis_plan.md's original refractory
# sourcing (UNITI-1/UNITI-2 are themselves 6-8 week induction / 44-week maintenance trials, not a
# lifetime-horizon design) -- not re-derived at the lifetime horizon in this pass.
refractory_scenario <- run_refractory_scenario()
print(refractory_scenario, row.names = FALSE)
write_table(refractory_scenario, "refractory_scenario_results.csv")

cat("\n=== 10/10: done ===\n")
cat("Primary (lifetime-horizon) outputs: base_case_results.csv, headroom_frontier.csv,",
    "headroom_at_sourced_price.csv, psa_draws.csv, psa_ce_plane.csv, psa_ceac.csv,",
    "evpi_surface.csv, evppi_by_subset.csv, evppi_convergence_subset_A.csv,",
    "evppi_convergence_subset_B.csv, ejp_deterministic.csv,",
    "ejp_frontier.csv, ejp_probabilistic.csv, headroom_frontier_by_duration.csv,",
    "pi_factorization_check.csv, evppi_prior_sensitivity.csv.\n")
cat("Comparability-scenario outputs (S5 horizon, S3 refractory):",
    "base_case_results_by_horizon.csv, headroom_frontier_by_horizon.csv,",
    "refractory_scenario_results.csv.\n")
cat("All tables in output/tables/.\n")
