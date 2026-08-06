# Probabilistic sensitivity analysis (analysis_plan.md §9.2, §10.2): 10,000 Monte Carlo draws
# over transition probabilities, cure parameters, costs and utilities (not just the 7 parameters
# the v6 workbook sampled — see docs/model_audit_v6.md A14). Summarise on the INMB scale, not
# per-iteration ICERs.
#
# FIRST PASS, 2026-08-05. Samples the two parameter groups that have an actual, traceable
# distribution sourced somewhere in this repo already -- the utility chain and Treg's acquisition
# price -- plus the cure fraction pi under Decision 4's own recorded uniform-prior fallback
# (analysis_plan.md §15: "the fallback is a uniform prior across the full 0-1 range for pi... That
# is less informative but entirely defensible; an invented point estimate is not"). Everything
# else this pass holds fixed at its deterministic base-case value, explicitly, rather than
# fabricating a distribution nobody has sourced yet -- see "Not yet implemented" below. This is
# the same incremental-pass discipline R/04 and R/05 were each built under.
#
# SECOND PASS, 2026-08-06 (peer review 2026-08-05, B3; docs/model_audit_v6.md A19). The post-cure
# relapse hazard h is now sampled too (point 4 below), not hardcoded to 0. **Every probabilistic
# number this module has ever produced is superseded by that change** -- see A19 and README.md's
# Status section.
#
# ---- What's sampled, and where its distribution actually comes from ------------------------------
#
# 1. **The utility chain** (Remission + the 4 multiplicative ratios; analysis_plan.md §7.1 items
#    19-20, §10.2's "must be sampled as a chain, as the current workbook correctly does"). Sourced
#    from `data/processed/model_psa_parameter_distributions.csv`'s own Uniform(lower_bound,
#    upper_bound) rows for each of Remission/Mild:Remission/M-SR:Mild/M-S:M-SR/Surgery:M-SR --
#    the exact bounds the v6 workbook itself used, not invented ones. Each is drawn independently
#    per iteration, then chained multiplicatively (Mild = Remission x ratio, M-SR = Mild x ratio,
#    M-S = M-SR x ratio, Surgery = M-SR x ratio) -- originally verified against
#    model_health_utilities.csv's own arithmetic while building this module, which is how A16
#    (docs/model_audit_v6.md) was first found: that file's Remission value, 0.9554, did not match
#    the 0.82 base value this same PSA file cites for the same quantity. **A16 is now resolved**
#    (2026-08-05) -- 0.9554 was a stray live PSA draw captured in the workbook snapshot, not a
#    parameter; model_health_utilities.csv is retired as a base-case input and
#    R/04_costs_utilities.R's load_health_state_utilities() now derives the deterministic vector
#    from this same PSA file's base values, via the identical chain this function already samples
#    around -- so the deterministic base case and this module's PSA draws are now the same central
#    value by construction, not two independently-maintained figures.
# 2. **Treg's acquisition price** (item 12). `model_psa_parameter_distributions.csv`'s
#    `psa_cost_treg_dose` row gives mean $19,917 and an `alpha_or_sd` column of 5080.87 -- verified
#    here to be a standard deviation, not a Gamma shape parameter despite the column's name:
#    method-of-moments (shape = (mean/sd)^2, scale = sd^2/mean) on (19917, 5080.87) gives shape
#    15.37, scale 1296.14, which reproduces analysis_plan.md item 12's own independently-stated
#    "Gamma (alpha ~= 15.4, scale ~= 1,297)" to the first decimal -- strong confirmation this is
#    the right reading of that column, not a guess.
# 3. **Cure fraction pi**: Uniform(0, 1), per Decision 4's own recorded fallback (above). Not this
#    pass's invention -- literally the recorded decision when elicitation didn't run.
# 4. **Post-cure relapse hazard h** (added 2026-08-06, B3): Gamma with mean fixed at
#    `duration_to_hazard(DEFAULT_RELAPSE_DURATION_YEARS)` = ln(2)/10 = 0.0693/yr -- the SAME
#    10-year median-SDR-duration anchor R/05_deterministic_results.R's headroom functions and
#    R/08_ejp.R's EJP functions now take as their deterministic default (same commit), so this
#    module's prior mean reproduces the deterministic base case by construction rather than being
#    a separately-maintained centre (the same discipline A16's resolution imposed on the utility
#    chain). Gamma, not Uniform or Beta, because h is a RATE: Briggs, Claxton & Sculpher,
#    *Decision Modelling for Health Economic Evaluation* (Oxford University Press, 2006), §4.3 --
#    Beta for a probability, Dirichlet for a multinomial transition row, Gamma for a rate or a
#    cost, all of which this module now uses in exactly those roles. Shape and the spread argument
#    are on sample_relapse_hazard_draws() below, which also states what the chosen shape implies
#    in T-years at each percentile and why that band is defensible against
#    RELAPSE_DURATION_GRID_YEARS.
#
# ---- Not yet implemented in this pass -------------------------------------------------------------
#
# - **Transition-probability (Dirichlet) sampling.** A14 (docs/model_audit_v6.md) calls this out as
#   the single biggest defect in the workbook's own PSA, and it is NOT fixed here either --
#   sample_dirichlet_row() below is a general, tested utility ready to use, but nothing calls it on
#   the real transition matrices, because no sourced concentration/precision parameter exists
#   anywhere in this repo for Aliyev's own published Table 3/4 rows (the only alpha/beta pseudo-
#   counts on file, in aliyev2019_appendixS1_table2_parameters.csv, are for the RETIRED from-
#   scratch DEALE endpoints -- R/00_derive_transition_probs.R's own history comment -- a different,
#   no-longer-used derivation path, not something that maps onto Table 3/4's rows). Wiring this in
#   needs either the underlying trial sample sizes or a documented, defensible precision assumption
#   -- neither exists yet.
# - ~~**Post-cure relapse hazard h.**~~ **Now sampled (2026-08-06, B3)** -- point 4 above. The
#   original reason for holding it fixed ("no numeric value, range, or fallback is recorded
#   anywhere for h, unlike pi") stopped being true on 2026-08-05, when the h sweep landed
#   `DEFAULT_RELAPSE_DURATION_YEARS <- 10` and `RELAPSE_DURATION_GRID_YEARS <- c(2, 5, 10, 20,
#   Inf)` (R/05_deterministic_results.R) as this project's recorded, PolTREG-anchored durability
#   assumption; this module simply never picked them up. That gap is what B3 named, and holding h
#   at 0 was not a neutral placeholder: zero relapse is the single MOST favourable assumption
#   available to the intervention (R/03_cure_fraction_module.R's duration_to_hazard() docstring),
#   so every PSA/EVPI/EVPPI/probabilistic-EJP figure produced before this change was computed
#   assuming Treg's cure is permanent for every one of its 10,000 draws.
# - **UST/IFX/ADA current drug prices** (item 14: "Re-extract... Gamma", status not yet done) and
#   **health-state monitoring/Surgery/CT-drug costs** (items 16-18): no PSA range is sourced for
#   any of these under the CURRENT CMS-pricing/Aliyev-inflation basis this project actually uses.
#   `model_dose_costs_and_psa_ranges.csv` does carry a Gamma range for an "IFX dose cost" --
#   deliberately NOT used here, for the same reason R/04's module header gives for not using that
#   file's cost figures at all: it's the old 8-week-cycle whole-dose composite, a different
#   quantity from the current $/mg CMS-ASP price this pipeline actually charges per vial-rounded
#   dose, not a re-parameterisation of it.
# - **The EJP posterior** (§10.2's fourth listed PSA output). Needs R/08_ejp.R's price-solving
#   algebra, which doesn't exist yet (R/08 is still a stub) -- psa_ejp_posterior() is not attempted
#   here; add it once R/08 lands.
# - **CEAF** (cost-effectiveness acceptability frontier). psa_ceac() below gives per-arm P(CE) by
#   WTP, which is what a CEAF is built from, but the frontier itself (the envelope + presentation
#   choice for ties) is left for whoever builds the actual figure -- a presentation decision, not
#   a computation this module is blocked on.

if (!exists("MAINTENANCE_STATES")) source("R/utils/transition_matrix.R")
if (!exists("run_base_case")) source("R/05_deterministic_results.R")

# ---- Method-of-moments Gamma parameterisation ------------------------------------------------

#' Gamma(shape, scale) via method of moments from a mean and standard deviation (mean = shape *
#' scale, variance = shape * scale^2). The standard, distribution-free way to turn "we know the
#' mean and roughly how spread out it is" into concrete Gamma parameters -- used here because
#' that's exactly the form `model_psa_parameter_distributions.csv`'s Gamma rows are recorded in
#' (a mean and an SD-like column), not as (shape, scale) directly.
gamma_shape_scale <- function(mean, sd) {
  stopifnot(mean > 0, sd > 0)
  shape <- (mean / sd)^2
  scale <- sd^2 / mean
  list(shape = shape, scale = scale)
}

# ---- Sourced PSA distributions (data/processed/model_psa_parameter_distributions.csv) ---------

#' The workbook's own PSA specification: one row per sampled parameter, with `lower_bound`/
#' `upper_bound` for Uniform rows and `alpha_or_sd` (a standard deviation -- see module header) for
#' Gamma rows.
load_psa_distributions <- function(proc_dir = "data/processed") {
  utils::read.csv(file.path(proc_dir, "model_psa_parameter_distributions.csv"), stringsAsFactors = FALSE)
}

#' One row's Uniform(lower_bound, upper_bound), by its `variable` label.
psa_uniform_bounds <- function(psa_df, variable_name) {
  row <- psa_df[psa_df$variable == variable_name, ]
  stopifnot(nrow(row) == 1, identical(row$distribution, "Uniform"))
  c(low = row$lower_bound, high = row$upper_bound)
}

#' One row's Gamma(shape, scale), by its `variable` label -- method-of-moments from `value` (the
#' mean) and `alpha_or_sd` (the SD, module header).
psa_gamma_params <- function(psa_df, variable_name) {
  row <- psa_df[psa_df$variable == variable_name, ]
  stopifnot(nrow(row) == 1, identical(row$distribution, "Gamma"))
  gamma_shape_scale(row$value, row$alpha_or_sd)
}

# ---- Utility chain draws (structurally correlated -- sampled as a chain, analysis_plan.md sec 10.2) --

#' `n_draws` independent Uniform draws for Remission and each of the 4 ratios, chained
#' multiplicatively into a full 6-state utility vector per draw (Death always 0). Returns an
#' `n_draws`-row matrix with one column per MAINTENANCE_STATES entry, ready to hand to
#' run_comparator_arm_lifetime()/run_treg_arm_lifetime()'s `utilities` override one row at a time.
sample_utility_chain_draws <- function(n_draws, proc_dir = "data/processed") {
  psa_df <- load_psa_distributions(proc_dir)

  draw_uniform <- function(variable_name) {
    b <- psa_uniform_bounds(psa_df, variable_name)
    stats::runif(n_draws, b[["low"]], b[["high"]])
  }

  remission <- draw_uniform("Remission")
  ratio_mild_remission <- draw_uniform("Mild:Remission Ratio")
  ratio_msr_mild <- draw_uniform("M-SR:Mild Ratio")
  ratio_ms_msr <- draw_uniform("M-S:M-SR Ratio")
  ratio_surgery_msr <- draw_uniform("Surgery:M-SR Ratio")

  mild <- remission * ratio_mild_remission
  msr <- mild * ratio_msr_mild
  ms <- msr * ratio_ms_msr
  surgery <- msr * ratio_surgery_msr

  matrix(
    c(ms, msr, mild, remission, surgery, rep(0, n_draws)),
    nrow = n_draws, ncol = length(MAINTENANCE_STATES),
    dimnames = list(NULL, MAINTENANCE_STATES)
  )
}

# ---- Treg acquisition price draws --------------------------------------------------------------

#' `n_draws` Gamma draws for Treg's acquisition price (module header point 2). This is the PRICE
#' component only -- run_treg_arm_lifetime()'s `price_usd` argument, which
#' treg_price_dependent_dose_cost() then adds administration/cyclophosphamide/observation-stay
#' costs on top of, exactly as the deterministic base case does.
sample_treg_price_draws <- function(n_draws, proc_dir = "data/processed") {
  psa_df <- load_psa_distributions(proc_dir)
  g <- psa_gamma_params(psa_df, "psa_cost_treg_dose")
  stats::rgamma(n_draws, shape = g$shape, scale = g$scale)
}

# ---- Post-cure relapse hazard draws (EVPPI subset B; peer review 2026-08-05's B3) ---------------

#' Shape parameter of the Gamma prior on the post-cure annual relapse hazard h. Exposed as a named
#' constant, not buried in sample_relapse_hazard_draws()'s signature, for the same reason
#' PI_PRIOR_SENSITIVITY_SPECS (R/07_evpi_evppi.R) exists: the LEVEL of a VOI result depends on the
#' prior's spread, so the spread has to be a visible, re-runnable modelling choice a reviewer can
#' vary, not a magic number. See sample_relapse_hazard_draws() for why 2.
RELAPSE_HAZARD_PSA_GAMMA_SHAPE <- 2

#' `n_draws` Gamma draws for the post-cure annual relapse hazard h -- run_treg_arm_lifetime()'s
#' `relapse_hazard_annual` argument, which R/03_cure_fraction_module.R then converts to a per-cycle
#' relapse probability. Parameterised by (mean, shape) rather than (shape, rate) because the mean
#' is the part with an actual anchor behind it and the shape is the part being ASSUMED; writing it
#' this way keeps the assumption in the signature instead of dissolved into a rate constant.
#'
#' ---- Why Gamma ------------------------------------------------------------------------------
#'
#' h is a rate (events per year), not a probability and not a multinomial row. The standard
#' health-economics PSA convention for a rate is a Gamma distribution -- Briggs, Claxton &
#' Sculpher, *Decision Modelling for Health Economic Evaluation* (Oxford University Press, 2006),
#' §4.3 -- for two reasons that both matter here: it has support on (0, Inf), matching the
#' parameter's own domain with no truncation or rejection step, and it is the conjugate form for
#' a Poisson event count, which is what a relapse process observed over person-time actually is.
#' The alternatives are wrong for the parameter type, not merely less conventional: a Beta would
#' impose an upper bound of 1 that has no meaning for a rate, and a Uniform over some [lo, hi]
#' would need two invented endpoints instead of one invented spread parameter.
#'
#' ---- Why the mean is fixed at ln(2)/10 -------------------------------------------------------
#'
#' `mean_hazard` defaults to `duration_to_hazard(DEFAULT_RELAPSE_DURATION_YEARS)` = ln(2)/10 =
#' 0.06931/yr, i.e. a median Sustained Deep Remission duration of 10 years. That is this project's
#' already-recorded base-case anchor (R/05_deterministic_results.R's own comment on
#' RELAPSE_DURATION_GRID_YEARS: the PolTREG T1D cohort supports a multi-year plateau in a subset
#' out to 7-12 years, but nothing published supports permanence), and as of the same commit it is
#' also the deterministic default in headroom_pi_star()/headroom_frontier() (R/05) and
#' ejp_deterministic()/ejp_frontier() (R/08). Fixing the PSA's MEAN there -- not its median or its
#' mode -- makes the probabilistic analysis the stochastic counterpart of the deterministic one by
#' construction, which is exactly the property A6 (docs/model_audit_v6.md) found the v6 workbook's
#' own PSA lacked. No new number is introduced by this function: the only thing being assumed here
#' that wasn't already assumed somewhere in this repository is the shape.
#'
#' ---- Why shape = 2, and what it implies in years ---------------------------------------------
#'
#' Shape k fixes the coefficient of variation at 1/sqrt(k) once the mean is pinned (Gamma sd =
#' mean/sqrt(k)); k = 2 gives CV = 0.71. Translated back into the units this project actually
#' reasons about via hazard_to_duration(), T = ln(2)/h, the implied median SDR duration at each
#' percentile of the sampled h is:
#'
#'   h percentile:   2.5th    25th    50th    75th    97.5th   99.5th   99.9th
#'   implied T:      82.6yr   20.8yr  11.9yr  7.4yr   3.6yr    2.7yr    2.2yr
#'
#' (h's LOW percentiles are T's HIGH ones -- a small hazard is a long remission.) Read against
#' RELAPSE_DURATION_GRID_YEARS = c(2, 5, 10, 20, Inf), the deterministic sweep this project
#' already committed to: the interquartile range [7.4, 20.8] brackets the 10-year anchor and its
#' two immediate neighbours on the grid; the sweep's pessimistic endpoint T = 2 sits at the
#' 99.95th percentile of h, so a 10,000-draw PSA draws roughly 5 iterations more pessimistic than
#' the pessimistic end of the sweep -- present in the tail, as it should be, without the prior
#' pretending 2-year durability is a central expectation. The median of 11.9 years is inside
#' PolTREG's own 7-12 year observation window, which is the closest thing to an external anchor
#' this parameter has (tests/testthat/test-external-validity.R asserts exactly that).
#'
#' Shape > 1 is doing real work at the other end and is not just a spread choice: the Gamma density
#' vanishes at h = 0 for any k > 1, so T = Inf (permanent cure, the Ovasave/CATS1-style upper
#' bound) has probability zero under this prior. An Exponential (k = 1) would instead put its
#' MODAL density at h = 0 -- i.e. treat near-permanence as the single most likely outcome, which is
#' precisely the anti-conservative framing `docs/treg-cd_decision_resolutions_2026-08-05.md` §3.3
#' corrected and B3 found still live in this module.
#'
#' ---- The tension this choice cannot resolve, stated rather than hidden ------------------------
#'
#' 26% of draws imply T > 20 years, more optimistic-tail mass than one would ideally assign given
#' that nothing published supports permanence. That is not a tuning failure; it is structural. In
#' hazard space the sweep grid is badly asymmetric about its own anchor -- T = 2 is 5x the anchor
#' hazard while T = 20 is only 0.5x it -- so NO distribution with mean h0 can put T = 2 at a
#' conventional upper tail percentile without simultaneously putting at least a quarter of its mass
#' below h0/2. Tightening the shape buys a smaller optimistic tail only by abandoning the
#' pessimistic endpoint entirely (at k = 4, T = 2 sits beyond the 99.99th percentile and is never
#' drawn at 10,000 iterations). k = 2 is the deliberate compromise: keep the pessimistic endpoint
#' inside the sampled support, accept the optimistic tail, and continue to report the DETERMINISTIC
#' T sweep (headroom_frontier_by_duration(), R/05) as the primary, prior-free presentation of
#' durability uncertainty -- the PSA characterises it, the sweep displays it. The same
#' level-vs-ranking caveat R/07_evpi_evppi.R's evppi_prior_sensitivity() states for pi applies
#' verbatim to subset B: the ranking of subsets is the finding, the level is prior-dependent.
sample_relapse_hazard_draws <- function(n_draws,
                                         mean_hazard = duration_to_hazard(DEFAULT_RELAPSE_DURATION_YEARS),
                                         shape = RELAPSE_HAZARD_PSA_GAMMA_SHAPE) {
  stopifnot(n_draws >= 0, mean_hazard > 0, shape > 0)
  stats::rgamma(n_draws, shape = shape, scale = mean_hazard / shape)
}

# ---- Dirichlet sampling (general utility, not yet wired to any real matrix -- see module header) --

#' Sample one row of a transition matrix from a Dirichlet distribution with mean `probs` and
#' total concentration `concentration` (larger = tighter around `probs`; concentration = sum of
#' the underlying pseudo-counts, e.g. a trial's sample size, if one is available). Built via
#' independent Gamma draws normalised to sum to 1 -- the standard construction, and the one that
#' needs no package beyond base R's `rgamma()`. Returns an `n_draws` x `length(probs)` matrix,
#' each row summing to exactly 1 by construction (Gamma draws are always > 0, so no zero-sum edge
#' case to guard against).
#'
#' NOT called anywhere in this module yet -- see "Not yet implemented" above for why (no sourced
#' concentration parameter for the real transition matrices). Kept here, tested, so wiring it in
#' later is a matter of supplying real `probs`/`concentration`, not writing this function under
#' time pressure once the missing data shows up.
sample_dirichlet_row <- function(n_draws, probs, concentration) {
  stopifnot(abs(sum(probs) - 1) < 1e-6, concentration > 0, all(probs >= 0))
  alpha <- probs * concentration
  alpha[alpha == 0] <- 1e-8  # a structurally-zero transition stays ~0 without a degenerate Gamma(0,.) draw
  raw <- vapply(alpha, function(a) stats::rgamma(n_draws, shape = a, scale = 1), numeric(n_draws))
  raw / rowSums(raw)
}

# ---- Orchestration ------------------------------------------------------------------------------

#' Run the full PSA: `n_draws` iterations, each sampling the utility chain, Treg's price, pi
#' (Uniform(0,1), Decision 4's fallback -- module header) and the post-cure relapse hazard h
#' (Gamma, sample_relapse_hazard_draws() -- added 2026-08-06, B3), then running all four arms
#' (UST/IFX/ADA/CT-backed comparators plus TREG) through R/05_deterministic_results.R's arm
#' runners with that draw's utility vector applied consistently across every arm. Everything else
#' (transition matrices, non-Treg costs) is held fixed at its deterministic base-case value -- see
#' module header for why each is deferred rather than sampled.
#'
#' The four sampled quantities are drawn in the order utility chain -> price -> pi -> h, and h was
#' appended LAST deliberately: for any fixed `seed`, every draw of the three parameters that were
#' already sampled before B3 is bit-for-bit what it was before h existed, so any difference between
#' a pre-B3 and post-B3 run at the same seed is attributable to h alone and not to the RNG stream
#' having been reshuffled underneath it.
#'
#' `seed` defaults to a fixed value for reproducibility -- analysis_plan.md doesn't specify one,
#' so this is a documented modelling choice, not an attempt to hide seed-sensitivity; pass `NULL`
#' to disable and get a fresh draw set each call.
#'
#' Returns one row per (draw, arm): draw index, intervention, qalys, total_cost, that draw's
#' sampled pi, treg_price and relapse_hazard_annual (all three NA for the three comparators, which
#' don't depend on any of them -- the relapse hazard reaches the model only through Treg's own
#' SDR track, R/03_cure_fraction_module.R, so a comparator arm's trace is untouched by it), and
#' that draw's 5 raw sampled utility-chain values (util_modsev/util_resp/util_mild/
#' util_remission/util_surgery, same variable_name convention as
#' data/processed/model_health_utilities.csv) -- recorded on EVERY row, comparators included,
#' since utility applies identically to every arm within a draw, not just Treg.
#' R/07_evpi_evppi.R's subset-E EVPPI needs these raw values directly; qalys/total_cost alone
#' don't let you regress net benefit against "how favourable was this draw's utility draw",
#' only against the CONSEQUENCE of that draw, which is exactly the effect you're trying to
#' explain, not a usable predictor.
#'
#' `n_cycles`/`baseline_age`/`life_table = NULL` (defaults, byte-identical to before these
#' comments existed): `n_cycles` still defaults to HORIZON_CYCLES_6YR and `baseline_age = NULL`
#' still means every arm runs on a fixed matrix with Aliyev's own embedded trial-cohort
#' mortality, exactly as this function always has. `analysis/run_full_analysis.R`'s primary PSA
#' now calls this with `n_cycles = HORIZON_CYCLES_LIFETIME, baseline_age =
#' ASSUMED_PATIENT_AGE_YEARS` (peer review 2026-08-05, B2/README.md's Status section) -- the
#' 6.15-year defaults below remain available for the S5 comparability-scenario PSA, unchanged.
#'
#' **Performance note (2026-08-06, same review):** the three comparator arms' occupancy traces
#' (UST/IFX/ADA) are simulated ONCE, before the draw loop, via
#' simulate_comparator_arm_lifetime() (R/05_deterministic_results.R) -- not once per draw. This
#' is exact, not an approximation: R/06 doesn't sample any transition probability (module
#' header's "Not yet implemented" list), so a comparator's induction split and maintenance matrix
#' -- and therefore its whole trace -- is identical across every draw for a fixed (n_cycles,
#' baseline_age) call; only that draw's utility vector changes what
#' attach_maintenance_costs_utilities() computes FROM the (unchanged) trace. This collapses what
#' would otherwise be `n_draws` x 3 full lifetime Markov simulations (the dominant cost of a
#' lifetime-horizon PSA, since each one re-derives an age-adjusted matrix pair every one of
#' HORIZON_CYCLES_LIFETIME cycles once `baseline_age` is supplied) down to 3, regardless of
#' `n_draws`. TREG is NOT hoisted this way: its occupancy trace genuinely varies per draw (pi and,
#' since 2026-08-06, h both vary), so it still simulates once per draw.
#'
#' **The pi-interpolation optimisation this used to flag is no longer available, and that is worth
#' recording rather than silently dropping.** `verify_pi_factorization()`
#' (R/05_deterministic_results.R) establishes that Treg's QALY gain is exactly linear in pi AT A
#' FIXED h -- g(h) is a per-h constant, so the shortcut was "simulate pi=0 and pi=1 once, then
#' interpolate every draw." Sampling h breaks that: each draw now needs its own g(h), so the
#' reference runs would have to be recomputed per distinct h, of which there are `n_draws`. The
#' factorisation itself is unaffected (it was never a claim about h), only its usefulness as a
#' shortcut here. A future pass wanting it back would have to interpolate g(h) over an h grid and
#' demonstrate the interpolation error is negligible against the QALY differences being reported --
#' a real piece of numerical work, not a code tidy-up.
run_psa <- function(n_draws = 10000, n_cycles = HORIZON_CYCLES_6YR, weight_kg = ASSUMED_PATIENT_WEIGHT_KG,
                     cycle_weeks = 2, annual_rate = 0.03, apply_cap = TRUE, cap_cycle = 52,
                     raw_dir = "data/raw", proc_dir = "data/processed", seed = 20260805,
                     baseline_age = NULL, life_table = NULL) {
  if (!is.null(seed)) set.seed(seed)

  matrices <- build_all_transition_matrices(raw_dir)
  utility_draws <- sample_utility_chain_draws(n_draws, proc_dir)
  price_draws <- sample_treg_price_draws(n_draws, proc_dir)
  pi_draws <- stats::runif(n_draws, 0, 1)
  # Drawn last on purpose -- see this function's own docstring on why the RNG order matters.
  relapse_hazard_draws <- sample_relapse_hazard_draws(n_draws)

  # Load once, reuse across all n_draws x 4-arm calls below (run_comparator_arm_lifetime()'s and
  # run_treg_arm_lifetime()'s own caching parameters, added alongside this module) -- none of
  # induction_data/schedule/prices varies across PSA draws in this pass, so re-reading them from
  # disk on every one of the ~40,000 arm-runner calls a 10,000-draw PSA makes would be pure
  # overhead, not a correctness requirement.
  induction_data <- load_published_induction(raw_dir)
  schedule <- load_dosing_schedule(raw_dir)
  prices <- load_drug_prices(raw_dir)
  if (is.null(life_table) && !is.null(baseline_age)) life_table <- load_life_table(raw_dir)

  # Comparator-trace hoist (this function's own docstring, "Performance note"): simulated ONCE
  # per therapy, outside the draw loop, since nothing PSA samples ever reaches this step.
  # `monitoring_costs`/`induction_costs` are likewise draw-invariant (perspective is always
  # "healthcare_sector" in this pass -- no scenario wrapper here varies it) and hoisted the same
  # way, for the same reason.
  comparator_arms <- stats::setNames(
    lapply(COMPARATOR_THERAPIES, function(tx) {
      simulate_comparator_arm_lifetime(
        tx, n_cycles, matrices, cap_cycle = cap_cycle, apply_cap = apply_cap, raw_dir = raw_dir,
        induction_data = induction_data, baseline_age = baseline_age, life_table = life_table,
        cycle_weeks = cycle_weeks
      )
    }),
    COMPARATOR_THERAPIES
  )
  monitoring_costs <- health_state_monitoring_costs()
  induction_costs <- stats::setNames(
    vapply(COMPARATOR_THERAPIES, induction_drug_cost, numeric(1), weight_kg = weight_kg,
           schedule = schedule, prices = prices),
    COMPARATOR_THERAPIES
  )

  rows <- vector("list", n_draws * (length(COMPARATOR_THERAPIES) + 1))
  k <- 0
  for (i in seq_len(n_draws)) {
    utilities_i <- utility_draws[i, ]
    util_cols <- list(
      util_modsev = unname(utilities_i[["Moderate-Severe"]]),
      util_resp = unname(utilities_i[["Moderate-Severe Responder"]]),
      util_mild = unname(utilities_i[["Mild"]]),
      util_remission = unname(utilities_i[["Remission"]]),
      util_surgery = unname(utilities_i[["Surgery"]])
    )

    for (tx in COMPARATOR_THERAPIES) {
      # Re-aggregates the SAME hoisted trace with this draw's own utility vector -- exactly what
      # run_comparator_arm_lifetime() computes downstream of its own (here-skipped) simulate step.
      attached <- attach_maintenance_costs_utilities(
        comparator_arms[[tx]], utilities_i, monitoring_costs, cycle_weeks, annual_rate,
        therapy = tx, weight_kg = weight_kg, schedule = schedule, prices = prices
      )
      s <- summarise_arm(attached, induction_cost = induction_costs[[tx]])
      k <- k + 1
      rows[[k]] <- c(list(draw = i, intervention = tx, qalys = s$qalys, total_cost = s$total_cost,
                           pi_sdr = NA_real_, treg_price = NA_real_,
                           relapse_hazard_annual = NA_real_), util_cols)
    }

    treg <- run_treg_arm_lifetime(n_cycles, pi_draws[i], relapse_hazard_draws[i], price_draws[i],
                                   matrices, weight_kg, cycle_weeks, annual_rate, cap_cycle = cap_cycle,
                                   apply_cap = apply_cap, raw_dir = raw_dir, proc_dir = proc_dir,
                                   utilities = utilities_i, induction_data = induction_data, prices = prices,
                                   baseline_age = baseline_age, life_table = life_table)
    k <- k + 1
    rows[[k]] <- c(list(draw = i, intervention = "TREG", qalys = treg$qalys, total_cost = treg$total_cost,
                         pi_sdr = pi_draws[i], treg_price = price_draws[i],
                         relapse_hazard_annual = relapse_hazard_draws[i]), util_cols)
  }

  do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE))
}

# ---- Summaries: cost-effectiveness plane and CEAC ------------------------------------------------

#' Incremental cost/QALY draws for `intervention` vs. `comparator`, one row per PSA draw --
#' exactly the data a cost-effectiveness-plane scatterplot needs (analysis_plan.md sec 10.2).
psa_cost_effectiveness_plane <- function(psa_results, intervention, comparator) {
  int_rows <- psa_results[psa_results$intervention == intervention, ]
  comp_rows <- psa_results[psa_results$intervention == comparator, ]
  stopifnot(nrow(int_rows) == nrow(comp_rows))
  int_rows <- int_rows[order(int_rows$draw), ]
  comp_rows <- comp_rows[order(comp_rows$draw), ]
  data.frame(
    draw = int_rows$draw,
    delta_qalys = int_rows$qalys - comp_rows$qalys,
    delta_cost = int_rows$total_cost - comp_rows$total_cost
  )
}

#' Cost-effectiveness acceptability curve: for each arm and each WTP threshold, the fraction of
#' PSA draws where that arm has the single highest NMB among ALL arms present in `psa_results` --
#' the standard CEAC definition, computed once per WTP value in `wtp_grid` (analysis_plan.md
#' sec 9.2's hygiene note: computed on the NMB scale, never by averaging ICERs).
psa_ceac <- function(psa_results, wtp_grid) {
  draws <- sort(unique(psa_results$draw))
  arms <- unique(psa_results$intervention)

  rows <- lapply(wtp_grid, function(wtp) {
    winners <- vapply(draws, function(d) {
      draw_rows <- psa_results[psa_results$draw == d, ]
      nmb <- draw_rows$qalys * wtp - draw_rows$total_cost
      draw_rows$intervention[which.max(nmb)]
    }, character(1))
    counts <- table(factor(winners, levels = arms))
    stats::setNames(as.list(as.numeric(counts) / length(draws)), arms)
  })

  cbind(wtp_usd = wtp_grid, do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE)))
}
