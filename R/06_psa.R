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
# ---- What's sampled, and where its distribution actually comes from ------------------------------
#
# 1. **The utility chain** (Remission + the 4 multiplicative ratios; analysis_plan.md §7.1 items
#    19-20, §10.2's "must be sampled as a chain, as the current workbook correctly does"). Sourced
#    from `data/processed/model_psa_parameter_distributions.csv`'s own Uniform(lower_bound,
#    upper_bound) rows for each of Remission/Mild:Remission/M-SR:Mild/M-S:M-SR/Surgery:M-SR --
#    the exact bounds the v6 workbook itself used, not invented ones. Each is drawn independently
#    per iteration, then chained multiplicatively (Mild = Remission x ratio, M-SR = Mild x ratio,
#    M-S = M-SR x ratio, Surgery = M-SR x ratio -- verified against model_health_utilities.csv's
#    own arithmetic; see A16, docs/model_audit_v6.md, for a real discrepancy found while doing
#    this verification: model_health_utilities.csv's actual Remission value, 0.9554, does not
#    match the 0.82 base value this same PSA file cites for the same quantity. That does NOT
#    block sampling here -- a Uniform only needs correct bounds, not a correct centre -- but it
#    is a real, unresolved gap between the deterministic base case and the cited literature
#    source, flagged, not silently picked around.
# 2. **Treg's acquisition price** (item 12). `model_psa_parameter_distributions.csv`'s
#    `psa_cost_treg_dose` row gives mean $19,917 and an `alpha_or_sd` column of 5080.87 -- verified
#    here to be a standard deviation, not a Gamma shape parameter despite the column's name:
#    method-of-moments (shape = (mean/sd)^2, scale = sd^2/mean) on (19917, 5080.87) gives shape
#    15.37, scale 1296.14, which reproduces analysis_plan.md item 12's own independently-stated
#    "Gamma (alpha ~= 15.4, scale ~= 1,297)" to the first decimal -- strong confirmation this is
#    the right reading of that column, not a guess.
# 3. **Cure fraction pi**: Uniform(0, 1), per Decision 4's own recorded fallback (above). Not this
#    pass's invention -- literally the recorded decision when elicitation didn't run.
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
# - **Post-cure relapse hazard h.** No numeric value, range, or fallback is recorded anywhere for
#   h (unlike pi, which has Decision 4's explicit fallback) -- held fixed at 0 in this pass, same
#   choice R/05_deterministic_results.R's headroom functions already make, for the same reason.
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

#' Run the full PSA: `n_draws` iterations, each sampling the utility chain, Treg's price, and pi
#' (Uniform(0,1), Decision 4's fallback -- module header), then running all four arms
#' (UST/IFX/ADA/CT-backed comparators plus TREG) through R/05_deterministic_results.R's arm
#' runners with that draw's utility vector applied consistently across every arm. Everything else
#' (transition matrices, relapse hazard h = 0, non-Treg costs) is held fixed at its deterministic
#' base-case value -- see module header for why each is deferred rather than sampled.
#'
#' `seed` defaults to a fixed value for reproducibility -- analysis_plan.md doesn't specify one,
#' so this is a documented modelling choice, not an attempt to hide seed-sensitivity; pass `NULL`
#' to disable and get a fresh draw set each call.
#'
#' Returns one row per (draw, arm): draw index, intervention, qalys, total_cost, plus that draw's
#' sampled pi and treg_price (NA for the three comparators, which don't depend on either).
run_psa <- function(n_draws = 10000, n_cycles = HORIZON_CYCLES_6YR, weight_kg = ASSUMED_PATIENT_WEIGHT_KG,
                     cycle_weeks = 2, annual_rate = 0.03, apply_cap = TRUE, cap_cycle = 52,
                     raw_dir = "data/raw", proc_dir = "data/processed", seed = 20260805) {
  if (!is.null(seed)) set.seed(seed)

  matrices <- build_all_transition_matrices(raw_dir)
  utility_draws <- sample_utility_chain_draws(n_draws, proc_dir)
  price_draws <- sample_treg_price_draws(n_draws, proc_dir)
  pi_draws <- stats::runif(n_draws, 0, 1)

  # Load once, reuse across all n_draws x 4-arm calls below (run_comparator_arm_lifetime()'s and
  # run_treg_arm_lifetime()'s own caching parameters, added alongside this module) -- none of
  # induction_data/schedule/prices varies across PSA draws in this pass, so re-reading them from
  # disk on every one of the ~40,000 arm-runner calls a 10,000-draw PSA makes would be pure
  # overhead, not a correctness requirement.
  induction_data <- load_published_induction(raw_dir)
  schedule <- load_dosing_schedule(raw_dir)
  prices <- load_drug_prices(raw_dir)

  rows <- vector("list", n_draws * (length(COMPARATOR_THERAPIES) + 1))
  k <- 0
  for (i in seq_len(n_draws)) {
    utilities_i <- utility_draws[i, ]

    for (tx in COMPARATOR_THERAPIES) {
      s <- run_comparator_arm_lifetime(tx, n_cycles, matrices, weight_kg, cycle_weeks, annual_rate,
                                        apply_cap, cap_cycle, raw_dir, proc_dir, utilities = utilities_i,
                                        induction_data = induction_data, schedule = schedule, prices = prices)
      k <- k + 1
      rows[[k]] <- data.frame(draw = i, intervention = tx, qalys = s$qalys, total_cost = s$total_cost,
                               pi_sdr = NA_real_, treg_price = NA_real_)
    }

    treg <- run_treg_arm_lifetime(n_cycles, pi_draws[i], relapse_hazard_annual = 0, price_draws[i],
                                   matrices, weight_kg, cycle_weeks, annual_rate, cap_cycle = cap_cycle,
                                   apply_cap = apply_cap, raw_dir = raw_dir, proc_dir = proc_dir,
                                   utilities = utilities_i, induction_data = induction_data, prices = prices)
    k <- k + 1
    rows[[k]] <- data.frame(draw = i, intervention = "TREG", qalys = treg$qalys, total_cost = treg$total_cost,
                             pi_sdr = pi_draws[i], treg_price = price_draws[i])
  }

  do.call(rbind, rows)
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
