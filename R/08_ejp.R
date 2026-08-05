# Economically justifiable price, deterministic and probabilistic, plus the headroom frontier
# (analysis_plan.md §9.1, Aim 1 and Aim 4). Report per-dose price P*, not total cost (see
# docs/model_audit_v6.md A7 for why the v6 EJP sheet is disconnected from the rest of the model).
#
# FIRST PASS, 2026-08-05.
#
# ---- The closed-form algebra (analysis_plan.md sec 9.1) -----------------------------------------
#
# INMB(P*) = 0 at the next-best non-dominated comparator gives:
#
#   P* = [ C_comp + lambda*(Q_treg - Q_comp) - C_treg,non-drug ] / D-bar
#
# D-bar (discounted expected doses per patient) = 1 under the single-dose base case this whole
# pipeline implements (R/04_costs_utilities.R's treg_dose_cost(), R/05_deterministic_results.R's
# run_treg_arm_lifetime() -- both enforce n_doses == 1; the 2-dose structural scenario needs per-
# cycle-trace-level logic neither has, same restriction, same reason, not re-litigated here).
# D_BAR below is hardcoded to 1 for exactly this reason, not left as a free parameter nobody has
# wired a real value for.
#
# Because Treg's total cost is EXACTLY linear in price (C_treg(P) = C_treg,non-drug + D-bar*P --
# treg_price_dependent_dose_cost(), R/05), C_treg,non-drug is available for free by calling
# run_treg_arm_lifetime() with price_usd = 0: whatever total_cost comes back at zero price IS
# C_treg,non-drug exactly, no arithmetic needed to back it out. This also means P* is a closed-form
# algebraic solve, not a root-find -- unlike R/05's headroom_pi_star(), which DOES need
# uniroot() because Q_treg is a genuinely nonlinear function of pi (the cure-branching mechanic),
# not of price.
#
# The comparator is identified via R/05's best_comparator_nmb() -- the same NMB-maximisation
# simplification for "next-best non-dominated option" that function's own docstring already
# earmarks for this module ("R/08_ejp.R's own comparator identification for solving P* should use
# this same function, not a separate implementation").
#
# ---- ejp_frontier() and R/05's headroom_frontier() are the same curve, traced from two sides ----
#
# R/05_deterministic_results.R's headroom_frontier() answers "at this price, what pi is required?"
# ejp_frontier() below answers "at this pi, what price is justified?" Both trace the same (pi,
# price) indifference curve INMB = 0 -- ejp_deterministic(pi, wtp) and headroom_pi_star(P*, wtp)
# should round-trip (feed one's output into the other and recover the input), which is exactly
# the cross-check tests/testthat/test-ejp.R uses to validate both implementations against each
# other, not just against hand-computed arithmetic.
#
# ---- Deterministic EJP has no single "the" pi to report at (same reason as R/05) ----------------
#
# Decision 4 (analysis_plan.md sec 15) means there is no numeric base-case pi -- ejp_deterministic()
# therefore takes `pi_sdr` as a required argument, same as headroom_pi_star() takes `price_usd`.
# The only pi value with an actual recorded number behind it is the Null scenario, pi = 0
# (analysis_plan.md sec 7.2) -- callers wanting "the" deterministic EJP for write-up purposes
# should pass that explicitly, not rely on a default this module would otherwise be inventing.
#
# ---- Probabilistic EJP: computed by algebra over R/06_psa.R's existing draws, no new simulation --
#
# ejp_probabilistic() solves the SAME closed-form P* per PSA draw, reusing R/06_psa.R's psa_results
# columns directly rather than re-running the Markov engine: C_treg,non-drug for draw k is
# total_cost_k - treg_price_k (the same linear-in-price algebra R/07_evpi_evppi.R's
# net_benefit_matrix() already exploits for its own price-override), and the next-best comparator
# for draw k is whichever of that draw's own UST/IFX/ADA rows has the highest NMB at lambda -- both
# quantities the PSA already computed, at zero extra simulation cost.
#
# analysis_plan.md sec 9.1 asks for two DISTINCT quantities, both implemented here independently
# rather than one derived from the other, per its own instruction ("these coincide only under
# symmetry"): the median (and 95% CI) of the P*_k density itself (ejp_probabilistic()'s own
# summary), and P_50 -- the specific price at which the fraction of draws still cost-effective
# equals 50% (ejp_p50(), built from the empirical survival function of P*_k, not from median()).
# Mathematically, for a continuous P*_k distribution these ARE identical (INMB(P) is strictly
# decreasing in P for every draw, so "cost-effective at P" is exactly "P <= P*_k", making
# P(cost-effective at P) = P(P*_k >= P), which crosses 0.5 exactly at P*_k's own median) -- but
# this project's per-draw comparator can SWITCH (whichever of UST/IFX/ADA has the best NMB varies
# draw to draw), which can introduce small kinks/ties a finite sample might expose. Reported and
# tested as two independently-computed numbers, expected to closely agree, not collapsed into one
# implementation on the strength of the theory alone.
#
# ---- Gross margin over COGS (analysis_plan.md sec 9.1, added 2026-08-05) -------------------------
#
# gross_margin_over_cogs(): (P* - COGS) / P*, COGS read directly from
# data/processed/model_tenham_derived_treg_dose_cost.csv's own "Total manufacturing (COGS) cost
# per treatment/dose" line ($4,979.19) -- NOT the $19,916.75 acquisition-cost assumption
# treg_dose_cost() uses, which is already COGS marked up 4x; this is a different, EJP-anchored
# markup question. EXPLICIT SCOPE LIMIT, carried forward verbatim from the plan text: this is NOT
# a manufacturer-profitability claim. It excludes R&D recoupment (including portfolio-level
# attrition), commercial infrastructure, cost of capital, and time-to-market risk -- state this
# limitation next to the figure itself whenever it's used in the manuscript, not only in a
# footnote. data/raw/market_comparator_cell_therapy_prices.csv (Ryoncil, tab-cel/Ebvallo) is
# Discussion-section context for this same figure -- reading and presenting it is a write-up task,
# not something this module computes.

if (!exists("MAINTENANCE_STATES")) source("R/utils/transition_matrix.R")
if (!exists("run_base_case")) source("R/05_deterministic_results.R")

#' Discounted expected doses per patient -- see module header for why this is 1, not a free
#' parameter, under every function in this file.
D_BAR <- 1

# ---- COGS (data/processed/model_tenham_derived_treg_dose_cost.csv) ------------------------------

#' Treg's manufacturing (COGS) cost per dose -- the pre-markup figure, read directly from the same
#' file load_treg_dose_acquisition_cost() (R/04_costs_utilities.R) reads its own (post-markup,
#' retail) line from, so the two stay reconciled to a single source file rather than one of them
#' drifting into a hardcoded number.
load_treg_cogs <- function(proc_dir = "data/processed") {
  lines <- readLines(file.path(proc_dir, "model_tenham_derived_treg_dose_cost.csv"))
  blank <- which(lines == "")
  stopifnot(length(blank) >= 1)
  second_block <- utils::read.csv(
    text = paste(lines[(blank[1] + 1):length(lines)], collapse = "\n"), stringsAsFactors = FALSE
  )
  val <- second_block$value_usd[
    second_block$final_buildup_item == "Total manufacturing (COGS) cost per treatment/dose"
  ]
  stopifnot(length(val) == 1)
  as.numeric(val)
}

#' Gross margin over COGS at a given price (module header: NOT a profitability claim -- state the
#' scope limit alongside this figure wherever it's reported).
gross_margin_over_cogs <- function(p_star, proc_dir = "data/processed") {
  cogs <- load_treg_cogs(proc_dir)
  (p_star - cogs) / p_star
}

# ---- Deterministic EJP (closed-form; module header) ----------------------------------------------

#' Solve for P* at a given, explicit pi (module header: no default -- Decision 4 leaves no numeric
#' base-case pi to fall back to). `matrices`/`comparator_summaries`, if NULL, are built once here;
#' pass precomputed ones (same objects R/05_deterministic_results.R's own functions use) when
#' calling this repeatedly, e.g. from ejp_frontier(), to avoid rebuilding them per pi.
#'
#' Returns the solved price plus every quantity the algebra used, so a caller can verify
#' INMB(p_star) == 0 directly rather than trusting the arithmetic blindly.
#'
#' `baseline_age`/`life_table = NULL` (default): same meaning as
#' R/05_deterministic_results.R's headroom_pi_star()'s equivalent parameters, threaded through
#' to both the comparator summaries and the Treg run -- preserves the round-trip cross-check
#' against headroom_pi_star() (test-ejp.R) at the lifetime horizon too, not just the two
#' pre-lifetime-horizon comparability scenarios.
ejp_deterministic <- function(pi_sdr, wtp_usd, relapse_hazard_annual = 0, n_cycles = HORIZON_CYCLES_6YR,
                               matrices = NULL, comparator_summaries = NULL,
                               weight_kg = ASSUMED_PATIENT_WEIGHT_KG, cycle_weeks = 2,
                               annual_rate = 0.03, apply_cap = TRUE, cap_cycle = 52,
                               raw_dir = "data/raw", proc_dir = "data/processed",
                               baseline_age = NULL, life_table = NULL) {
  if (is.null(matrices)) matrices <- build_all_transition_matrices(raw_dir)
  if (is.null(comparator_summaries)) {
    comparator_summaries <- stats::setNames(
      lapply(COMPARATOR_THERAPIES, function(tx) {
        run_comparator_arm_lifetime(tx, n_cycles, matrices, weight_kg, cycle_weeks, annual_rate,
                                     apply_cap, cap_cycle, raw_dir, proc_dir,
                                     baseline_age = baseline_age, life_table = life_table)
      }),
      COMPARATOR_THERAPIES
    )
  }

  # price_usd = 0: treg_price_dependent_dose_cost()'s exact linearity in price means whatever
  # total_cost comes back here IS C_treg,non-drug, with no need to back it out algebraically.
  treg <- run_treg_arm_lifetime(n_cycles, pi_sdr, relapse_hazard_annual, price_usd = 0, matrices,
                                 weight_kg, cycle_weeks, annual_rate, cap_cycle = cap_cycle,
                                 apply_cap = apply_cap, raw_dir = raw_dir, proc_dir = proc_dir,
                                 baseline_age = baseline_age, life_table = life_table)

  best <- best_comparator_nmb(comparator_summaries, wtp_usd)
  comp <- comparator_summaries[[best$comparator]]

  # analysis_plan.md's own notation calls this term "C_treg,non-drug" -- but per THIS codebase's
  # cost decomposition (summarise_arm(), R/04_costs_utilities.R) that is NOT the `non_drug_cost`
  # field: treg_price_dependent_dose_cost() always adds the $57.90 IV administration fee (plus
  # cyclophosphamide/observation-stay if supplied) regardless of price, and summarise_arm() folds
  # that into `drug_cost`/`total_cost`, not `non_drug_cost`. A real bug caught while validating
  # this function against R/05's headroom_pi_star() (round-trip test): using `non_drug_cost` here
  # silently dropped the $57.90 admin fee, shifting every P* by exactly that amount. `total_cost`
  # AT price_usd = 0 is the correct quantity -- by treg_price_dependent_dose_cost()'s exact
  # linearity, it equals the whole "everything except D-bar*P" term, administration fee included.
  p_star <- (comp$total_cost + wtp_usd * (treg$qalys - comp$qalys) - treg$total_cost) / D_BAR

  list(p_star = p_star, comparator = best$comparator, q_treg = treg$qalys,
       c_treg_ex_price = treg$total_cost, q_comp = comp$qalys, c_comp = comp$total_cost,
       wtp_usd = wtp_usd, pi_sdr = pi_sdr)
}

#' The (pi, price) indifference curve, swept over `pi_grid` -- the same curve
#' R/05_deterministic_results.R's headroom_frontier() traces from the price side (module header).
#' Builds matrices/comparator summaries ONCE for the whole grid, same pattern as
#' headroom_frontier() itself.
#' `baseline_age`/`life_table = NULL` (default): same meaning as ejp_deterministic()'s
#' equivalent parameters, threaded through unchanged.
ejp_frontier <- function(pi_grid, wtp_usd, relapse_hazard_annual = 0, n_cycles = HORIZON_CYCLES_6YR,
                          weight_kg = ASSUMED_PATIENT_WEIGHT_KG, cycle_weeks = 2, annual_rate = 0.03,
                          apply_cap = TRUE, cap_cycle = 52, raw_dir = "data/raw",
                          proc_dir = "data/processed", baseline_age = NULL, life_table = NULL) {
  matrices <- build_all_transition_matrices(raw_dir)
  comparator_summaries <- stats::setNames(
    lapply(COMPARATOR_THERAPIES, function(tx) {
      run_comparator_arm_lifetime(tx, n_cycles, matrices, weight_kg, cycle_weeks, annual_rate,
                                   apply_cap, cap_cycle, raw_dir, proc_dir,
                                   baseline_age = baseline_age, life_table = life_table)
    }),
    COMPARATOR_THERAPIES
  )

  rows <- lapply(pi_grid, function(pi_sdr) {
    res <- ejp_deterministic(pi_sdr, wtp_usd, relapse_hazard_annual, n_cycles, matrices,
                              comparator_summaries, weight_kg, cycle_weeks, annual_rate, apply_cap,
                              cap_cycle, raw_dir, proc_dir, baseline_age, life_table)
    data.frame(pi_sdr = pi_sdr, wtp_usd = wtp_usd, p_star = res$p_star, comparator = res$comparator)
  })
  do.call(rbind, rows)
}

# ---- Probabilistic EJP (algebra over R/06_psa.R's existing draws; module header) -----------------

#' P*_k for every PSA draw, at a fixed WTP threshold -- pure post-hoc algebra over
#' `psa_results` (R/06_psa.R's run_psa() output), no new Markov simulation (module header).
#' Returns one row per draw: draw index, p_star, and which comparator was "next-best" for that
#' draw (which can, and does, vary draw to draw -- module header's note on why P_50 and the
#' median aren't collapsed into one implementation).
ejp_probabilistic_draws <- function(psa_results, wtp_usd) {
  draws <- sort(unique(psa_results$draw))
  treg <- psa_results[psa_results$intervention == "TREG", ]
  treg <- treg[match(draws, treg$draw), ]
  c_treg_non_drug <- treg$total_cost - treg$treg_price  # module header: exact linear algebra

  rows <- lapply(seq_along(draws), function(i) {
    comp_summaries <- stats::setNames(
      lapply(COMPARATOR_THERAPIES, function(tx) {
        r <- psa_results[psa_results$intervention == tx & psa_results$draw == draws[i], ]
        list(qalys = r$qalys, total_cost = r$total_cost)
      }),
      COMPARATOR_THERAPIES
    )
    best <- best_comparator_nmb(comp_summaries, wtp_usd)
    comp <- comp_summaries[[best$comparator]]
    p_star <- (comp$total_cost + wtp_usd * (treg$qalys[i] - comp$qalys) - c_treg_non_drug[i]) / D_BAR
    data.frame(draw = draws[i], p_star = p_star, comparator = best$comparator)
  })
  do.call(rbind, rows)
}

#' Summary of the P*_k density (analysis_plan.md sec 9.1): median with a 95% credible interval,
#' plus the full vector of draws for anyone who wants the density itself (e.g. for a histogram/
#' violin plot -- not reduced to summary statistics only).
ejp_probabilistic <- function(psa_results, wtp_usd) {
  draws <- ejp_probabilistic_draws(psa_results, wtp_usd)
  ci <- stats::quantile(draws$p_star, c(0.025, 0.5, 0.975))
  list(median = unname(ci[2]), ci_low = unname(ci[1]), ci_high = unname(ci[3]),
       draws = draws, wtp_usd = wtp_usd)
}

#' P_50 (analysis_plan.md sec 9.1): the price at which the fraction of draws still cost-effective
#' equals 50%, built from the empirical survival function of `p_star_draws` directly (module
#' header) -- NOT from median(), even though the two are expected to closely agree, because the
#' plan explicitly treats them as distinct quantities to be reported and compared, not one
#' derived from the other. The survival function P(p_star >= P) is a right-continuous step
#' function (module header: INMB(P) is strictly decreasing in P for every draw, so it steps down
#' by exactly 1/n at each sorted p_star value) -- P_50 is the smallest observed p_star at which
#' that step function has dropped to <= 0.5.
ejp_p50 <- function(p_star_draws) {
  sorted <- sort(p_star_draws)
  n <- length(sorted)
  survival <- rev(seq_len(n)) / n  # survival[i] = P(p_star >= sorted[i])
  idx <- which(survival <= 0.5)[1]
  stopifnot(!is.na(idx))
  sorted[idx]
}
