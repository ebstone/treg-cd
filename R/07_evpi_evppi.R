# EVPI/EVPPI (analysis_plan.md §9.2, §9.3, Aims 2-3): population and per-patient EVPI as a
# surface over price x WTP threshold; EVPPI by parameter subset (A-G) via the `voi` package
# (GAM for low-dimensional subsets, SPDE-INLA for higher-dimensional blocks), cross-checked
# against `BCEA::evppi()`. Report a convergence check across draw counts (Gate 4, §14).
#
# FIRST PASS, 2026-08-05. Computes EVPI and EVPPI for exactly the parameters R/06_psa.R actually
# samples -- pi (subset A), Treg's acquisition price (subset C), and the utility chain (subset E)
# -- plus their unions. Subsets B/D/F/G are not attempted: EVPPI for a parameter the PSA holds
# fixed is undefined by construction (A14, docs/model_audit_v6.md, says exactly this about the
# workbook's own PSA), not a choice this module is deferring on its own initiative.
#
# SECOND PASS, 2026-08-06 (peer review 2026-08-05's B3; docs/model_audit_v6.md A19). **Subset B
# (the post-cure relapse hazard h) and the durability block A u B now exist**, because R/06_psa.R
# now samples h instead of hardcoding it to 0. This is the whole point of that fix rather than a
# side effect of it: analysis_plan.md's Aim 3 is stated as a falsifiable comparison -- "the
# durability block (cure fraction + relapse hazard) does or does not carry higher EVPPI than the
# cost block" -- and with h unsampled, the durability block was literally just subset A wearing a
# different label, so Aim 3 could not be answered at all, only half-answered without saying so.
# Every EVPI/EVPPI number this module produced before this pass is superseded, since the draws
# underneath them all assumed a permanent cure.
#
# ---- A real environment gap, not a design choice: voi/BCEA/INLA are unavailable here -----------
#
# analysis_plan.md §9.3 specifies the `voi` package (GAM for low-dimensional subsets, SPDE-INLA
# for higher-dimensional ones) with a `BCEA::evppi()` cross-check. All three attempted and
# recorded here: `voi` and `BCEA` fail to install in this sandbox because their own dependencies
# (`earth`, `mvtnorm`) need a Fortran compiler (`gfortran`) that isn't present, and `INLA` isn't
# even on CRAN. This is an environment limitation, not a judgement that the specified method is
# wrong -- flagged so a future session with a fuller toolchain can actually run the cross-check
# analysis_plan.md asks for, rather than assuming this module already did.
#
# In its place: evppi_gam() below implements the SAME underlying method voi's own low-dimensional
# path uses -- Strong, Oakley & Brennan (2014)'s nonparametric regression estimator -- directly
# against `mgcv::gam()`, which IS available (`mgcv` ships as part of every standard R
# installation, unlike voi/BCEA/INLA, so this has no real dependency risk). Every subset actually
# computable from R/06's current PSA (A, B, C, E, and their unions) is low-dimensional enough (<=7
# raw parameters, PCA-reduced to 2 components above 2 -- see pca_reduce()) that GAM is the
# right tool regardless of voi's own availability; SPDE-INLA's higher-dimensional case (subsets
# D/F, comparator/induction transition probabilities) doesn't arise yet anyway, since R/06 doesn't
# sample those parameters at all (module header, R/06_psa.R). If/when this environment gets a
# working Fortran toolchain, cross_check_voi() below is written and ready -- it just isn't
# exercised by anything else in this module. Note that analysis_plan.md §9.3's own method
# specification names "A, B, C, and the A u B pair" as exactly the low-dimensional GAM cases, so
# adding B and A u B here is running the plan's specified method on the plan's specified subsets,
# not extending it.
#
# ---- Two distinct EVPI framings, both requested, both implemented differently -------------------
#
# analysis_plan.md asks for two different things that both use the word "EVPI":
#
# 1. **section 9.2's price-conditioned EVPI surface**: "EVPI here is conditional on the Treg
#    price, because the decision set itself changes with price... report EVPI as a surface over
#    price x lambda." This is the value of resolving REMAINING uncertainty (pi, utilities) to a
#    decision-maker who has already fixed a specific candidate price P -- price itself is treated
#    as KNOWN, not uncertain, at each point on the surface. evpi_surface() implements this by
#    overriding every draw's Treg cost to a fixed price (treg_price_dependent_dose_cost() is
#    exactly linear in price with slope 1, so this is just algebra: cost_at(P) = cost_at(sampled
#    price) - sampled_price + P -- no need to re-run the Markov engine per price point).
# 2. **section 9.3's subset-C EVPPI**: "is further manufacturing cost-engineering worth funding?"
#    -- this treats price as ONE OF the genuinely uncertain parameters (using its actual sampled
#    Gamma distribution, no override) and asks how much of the TOTAL reference-case decision
#    uncertainty (pi + price + utilities all varying together, exactly as R/06 ran it) is
#    attributable to not knowing the eventual price specifically. evppi_by_subset() implements
#    this using psa_results exactly as run, no price override.
#
# These are not in tension -- they answer different questions ("what if we fix the price" vs.
# "how much does not-yet-knowing-the-price cost us") -- but both are worth stating explicitly so
# a reader doesn't assume price is being treated the same way in both outputs.
#
# ---- pi prior-sensitivity on EVPPI (added 2026-08-05, docs/treg-cd_decision_resolutions_2026-08-05.md sec 3.2) --
#
# U(0,1) (Decision 4's own recorded fallback for pi) maximises prior variance on pi, which makes
# pi come out as the dominant EVPPI parameter "more or less by construction" -- a competent
# reviewer would read that as an artefact of the prior, not a finding, unless it's checked.
# prior_reweight() + evppi_prior_sensitivity() below re-estimate evppi_by_subset() under two
# alternative priors for pi (Beta(1,3): mass toward low cure fractions, consistent with the
# Ovasave/CATS1 experience; Beta(2,2): symmetric, informative) via importance-weighting the SAME
# 10,000 PSA draws R/06_psa.R already produced -- no re-simulation, exactly as the memo specifies.
# The stated finding to report is the RANKING of subsets by EVPPI (expected to be stable across
# all three priors) vs. the LEVEL (expected to shift with the prior) -- evppi_prior_sensitivity()
# computes both; `rank_within_prior` makes the ranking-stability claim directly checkable from the
# returned table rather than requiring a second pass over it.
#
# ---- Not yet implemented in this pass -------------------------------------------------------------
#
# - **Population EVPI** (analysis_plan.md sec 9.4, Decision 6). population_evpi() implements the
#   FORMULA (per-patient EVPI x effective population size) but takes effective_population as a
#   required argument with no default -- Decision 6's own status table marks the moderate-to-
#   severe-eligible and biologic-naive-vs-refractory fractions "Open -- must not be assumed", so
#   there is no sourced number to default to. Not a gap in this module; a gap in Gate 3's
#   parameterisation work.
# - **EVPPI for subsets D/F/G** (comparator maintenance transitions; induction response
#   probabilities; surgery transitions/cost) -- undefined until R/06_psa.R actually samples those
#   parameters (module header, R/06_psa.R, explains why each isn't sampled yet). ~~Subset B~~ is no
#   longer on this list as of 2026-08-06: h IS sampled now. The distinction that kept it here is
#   still worth keeping in mind when reading the two side by side, though --
#   R/05_deterministic_results.R's headroom_frontier_by_duration() is a DETERMINISTIC sweep over T
#   (it shows how pi* moves as durability is varied by hand), while subset B's EVPPI is a
#   probabilistic quantity built on R/06's Gamma prior over h. They answer different questions and
#   the sweep is NOT a substitute for the EVPPI, which is why B stayed genuinely undefined until
#   the PSA sampled h rather than being "already covered by the sweep."
# - **The voi/BCEA cross-check** analysis_plan.md sec 9.3 asks for -- see above; written
#   (cross_check_voi()) but inert in this environment.

if (!exists("MAINTENANCE_STATES")) source("R/utils/transition_matrix.R")
if (!exists("run_psa")) source("R/06_psa.R")

# ---- Core EVPI arithmetic -------------------------------------------------------------------

#' Per-patient EVPI from a draws x arms net-benefit matrix (analysis_plan.md sec 9.2's
#' definition): E_theta[max_j NB_j(theta)] - max_j E_theta[NB_j(theta)] -- the mean benefit of
#' the decision made with perfect foresight, minus the benefit of the single best decision made
#' under current information. Always >= 0 (Jensen's inequality on max()); 0 exactly when one arm
#' dominates in every draw, i.e. there is no decision uncertainty left to resolve.
#'
#' `weights = NULL` (default): every draw counted equally, i.e. the expectation is taken under
#' whatever prior actually generated the draws (R/06_psa.R: pi ~ U(0,1), Decision 4's fallback) --
#' exactly the prior arithmetic this function always had, byte-identical to before `weights`
#' existed. A non-NULL vector (same length as `nrow(nb)`, needn't already sum to 1 -- renormalised
#' internally) re-weights the SAME draws to approximate the expectation under a DIFFERENT prior
#' via importance sampling (prior_reweight() below computes exactly such a vector) -- no new PSA
#' draws, no re-simulation, per `docs/treg-cd_decision_resolutions_2026-08-05.md` §3.2.
evpi_from_nb <- function(nb, weights = NULL) {
  nb <- as.matrix(nb)
  if (is.null(weights)) weights <- rep(1, nrow(nb))
  stopifnot(length(weights) == nrow(nb), all(weights >= 0))
  w <- weights / sum(weights)
  sum(apply(nb, 1, max) * w) - max(colSums(nb * w))
}

#' Wide draws x arms net-benefit matrix at a given WTP threshold, built from R/06_psa.R's
#' long-format run_psa() output. `treg_price_override`, if supplied, replaces every draw's
#' sampled Treg price with this fixed value (module header, framing 1) by exploiting
#' treg_price_dependent_dose_cost()'s exact linearity in price (R/05_deterministic_results.R):
#' adjusted_cost = original_cost - original_sampled_price + override. Left NULL (default), the
#' natural sampled price is used as-is (module header, framing 2).
net_benefit_matrix <- function(psa_results, wtp_usd, treg_price_override = NULL) {
  arms <- unique(psa_results$intervention)
  draws <- sort(unique(psa_results$draw))

  cost <- psa_results$total_cost
  if (!is.null(treg_price_override)) {
    is_treg <- psa_results$intervention == "TREG"
    cost[is_treg] <- cost[is_treg] - psa_results$treg_price[is_treg] + treg_price_override
  }
  nmb <- wtp_usd * psa_results$qalys - cost

  wide <- matrix(NA_real_, nrow = length(draws), ncol = length(arms), dimnames = list(NULL, arms))
  for (a in arms) {
    rows <- psa_results$intervention == a
    wide[match(psa_results$draw[rows], draws), a] <- nmb[rows]
  }
  stopifnot(!anyNA(wide))
  wide
}

#' The price x WTP EVPI surface analysis_plan.md sec 9.2 asks for (framing 1, module header):
#' one row per (price, wtp) grid point, holding Treg's price fixed at that price for every draw
#' and computing per-patient EVPI over the remaining sampled uncertainty (pi, utilities).
evpi_surface <- function(psa_results, price_grid_usd, wtp_grid_usd) {
  grid <- expand.grid(price_usd = price_grid_usd, wtp_usd = wtp_grid_usd)
  grid$evpi_per_patient <- mapply(function(price, wtp) {
    evpi_from_nb(net_benefit_matrix(psa_results, wtp, treg_price_override = price))
  }, grid$price_usd, grid$wtp_usd)
  grid
}

#' Population EVPI = per-patient EVPI x effective population size (analysis_plan.md sec 9.4).
#' `effective_population` has NO default -- see module header: the fractions Decision 6 needs to
#' compute it (moderate-to-severe-eligible, biologic-naive-vs-refractory) are not sourced anywhere
#' in this repository yet, and inventing one would be exactly the fabrication this project's
#' methodology exists to avoid. Exists so the arithmetic is ready the moment Gate 3 supplies a
#' real number, not to encourage supplying a placeholder now.
population_evpi <- function(per_patient_evpi, effective_population) {
  stopifnot(effective_population >= 0)
  per_patient_evpi * effective_population
}

# ---- Dimension reduction for higher-dimensional subsets (analysis_plan.md sec 9.3) --------------

#' First `n_components` principal-component scores of `draws_matrix` (analysis_plan.md sec 9.3:
#' "where subsets are dimension-heavy, consider a principal-components pre-reduction and report
#' it"). Used for the utility-chain subset (5 raw sampled values) before handing them to
#' evppi_gam() -- a 5-dimensional tensor smooth is both statistically unstable at typical PSA
#' sample sizes and slow to fit; 2 components capture the great majority of the chain's variation
#' since all 5 values are driven by only 5 independent Uniform draws combined multiplicatively,
#' not 5 truly free dimensions of variation.
pca_reduce <- function(draws_matrix, n_components = 2) {
  pr <- stats::prcomp(draws_matrix, center = TRUE, scale. = TRUE)
  n_components <- min(n_components, ncol(pr$x))
  pr$x[, seq_len(n_components), drop = FALSE]
}

# ---- GAM-based EVPPI (Strong, Oakley & Brennan 2014; voi's own low-dimensional method) -----------

#' EVPPI for one parameter subset, estimated by nonparametric regression (module header): fit
#' each arm's net benefit as a smooth function of `param_matrix`'s columns via `mgcv::gam()`,
#' then EVPPI = mean(rowMax(fitted)) - max(colMeans(actual nb)) -- the fitted values estimate
#' E[NB_j | parameter subset], and swapping them in for the raw draws in evpi_from_nb()'s formula
#' is exactly the Strong-Oakley-Brennan construction. A single parameter uses a thin-plate
#' smooth (`s()`); 2+ use a tensor product smooth (`te()`) so genuine interactions between the
#' parameters are captured, not just their marginal effects.
#'
#' `te()`'s default per-margin basis size means its TOTAL basis grows roughly as k^d -- harmless
#' at 1-2 dimensions (mgcv's own defaults there are already well-tuned; an explicit override at
#' d<=2 was tried and measured SLOWER than leaving it alone -- 2D: ~1s at mgcv's default vs ~20s
#' forcing k=14 -- REML's cost isn't linear in basis size, so a bigger explicit k is a pure loss
#' when the default already works) but a real performance cliff at 3+ (measured on 2,000 draws,
#' this machine: 3D ~12s; 4D still running after 5 minutes). `k_per_margin` only overrides the
#' default once d >= 3, as a defensive backstop for a direct caller who skips pre-reduction --
#' evppi_by_subset() below always pre-reduces via reduce_for_gam() first, so in practice this
#' module never actually exercises the d>=3 branch itself; it exists for evppi_gam() to still be
#' safe to call directly with a wider parameter matrix.
#'
#' `weights = NULL` (default): unweighted `mgcv::gam()` fit and a plain-mean EVPPI formula,
#' byte-identical to before `weights` existed. A non-NULL vector is passed straight through as
#' `mgcv::gam()`'s own `weights` argument (its standard prior-weights mechanism, the same one
#' `glm()` uses) AND into the final EVPPI arithmetic via `evpi_from_nb()`'s own weights argument --
#' both stages of the estimator need to agree on the prior, not just the regression fit, or the
#' fitted E[NB|subset] curve and the E[max(...)] it's compared against would be computed under two
#' different priors. Existing solely so evppi_by_subset() can implement
#' `docs/treg-cd_decision_resolutions_2026-08-05.md` §3.2's prior-sensitivity analysis by
#' re-weighting the SAME draws under alternative priors for pi, with no re-simulation.
evppi_gam <- function(nb, param_matrix, max_basis_size = 200, weights = NULL) {
  param_matrix <- as.matrix(param_matrix)
  stopifnot(nrow(param_matrix) == nrow(nb))
  d <- ncol(param_matrix)
  colnames(param_matrix) <- paste0("p", seq_len(d))
  df <- data.frame(param_matrix)

  if (d == 1) {
    smooth_term <- sprintf("s(%s)", colnames(param_matrix)[1])
  } else if (d == 2) {
    smooth_term <- sprintf("te(%s)", paste(colnames(param_matrix), collapse = ","))
  } else {
    k_per_margin <- max(3, floor(max_basis_size^(1 / d)))
    smooth_term <- sprintf("te(%s, k = %d)", paste(colnames(param_matrix), collapse = ","), k_per_margin)
  }

  fitted <- vapply(colnames(nb), function(arm) {
    df$nb_arm <- nb[, arm]
    fit <- mgcv::gam(stats::as.formula(paste("nb_arm ~", smooth_term)), data = df, method = "REML",
                      weights = weights)
    as.numeric(stats::predict(fit))
  }, numeric(nrow(nb)))

  w <- if (is.null(weights)) rep(1, nrow(nb)) else weights
  w <- w / sum(w)
  sum(apply(fitted, 1, max) * w) - max(colSums(nb * w))
}

#' PCA-reduce `m` to `n_components` columns only if it has more than `max_raw_dims` columns to
#' begin with -- e.g. leaves a 1- or 2-column subset (A, C, A u C) untouched, exactly as
#' evppi_gam()'s own `te()` handles those cheaply and exactly, and only reduces (and reports the
#' variance retained) once a subset actually needs it (E's 5 raw columns, and any union
#' involving E). This is what keeps evppi_by_subset() itself from hitting the te()-dimension
#' blow-up evppi_gam()'s own `k_per_margin` cap is only a backstop for (module header).
reduce_for_gam <- function(m, max_raw_dims = 2, n_components = 2) {
  m <- as.matrix(m)
  n_raw <- ncol(m)
  if (n_raw <= max_raw_dims) {
    return(list(matrix = m, n_raw = n_raw, n_used = n_raw, variance_explained = 1))
  }
  pr <- stats::prcomp(m, center = TRUE, scale. = TRUE)
  n_used <- min(n_components, ncol(pr$x))
  variance_explained <- sum((pr$sdev[seq_len(n_used)])^2) / sum(pr$sdev^2)
  list(matrix = pr$x[, seq_len(n_used), drop = FALSE], n_raw = n_raw, n_used = n_used,
       variance_explained = variance_explained)
}

#' EVPPI for every subset actually computable from R/06_psa.R's current PSA output (module
#' header): A (pi), B (post-cure relapse hazard h), C (Treg price), E (the utility chain), and
#' their decision-relevant unions -- each evaluated at the REFERENCE scenario (no price override,
#' module header framing 2), against the same total EVPI so every row's share of it is directly
#' comparable.
#'
#' **`A u B` is the row Aim 3 is actually about**, and it is in analysis_plan.md §9.3's own subset
#' table under that name ("Durability block -- combined value of the whole durability question;
#' report jointly, A and B are decision-relevant together"). Aim 3's falsifiable claim is the
#' comparison of that row against C (the cost block): a single trial extension resolving both the
#' cure fraction and its durability is one fundable research action, so its value is EVPPI(A u B),
#' not EVPPI(A) + EVPPI(B) -- EVPPI is not additive over subsets, and reporting the sum in its
#' place would be a methodological error, not a shortcut. `A u B u C u E` (added 2026-08-06,
#' replacing the old `A u C u E`) is the all-sampled-parameters row: it is a useful upper reference
#' against total EVPI, and it necessarily moved when h joined the sampled set, since "all
#' parameters" now means four blocks rather than three.
#'
#' Any subset with more than 2 raw
#' parameters (E alone, and every union involving it) is PCA-reduced to 2 components first
#' (reduce_for_gam(), analysis_plan.md sec 9.3's own "consider a principal-components pre-
#' reduction and report it") -- `n_raw_params`/`n_params_used`/`variance_explained` in the
#' returned table IS that reporting.
#'
#' **Known artifact of the PCA shortcut, found while validating this function:** true EVPPI is
#' monotone in the subset (resolving MORE parameters can only raise or hold the value of perfect
#' information about them, never lower it) -- but a PCA-reduced superset can score BELOW an
#' unreduced subset it contains, because cramming more raw parameters into the same 2 components
#' dilutes how much of any one of them survives the reduction. Observed directly: `A u C` (2 raw
#' dims, no reduction needed) can score higher than `A u C u E` (7 raw dims reduced to the same 2
#' components) even though the latter is a strict superset. This is a property of the
#' approximation, not a bug in the arithmetic -- don't read the subset table as a clean monotone
#' ladder; `variance_explained` tells you how much to trust a given reduced row's number.
#'
#' Parameter values are pulled from the TREG rows only: pi_sdr, treg_price and
#' relapse_hazard_annual are all NA for the three comparators (R/06_psa.R's run_psa(), by
#' construction -- they don't depend on any of them), while the 5 util_* columns are identical
#' across all 4 arms within a draw, so reading them from TREG's rows alone is not a loss of
#' information, just a convenient single source with no NAs. That the relapse hazard is a
#' Treg-only column is a modelling fact, not a data-layout convenience: h enters only through the
#' SDR track (R/03_cure_fraction_module.R), which no comparator arm has.
#'
#' `weights = NULL` (default): every draw counted equally under the prior that actually generated
#' it (U(0,1) for pi, module header) -- unchanged from before `weights` existed. A non-NULL vector
#' (one weight per draw, in ascending-draw-index order -- prior_reweight()'s own output, or
#' evppi_prior_sensitivity()'s caller of it, already produces weights in this order) is threaded
#' into both `evpi_from_nb()` (the total-EVPI denominator) and every `evppi_gam()` call (the
#' subset numerators), so the whole table -- total EVPI, each subset's EVPPI, and each
#' `proportion_of_total_evpi` -- is internally consistent under the SAME re-weighted prior, not a
#' mix of the reference prior in one place and the alternative in another.
evppi_by_subset <- function(psa_results, wtp_usd, weights = NULL) {
  nb <- net_benefit_matrix(psa_results, wtp_usd)
  total_evpi <- evpi_from_nb(nb, weights)

  treg <- psa_results[psa_results$intervention == "TREG", ]
  treg <- treg[order(treg$draw), ]
  util_cols <- c("util_modsev", "util_resp", "util_mild", "util_remission", "util_surgery")

  # Fail loudly, here, on a psa_results produced before 2026-08-06: such a table has no
  # relapse_hazard_annual column at all, and cbind()-ing a NULL column would silently drop subset
  # B's parameter rather than erroring, leaving a "B" row in the output table that was fitted
  # against nothing. A stale draws file is exactly the failure mode this project has been bitten by
  # before (A16's stray PSA draw), so it gets an explicit check rather than a downstream NA.
  stopifnot("relapse_hazard_annual" %in% names(treg), !anyNA(treg$relapse_hazard_annual))

  raw_subsets <- list(
    A = cbind(pi = treg$pi_sdr),
    B = cbind(h = treg$relapse_hazard_annual),
    C = cbind(price = treg$treg_price),
    E = as.matrix(treg[, util_cols]),
    `A u B` = cbind(pi = treg$pi_sdr, h = treg$relapse_hazard_annual),
    `A u C` = cbind(pi = treg$pi_sdr, price = treg$treg_price),
    `A u E` = cbind(pi = treg$pi_sdr, as.matrix(treg[, util_cols])),
    `C u E` = cbind(price = treg$treg_price, as.matrix(treg[, util_cols])),
    `A u B u C u E` = cbind(pi = treg$pi_sdr, h = treg$relapse_hazard_annual,
                             price = treg$treg_price, as.matrix(treg[, util_cols]))
  )

  reduced <- lapply(raw_subsets, reduce_for_gam)
  evppi_vals <- vapply(reduced, function(r) evppi_gam(nb, r$matrix, weights = weights), numeric(1))

  data.frame(
    subset = names(reduced), evppi = evppi_vals,
    proportion_of_total_evpi = evppi_vals / total_evpi,
    n_raw_params = vapply(reduced, `[[`, numeric(1), "n_raw"),
    n_params_used = vapply(reduced, `[[`, numeric(1), "n_used"),
    variance_explained = vapply(reduced, `[[`, numeric(1), "variance_explained"),
    total_evpi = total_evpi, wtp_usd = wtp_usd,
    row.names = NULL
  )
}

# ---- pi prior-sensitivity analysis on EVPPI (resolutions memo §3.2) --------------------------

#' Importance weights turning the existing pi ~ U(0,1) PSA draws (R/06_psa.R, Decision 4's own
#' recorded fallback) into a sample representative of a DIFFERENT prior `dprior` (a density
#' function taking a vector of pi values, e.g. `function(pi) stats::dbeta(pi, 1, 3)`), by the
#' standard importance-sampling identity w_k ~ f_new(pi_k) / f_unif(pi_k). f_unif(pi_k) = 1
#' everywhere on [0, 1] (the Uniform density), so the ratio collapses to just f_new(pi_k) itself
#' -- renormalised here to mean 1 so the result can be passed straight into evpi_from_nb()'s or
#' evppi_gam()'s `weights` argument without the caller having to rescale anything else downstream.
#' This is `docs/treg-cd_decision_resolutions_2026-08-05.md` §3.2's whole point: "this is cheap.
#' It requires no re-simulation" -- the SAME 10,000 draws R/06_psa.R already produced are reused,
#' just re-weighted, to approximate what EVPPI would have looked like under a prior the PSA never
#' actually sampled from.
#'
#' **Checked when h joined the sampled set (2026-08-06, B3): this function still needs to know
#' about pi only, and that is a derivation, not an assumption.** The concern is real -- these
#' weights are now applied to a table whose draws vary in four parameters, not three, so it is
#' worth writing down why a pi-only weight is still the correct importance weight. R/06_psa.R
#' draws pi and h independently of each other and of everything else, so the joint sampling density
#' factorises, f_old(pi, h, price, u) = f_unif(pi) f_gamma(h) f(price) f(u). The alternative prior
#' being asked about changes the pi factor and nothing else, so the importance-weight ratio is
#' w = f_new(pi) f_gamma(h) f(price) f(u) / [f_unif(pi) f_gamma(h) f(price) f(u)] = f_new(pi) /
#' f_unif(pi) -- every non-pi factor cancels exactly, whatever it is and however many there are.
#' The weights below are therefore unchanged and still exact. (Had pi and h been sampled with any
#' dependence between them -- a correlated durability prior, say -- this would NOT hold and the
#' weight would need the joint density.)
#'
#' The same identity also says what an h-prior sensitivity analysis would look like, if one is
#' wanted: weights of dgamma(h, new_shape, new_rate) / dgamma(h, RELAPSE_HAZARD_PSA_GAMMA_SHAPE,
#' shape/mean) evaluated on the sampled h column, passed through the identical machinery. Not
#' implemented here -- the shape is a single documented modelling choice
#' (sample_relapse_hazard_draws(), R/06_psa.R) rather than a contested one like pi's U(0,1), and
#' the deterministic T sweep already reports durability sensitivity in units a reader can read --
#' but it is a two-line addition on top of this function if a reviewer asks for it.
prior_reweight <- function(pi_values, dprior) {
  w <- dprior(pi_values)
  stopifnot(all(is.finite(w)), all(w >= 0), any(w > 0))
  w / mean(w)
}

#' The three priors for pi analysis_plan.md's own uniform-prior fallback (Decision 4) needs
#' checked against, per the resolutions memo §3.2: U(0,1) maximises prior variance on pi, so pi
#' coming out as the dominant EVPPI parameter "more or less by construction" under it is exactly
#' the artefact a competent reviewer would flag -- Beta(1,3) (mass toward LOW cure fractions,
#' consistent with the Ovasave/CATS1 experience of a cure that mostly doesn't take) and Beta(2,2)
#' (symmetric, informative, still centred on 0.5) are the memo's own named alternatives. U(0,1)'s
#' entry is included as the reference case (weight 1 for every draw, i.e. no reweighting at all)
#' so evppi_prior_sensitivity() below can report all three on one consistent table without a
#' separate unweighted code path.
PI_PRIOR_SENSITIVITY_SPECS <- list(
  `U(0,1) [reference]` = function(pi) stats::dunif(pi, 0, 1),
  `Beta(1,3)` = function(pi) stats::dbeta(pi, 1, 3),
  `Beta(2,2)` = function(pi) stats::dbeta(pi, 2, 2)
)

#' evppi_by_subset() re-estimated once per prior in `prior_specs` (default:
#' PI_PRIOR_SENSITIVITY_SPECS), via prior_reweight()'s importance weights on the SAME PSA draws --
#' no new simulation. Returns one evppi_by_subset() table per prior, row-bound with a `prior`
#' column, plus `rank_within_prior` (1 = highest EVPPI, computed separately within each prior) so
#' the memo's own stated finding is directly readable off the table: "report the RANKING of EVPPI
#' by parameter subset as the finding, and the LEVEL as prior-dependent" (§3.2) -- a stable rank
#' across all three `prior` groups alongside a shifting `evppi` column IS that finding; this
#' function computes the ingredients, the manuscript states the conclusion.
#'
#' Since 2026-08-06 the table this reweights includes subsets B and A u B, so the memo's own
#' ranking-vs-level caution now covers Aim 3's durability-block-vs-cost-block comparison directly:
#' if A u B outranks C under all three priors for pi, that ranking is the reportable finding and
#' its dollar level is not. Only pi's prior is varied here (prior_reweight()'s docstring derives
#' why pi-only weights remain exact now that h is also sampled).
evppi_prior_sensitivity <- function(psa_results, wtp_usd, prior_specs = PI_PRIOR_SENSITIVITY_SPECS) {
  treg <- psa_results[psa_results$intervention == "TREG", ]
  treg <- treg[order(treg$draw), ]
  pi_values <- treg$pi_sdr

  rows <- lapply(names(prior_specs), function(prior_name) {
    w <- prior_reweight(pi_values, prior_specs[[prior_name]])
    res <- evppi_by_subset(psa_results, wtp_usd, weights = w)
    res$prior <- prior_name
    res
  })
  combined <- do.call(rbind, rows)
  # ties.method = "min": GAM-estimation noise can put two subsets' EVPPI within floating-point
  # distance of each other (observed directly among small-effect subsets), where the default
  # tie-averaging rank() would return non-integer ranks like 4.5 -- "min" gives both tied subsets
  # the better (lower, i.e. more favourable) integer rank instead, which is what "rank 1 = highest
  # EVPPI" should mean for a genuine (if close) tie.
  combined$rank_within_prior <- stats::ave(-combined$evppi, combined$prior,
                                            FUN = function(x) rank(x, ties.method = "min"))
  combined
}

# ---- Convergence check across draw counts (analysis_plan.md sec 9.3's explicit requirement) -----

#' EVPPI estimate vs. number of draws, for one parameter subset (`param_cols`: 1-2 raw column
#' names from `psa_results`, e.g. "pi_sdr" or c("pi_sdr", "treg_price") -- for a higher-dimensional
#' subset, PCA-reduce first via reduce_for_gam()/pca_reduce() and pass the resulting component
#' columns instead, same as evppi_by_subset() does), by re-estimating on successive PREFIXES of
#' `psa_results`' draws (1:n for each n in `draw_counts`) rather than re-running run_psa() at each
#' size. Valid because PSA draws are i.i.d. and unordered by construction -- the first n draws of
#' a large, properly random set are statistically equivalent to a fresh run of n draws -- and far
#' cheaper than re-running the Markov engine from scratch at every grid point. `draw_counts`
#' should not exceed the number of draws actually in `psa_results`.
evppi_convergence <- function(psa_results, param_cols, wtp_usd, draw_counts) {
  max_draw <- max(psa_results$draw)
  stopifnot(all(draw_counts <= max_draw))

  estimates <- vapply(draw_counts, function(n) {
    subset_results <- psa_results[psa_results$draw <= n, ]
    nb <- net_benefit_matrix(subset_results, wtp_usd)
    treg <- subset_results[subset_results$intervention == "TREG", ]
    treg <- treg[order(treg$draw), ]
    evppi_gam(nb, as.matrix(treg[, param_cols, drop = FALSE]))
  }, numeric(1))

  data.frame(n_draws = draw_counts, evppi = estimates)
}

# ---- voi/BCEA cross-check (written, inert in this environment -- module header) -----------------

#' Cross-check evppi_gam()'s estimate against `voi::evppi()` for one subset, IF `voi` is
#' installed (analysis_plan.md sec 9.3: "cross-check at least two subsets with BCEA::evppi() and
#' report agreement"). Returns NULL with a message, rather than an error, when `voi` isn't
#' available -- true in this sandbox (module header), but this function is written so a future
#' session with a working Fortran toolchain can actually run the check without editing this file.
cross_check_voi <- function(psa_results, param_matrix, wtp_usd) {
  if (!requireNamespace("voi", quietly = TRUE)) {
    message("voi is not installed in this environment (module header: gfortran unavailable, ",
            "voi/BCEA/mvtnorm/earth fail to build) -- skipping cross-check.")
    return(invisible(NULL))
  }
  nb <- net_benefit_matrix(psa_results, wtp_usd)
  voi::evppi(outputs = nb, inputs = data.frame(param_matrix), pars = colnames(param_matrix))
}
