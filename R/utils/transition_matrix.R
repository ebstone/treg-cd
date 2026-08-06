# Shared transition-matrix helpers used by both R/00_derive_transition_probs.R (loading and
# validating Aliyev et al. 2019's published matrices) and R/02_markov_engine.R (simulating a
# cohort through them). Kept here rather than duplicated, or sourced from R/00, so R/02 doesn't
# have to run R/00's whole load/validate/write-CSV pipeline as a side effect just to get a
# matrix-construction helper.

#' Check that a long-format transition table (from_state, to_state, probability, grouped by
#' `by`) has, for every group and from_state, probabilities in [0,1] summing to 1.
validate_row_sums <- function(df, by, tol = 1e-6) {
  key <- do.call(paste, c(df[c(by, "from_state")], sep = ""))
  sums <- tapply(df$probability, key, sum)
  bad <- sums[abs(sums - 1) > tol]
  if (length(bad) > 0) {
    stop("Row sums not equal to 1 for: ", paste(names(bad), "=", round(bad, 6), collapse = "; "))
  }
  if (any(df$probability < 0 | df$probability > 1)) {
    stop("Probabilities outside [0,1] found in transition table")
  }
  invisible(TRUE)
}

#' Build a square transition matrix (states x states) from a long-format
#' (from_state, to_state, probability) data frame for one therapy. Any from_state present in
#' `states` but absent from the data (biologic arms have no "Moderate-Severe" row: patients who
#' deteriorate to Moderate-Severe exit to the CT track rather than continuing on this matrix, by
#' this study's own model design, analysis_plan.md §6.1) is padded as an absorbing self-loop --
#' the correct convention for computing "this matrix's own behaviour in isolation"; the real
#' cross-matrix CT-switch is R/02_markov_engine.R's job, not this helper's.
#'
#' Rows are then normalised to sum to exactly 1. Aliyev's published values are rounded to 3-4
#' significant figures, so raw row sums are off by up to ~0.1% (e.g. Table 4's UST Surgery row
#' sums to 1.00096, not 1) -- validate_row_sums()'s tolerance already accepts this as expected
#' publication rounding, not a bug. Left uncorrected, that ~0.1% compounds *multiplicatively*
#' every cycle a Markov engine multiplies through this matrix: over a 60-cycle run it inflates
#' total cohort size by over 1%, and over a lifetime horizon (hundreds of cycles) by several
#' percent -- large enough to matter for a research result. Normalising here fixes it at the
#' one place it can be fixed without touching the verbatim-transcribed source CSVs (which stay
#' as published, rounding warts included, for provenance).
build_transition_matrix <- function(df, states) {
  m <- matrix(0, nrow = length(states), ncol = length(states), dimnames = list(states, states))
  for (i in seq_len(nrow(df))) {
    m[df$from_state[i], df$to_state[i]] <- df$probability[i]
  }
  missing_from <- states[rowSums(m) == 0]
  for (s in missing_from) m[s, s] <- 1
  sweep(m, 1, rowSums(m), "/")
}

#' Convert a square transition matrix to its n-cycle equivalent (Chapman-Kolmogorov). Not used
#' in the default pipeline as of the 2026-08-04 switch to Aliyev's native 2-week cycle (nothing
#' needs cycle-length conversion any more), but kept as a general, tested utility -- e.g. for a
#' future N-cycle-ahead projection without needing the full per-cycle trace.
matrix_power <- function(m, n) {
  stopifnot(n >= 1, nrow(m) == ncol(m))
  result <- m
  if (n > 1) for (i in seq_len(n - 1)) result <- result %*% m
  result
}

#' Convert a single cumulative event probability from one cycle length to another under the
#' constant-hazard (exponential-survival) assumption -- the same DEALE-style conversion Aliyev's
#' own Appendix S2 worked example uses (`R/00_derive_transition_probs.R`'s module header: UST
#' week-6 remission 0.349 -> rate = -ln(1-0.349)/(6/52) -> 2-week probability
#' 1-exp(-rate*2/52) = 0.133), generalised here to an arbitrary FROM/TO cycle length rather than
#' hardcoded to a 2-week target. `p_new = 1 - (1-p_old)^(weeks_new/weeks_old)` is algebraically the
#' same two-step (probability -> rate -> probability) computation collapsed into one line, without
#' the intermediate rate ever needing to be materialised (and without an intermediate rounding
#' step a two-line version would introduce).
#'
#' This is a single-event conversion, not a full matrix operation: it's only exact when `p_old` is
#' the probability of ONE event (e.g. "leaves state X by any route"), not when applied
#' independently to several competing destination probabilities from the same row (which would
#' need a proper competing-risks or matrix-fractional-power treatment this project has no existing
#' utility for -- `matrix_power()` above only does forward INTEGER powers). Used by
#' `R/utils/surgery_row_sensitivity.R` (the R2 sensitivity check, `docs/model_audit_v6.md` A17's
#' own follow-up note) exactly this way: converts the aggregate "leaves Surgery" probability, then
#' re-splits the converted mass across destinations in their ORIGINAL relative proportions, rather
#' than converting each destination independently.
convert_cycle_length_probability <- function(p_old, weeks_old, weeks_new) {
  stopifnot(p_old >= 0, p_old <= 1, weeks_old > 0, weeks_new > 0)
  1 - (1 - p_old) ^ (weeks_new / weeks_old)
}

#' Replace every non-Death row's Death-column probability with `death_prob`, uniformly (the "no
#' CD excess mortality" assumption, analysis_plan.md §7.1 item 7 -- background mortality applies
#' identically regardless of current CD health state, so one scalar per cycle is the right
#' granularity, not one per state), rescaling that row's remaining entries proportionally so it
#' still sums to exactly 1. The Death row itself (an absorbing self-loop, probability 1) is left
#' untouched -- rescaling "probability of death, given already dead" isn't a meaningful
#' operation. Used by R/02_markov_engine.R's lifetime-horizon path (R/utils/life_table.R) to
#' substitute Aliyev's own small trial-cohort mortality figure for a life-table one at each
#' cycle's attained age -- see that module's header for why REPLACE, not ADD.
#'
#' A row whose original Death probability was already 1 (a padded absorbing self-loop from
#' build_transition_matrix(), e.g. the biologic-arm Moderate-Severe row -- see that function's
#' own header) is left alone: an absorbing row has nowhere else to redistribute mass, and the
#' rescaling below would divide by a zero denominator. A row with NO Death entry at all (the
#' same padded self-loops, `old_death == 0`) is handled by the general formula correctly --
#' `scale = (1 - death_prob) / (1 - 0)` shrinks the self-loop and installs `death_prob` in its
#' place, which is exactly right: even a state this matrix treats as absorbing in isolation
#' (M-S under a biologic arm, before R/02's own CT-switch mechanic sweeps it away) is not
#' actually immune to background mortality for the one step it's evaluated.
age_adjust_matrix <- function(m, death_prob) {
  stopifnot(death_prob >= 0, death_prob <= 1, "Death" %in% colnames(m))
  states <- rownames(m)
  death_idx <- match("Death", states)
  out <- m
  for (i in seq_len(nrow(m))) {
    if (i == death_idx) next
    old_death <- m[i, death_idx]
    if (old_death >= 1) next  # already fully absorbing; nothing left to rescale
    scale <- (1 - death_prob) / (1 - old_death)
    out[i, ] <- m[i, ] * scale
    out[i, death_idx] <- death_prob
  }
  out
}

#' Simulate a cohort through `n_cycles` of a single transition matrix, returning the full
#' per-cycle occupancy trace (not just the endpoint). Moved here from R/02_markov_engine.R on
#' 2026-08-06 (B1 fix, see R/00_derive_transition_probs.R's own history comment): R/00 needs this
#' to correctly re-derive the published induction table's terminal end-of-induction split (running
#' Aliyev's own 2-week induction matrix forward `induction_cycles` times), and R/00 is sourced
#' before R/02 in this project's dependency order -- putting the primitive here, alongside the
#' other matrix-construction/validation helpers both R/00 and R/02 already share, avoids R/00
#' having to source "downstream" of itself. R/02's own maintenance engine still uses this
#' identically (unchanged behaviour, just relocated) via its own top-of-file source() of this file.
#'
#' Row 1 of the returned matrix is `initial_state` (cycle 0, before any transition); row i+1 is
#' the occupancy after i cycles. Cohort conservation (each row summing to whatever
#' `initial_state` summed to) follows directly from `m`'s rows summing to 1 -- validated at load
#' time by validate_row_sums(), not re-checked here on every multiply.
simulate_cohort <- function(m, initial_state, n_cycles) {
  states <- rownames(m)
  stopifnot(!is.null(states), length(initial_state) == length(states), n_cycles >= 0)
  trace <- matrix(0, nrow = n_cycles + 1, ncol = length(states), dimnames = list(NULL, states))
  trace[1, ] <- initial_state
  for (t in seq_len(n_cycles)) {
    trace[t + 1, ] <- trace[t, ] %*% m
  }
  trace
}

#' Long-format (from_state, to_state, probability) rows from a square matrix.
matrix_to_long <- function(m, therapy) {
  states <- rownames(m)
  do.call(rbind, lapply(states, function(from) {
    data.frame(therapy = therapy, from_state = from, to_state = states,
               probability = unname(m[from, ]), stringsAsFactors = FALSE)
  }))
}
