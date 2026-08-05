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

#' Long-format (from_state, to_state, probability) rows from a square matrix.
matrix_to_long <- function(m, therapy) {
  states <- rownames(m)
  do.call(rbind, lapply(states, function(from) {
    data.frame(therapy = therapy, from_state = from, to_state = states,
               probability = unname(m[from, ]), stringsAsFactors = FALSE)
  }))
}
