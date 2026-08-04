# Induction decision tree (analysis_plan.md §6.1): partitions each arm's cohort at end of
# induction into Remission, Mild, M-SR, or non-response (-> CT Markov track).
#
# The state-transition diagram (§6.3) makes a detail explicit that R/02_markov_engine.R's
# run_maintenance_arm() alone doesn't handle: induction non-responders don't enter the
# biologic's maintenance matrix at all -- they go straight to the CT track. So this script's
# job is to turn one row of R/00's load_published_induction() output into the *two* initial
# occupancy vectors run_maintenance_arm() needs: what starts on the biologic (everyone who
# responded), and what starts already on CT (everyone who didn't).

if (!exists("MAINTENANCE_STATES")) source("R/utils/transition_matrix.R")
if (!exists("load_published_induction")) source("R/00_derive_transition_probs.R")

#' Split one therapy's induction row (from load_published_induction()) into the initial
#' occupancy vectors for run_maintenance_arm()'s two tracks.
#'
#' Both vectors are normalised by their combined sum: Aliyev's published induction proportions
#' are rounded to 3-4 significant figures and don't sum to exactly 1 (e.g. UST's row sums to
#' 0.9994, not 1 -- see load_published_induction()'s own tolerance note in R/00). Left
#' unnormalised, every downstream aggregate (QALYs, costs, population counts) would carry that
#' same fraction-of-a-percent error for no reason -- the same rationale build_transition_matrix()
#' already applies to the maintenance matrices themselves.
run_decision_tree <- function(induction_row) {
  stopifnot(nrow(induction_row) == 1)

  on_biologic <- stats::setNames(rep(0, length(MAINTENANCE_STATES)), MAINTENANCE_STATES)
  on_biologic["Moderate-Severe Responder"] <- induction_row$to_moderate_severe_responder
  on_biologic["Mild"] <- induction_row$to_mild
  on_biologic["Remission"] <- induction_row$to_remission
  # Moderate-Severe (non-response), Surgery, Death all stay 0: non-responders never start on
  # the biologic, and induction has no route to Surgery/Death.

  on_ct <- stats::setNames(rep(0, length(MAINTENANCE_STATES)), MAINTENANCE_STATES)
  on_ct["Moderate-Severe"] <- induction_row$to_moderate_severe

  total <- sum(on_biologic) + sum(on_ct)
  list(initial_on_biologic = on_biologic / total, initial_on_ct = on_ct / total)
}

#' run_decision_tree() for every therapy in load_published_induction()'s output (UST/IFX/ADA --
#' Table 3 has no CT row, since CT isn't an induction arm in its own right; TREG isn't in Table 3
#' either, since its efficacy is a modelling assumption rather than a sourced trial endpoint,
#' see analysis_plan.md §15 Decision 4).
run_decision_tree_all <- function(raw_dir = "data/raw") {
  induction <- load_published_induction(raw_dir)
  therapies <- induction$therapy
  stats::setNames(
    lapply(therapies, function(tx) run_decision_tree(induction[induction$therapy == tx, ])),
    therapies
  )
}
