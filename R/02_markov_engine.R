# Cohort Markov maintenance engine (analysis_plan.md §6.1, §12.1): custom matrix-based engine
# (list of transition matrices multiplied through a state-occupancy vector), not heemod, since
# the week-56 cure landmark and arm-specific CT-switch rules are awkward in heemod's formula
# interface. Reproduce the v6 base case as a regression test before adding the cure structure.
