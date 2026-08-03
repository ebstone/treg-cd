# Model Structure

Status: not yet written. This will become the technical specification of the decision-tree + Markov + mixture-cure engine as it is actually implemented in `R/01_decision_tree.R` through `R/04_costs_utilities.R` — health states, transition structure, the week-56 cure-split rule, and cost/utility attachment — kept in sync with the code rather than with the design rationale.

Until then, the working specification is `docs/analysis_plan.md` Section 6 (Model structure), including the state-transition diagram and the mixture-cure extension (Sustained Deep Remission state, week-56 landmark split, relapse hazard). Populate this file once the engine in `R/` is built, so it reflects what the code does rather than what was proposed.
