# Early Economic Evaluation of Allogeneic Treg Therapy for Crohn's Disease

Cost-effectiveness analysis of a hypothetical allogeneic T-regulatory cell ("Treg") therapy
versus ustekinumab and infliximab (± adalimumab) for moderate-to-severe Crohn's disease.
Prepared by D. Jadambaa, E. Stone and B. Abraham (Johns Hopkins Bloomberg School of Public
Health) for submission to *PharmacoEconomics*.

Because no Treg efficacy data yet exist, the paper's primary outputs are the **economically
justifiable price (EJP)**, the **minimum durable cure fraction required for cost-effectiveness**
(a "headroom" analysis), and **expected value of (partial) perfect information** (EVPI/EVPPI) —
not an ICER built on an invented efficacy assumption. See `docs/analysis_plan.md` for the full
rationale, model design and open decisions.

## Status

Gate 0 closed 2026-08-04 (`docs/analysis_plan.md` §15 — Decisions 1, 4, 5 final; Decisions 2, 3,
6 recorded provisionally pending full co-author sign-off). Gate 1 transition probabilities closed
2026-08-04 — induction and maintenance transition probabilities for all four therapies
(UST/IFX/ADA/CT) are sourced directly from Aliyev et al. 2019 Appendix S2, used at Aliyev's
native 2-week cycle (`R/00_derive_transition_probs.R`, `data/processed/DERIVATION_NOTES.md`).
Gate 2 engine in progress: `R/01_decision_tree.R` (induction split into the biologic/CT initial
occupancy vectors), `R/02_markov_engine.R` (cohort Markov core, with the M-S-to-CT switch and
2-year cap from analysis_plan.md §6.1/§6.4), and `R/03_cure_fraction_module.R` (the Treg
mixture-cure extension: week-56 landmark split, Sustained Deep Remission state, relapse hazard,
cap-aware relapse re-entry) are built, wired together, and tested end to end (induction ->
maintenance [-> cure branching for Treg], cohort-conserving over a full lifetime horizon).
Cost/utility attachment (`R/04_costs_utilities.R`): utility attachment and the arm-independent
per-cycle health-state monitoring cost landed 2026-08-04 (first pass). UST/IFX drug
acquisition + administration costs (induction and maintenance, both dose size and cycle-aligned
frequency, sourced from the STELARA/REMICADE Prescribing Information —
`data/raw/biologic_dosing_schedule.csv`) landed 2026-08-04 (second pass), deliberately excluding
Treg's non-cured track (efficacy-equivalent to UST, but never actually charged UST's drug cost —
see the module's header comment). ADA's dosing schedule was sourced and its drug-cost layer wired
in the same day (third/fourth pass, 2026-08-04 — HUMIRA Prescribing Information; priced via CMS
Medicaid's NADAC rather than ASP, since adalimumab is entirely self-administered SC with no
applicable Part-B ASP price for any of its doses). All three biologic comparators (UST/IFX/ADA)
now have both dosing and cost fully wired in. Treg's own one-time/two-dose acquisition cost is
still not wired in.

## Repository structure

- `data/raw/` — verbatim source extracts (Aliyev et al. 2019, ten Ham et al. 2020, CMS, HCUP,
  the manuscript's own tables). See `data/data_dictionary.md`.
- `data/processed/` — parameter values as the current `IBD_CEA_v6_PSA.xlsm` model actually uses
  them (a snapshot of that workbook, not independently re-derived).
- `docs/analysis_plan.md` — the health economic analysis plan (CHEERS 2022 item 4): objectives,
  model design, parameterisation, VOI plan, and the six decisions requiring sign-off.
- `docs/model_audit_v6.md` — audit findings (14 defects/undocumented choices) from the existing
  Excel workbooks.
- `docs/CHEERS_2022_checklist.md` — working CHEERS 2022 submission checklist.
- `docs/model_structure.md` — technical model spec, to be written alongside the `R/` engine.
- `R/` — analysis pipeline (decision tree, Markov engine, cure module, costs/utilities,
  deterministic results, PSA, EVPI/EVPPI, EJP), currently stubs only.
- `analysis/` — top-level entry points that run the `R/` pipeline.
- `tests/testthat/` — regression and invariant tests (transition rows sum to 1, cohort
  conservation, discounting checks).
- `output/` — generated figures and tables (not committed as final until a tagged release).

## Reproducing the analysis

Not yet runnable end-to-end. Once the `R/` pipeline is implemented:

```r
renv::restore()               # restore the locked package environment
source("analysis/run_full_analysis.R")   # regenerates every number, figure and table
```

`renv.lock` has not been generated yet (no R environment has been initialized in this
repository). Run `renv::init()` once the `R/` scripts have real package dependencies, and commit
the resulting lockfile — required for the journal's Data Availability Statement
(`docs/analysis_plan.md` §12.3).

## Licence

Code (`R/`, `analysis/`, `tests/`) is MIT-licensed. Data and documentation (`data/`, `docs/`,
this file) are CC-BY-4.0. See `LICENSE`.

## Citation

See `CITATION.cff`. Citation metadata is a placeholder pending the first Zenodo release tag.
