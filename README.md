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
Gate 2 engine in progress: `R/02_markov_engine.R` (cohort Markov core, with the M-S-to-CT switch
and 2-year cap from analysis_plan.md §6.1/§6.4) is built and tested. The induction decision tree
(`R/01`), mixture-cure module (`R/03`), and cost/utility attachment (`R/04`) are not yet written.

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
