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
now have both dosing and cost fully wired in. Treg's own one-time dose cost was sourced and wired
in the same day (fifth pass, 2026-08-04 — `treg_dose_cost()`): acquisition cost (ten Ham-derived,
$19,916.75) is fully traceable and wired in; infusion administration reuses the existing $57.90
rate (a shorter infusion doesn't change the billing code); cyclophosphamide preconditioning has a
sourced price but no sourced dose (defaults to $0, a known-incomplete placeholder); overnight
observation-stay cost **price is now fully sourced** — CMS's own CY2026 Addendum A gives C-APC
8011 (Comprehensive Observation Services) directly at **$2,672.15** (RW 29.2310, national
unadjusted), user-supplied 2026-08-05 and cross-validated to the cent against a $91.415
conversion factor independently back-calculated from the companion Addendum B file the same
day. Still wired in as an explicit opt-in `treg_dose_cost(observation_stay_cost_usd = ...)`
argument (default 0) rather than a base-case value — not because the price is uncertain, but
because whether a Treg infusion actually qualifies for this billing category (≥8 hours of
observation) is a separate, still-open clinical/modelling question. Single-dose base case only (analysis_plan.md §4.1) —
the 2-dose structural scenario needs per-cycle-trace logic not implemented here. A genuine,
previously-unflagged workbook inconsistency (A15, `docs/model_audit_v6.md`) was found and
documented, not silently resolved, while wiring the acquisition cost in.

`R/05_deterministic_results.R` (first pass, 2026-08-05): wires `R/00`–`R/04` together end to end
for all four arms (UST/IFX/ADA comparators, TREG intervention) at a 6.15-year horizon (Aliyev's
own 40×8-week span, converted to the native 2-week cycle) — a true lifetime horizon needs US
life-table mortality extrapolation, not sourced anywhere in this repo, and is a flagged gap rather
than silently approximated. Implements `run_base_case()` (UST/IFX/ADA plus TREG at π=0, the only
cure fraction with an actual number behind it per Decision 4 — labelled the Null floor scenario,
not a Treg efficacy claim) and, prioritised over the one-way tornado per analysis_plan.md §10.1's
own recommendation, the (π, price) **headroom frontier** (Aim 4): `headroom_pi_star()` solves via
root-finding for the minimum cure fraction at which Treg's NMB matches the best comparator's, and
`headroom_frontier()` sweeps it across a price grid. `analysis/run_base_case.R` is now a real, not
stub, entry point, writing `output/tables/base_case_results.csv`. 10 new tests (64 passing total).
Explicitly out of scope for this pass: the full one-way tornado (most Section 7.1 parameter
ranges aren't sourced yet), discount-rate/societal/refractory scenarios, and S1–S12 structural
scenarios — see the module's own header for the complete list. Also recorded there: Gate 2's
"reproduce Aliyev's published results" criterion (Aim 5) can't be fully closed from data already
in this repo — only Aliyev's transition-probability appendix was transcribed, not his own cost/
utility parameters or a results table to diff against, since this project deliberately re-sources
costs from current CMS pricing rather than his 2017 figures.

`R/06_psa.R` (first pass, 2026-08-05): 10,000-draw Monte Carlo PSA, sampling the two parameter
groups with an actual sourced distribution on file — the utility chain (Uniform bounds from
`data/processed/model_psa_parameter_distributions.csv`, sampled and chained multiplicatively per
analysis_plan.md §10.2's requirement) and Treg's acquisition price (Gamma via method-of-moments,
verified to reproduce analysis_plan.md item 12's own cited α≈15.4/scale≈1,297 figures exactly) —
plus the cure fraction π under Decision 4's own recorded Uniform(0,1) fallback. Everything else
(transition probabilities, relapse hazard h, comparator drug/monitoring costs) is held fixed at
its base-case value this pass; `sample_dirichlet_row()` is a general, tested utility ready for
transition-probability sampling once a sourced concentration parameter exists, but nothing calls
it yet — none of Aliyev's own published Table 3/4 rows come with a usable sample size in this
repo. Outputs: `psa_cost_effectiveness_plane()` and `psa_ceac()`. Deferred: the EJP posterior
(needs `R/08_ejp.R`, still a stub) and the CEAF (a presentation decision, not a computation this
module is blocked on). **New finding while building this: A16** (`docs/model_audit_v6.md`) — the
Remission utility `R/04`/`R/05` actually run on (0.9554, from `model_health_utilities.csv`)
doesn't match the 0.82 base value analysis_plan.md §7.1 item 19 and the PSA distribution file
itself cite for the same quantity. Doesn't block PSA sampling (a Uniform only needs correct
bounds, not a correct centre) but is unreconciled and flagged for co-author attention. Added
caching parameters (`induction_data`/`schedule`/`prices`, all optional) to
`R/05_deterministic_results.R`'s arm-runners to keep a 10,000-draw PSA's runtime reasonable (cut
wall-clock time by roughly a third by not re-reading the same three CSVs on every one of the
~40,000 arm-runner calls a full PSA makes); fully backward-compatible, all existing R/05 tests
still pass unchanged. 14 new tests (78 passing total).

`R/07_evpi_evppi.R` (first pass, 2026-08-05): EVPI and EVPPI for exactly the parameters
`R/06_psa.R` samples (π, Treg price, the utility chain) and their unions — EVPPI for a held-fixed
parameter is undefined by construction (A14), not a choice this module is deferring. Two distinct
EVPI framings, both requested by analysis_plan.md and kept explicitly separate: `evpi_surface()`
(§9.2) fixes Treg's price at each grid point (exploiting `treg_price_dependent_dose_cost()`'s
exact linearity in price — an algebraic cost swap, no re-simulation needed) and reports per-patient
EVPI over the remaining uncertainty, as a price × λ surface; `evppi_by_subset()` (§9.3) instead
treats price as one of the genuinely uncertain parameters (its actual sampled Gamma draws, no
override) and decomposes the REFERENCE-case total EVPI by subset (A/C/E and their unions — B/D/F/G
are undefined, same reason as above). `population_evpi()` implements the §9.4 formula but takes
`effective_population` as a required argument with no default — the fractions Decision 6 needs to
compute a real number aren't sourced yet, so none is invented. **Environment gap, not a design
choice:** analysis_plan.md §9.3 specifies the `voi` package (cross-checked against `BCEA::evppi()`,
SPDE-INLA for higher-dimensional subsets); neither `voi` nor `BCEA` can be installed in this
sandbox — their own dependencies (`earth`, `mvtnorm`) need `gfortran`, which isn't present, and
`INLA` isn't on CRAN at all. `evppi_gam()` implements the same underlying method (Strong, Oakley &
Brennan 2014 nonparametric regression) directly against `mgcv` instead, which ships with every
standard R install; `cross_check_voi()` is written and ready for an environment with a working
Fortran toolchain, but inert here. A real performance cliff was found and fixed while building
this: `mgcv::te()`'s default per-margin basis makes its total basis grow as roughly k^d — fine at
1-2 parameters, but a direct 4-parameter fit was still running after 5 minutes on 2,000 draws.
`reduce_for_gam()` PCA-reduces any subset over 2 raw parameters to 2 components first (per
analysis_plan.md §9.3's own suggestion), and `evppi_gam()` caps `te()`'s basis size as a backstop
for any direct caller that skips pre-reduction. **Known artifact of that shortcut, documented
where it's used:** a PCA-reduced superset can score BELOW an unreduced subset it contains (more
raw parameters diluted into the same 2 components can lose more signal than fewer parameters kept
at full resolution) — true EVPPI is monotone in the subset, this approximation isn't guaranteed to
be. 15 new tests (93 passing total).

`R/08_ejp.R` (first pass, 2026-08-05): deterministic and probabilistic economically justifiable
price (EJP, Aim 1), plus gross margin over COGS. `ejp_deterministic()` solves P\* in closed form
(analysis_plan.md §9.1: P\* = [C_comp + λ(Q_treg − Q_comp) − C_treg,non-drug] / D̄, D̄ = 1 under the
single-dose base case) rather than by root-finding — Treg's total cost is exactly linear in price,
so C_treg,non-drug is just `run_treg_arm_lifetime()` called at `price_usd = 0`. `ejp_frontier()`
sweeps this over a π grid; it and `R/05_deterministic_results.R`'s `headroom_frontier()` trace the
same (π, price) indifference curve from opposite directions, and round-trip through each other
exactly — the test suite's strongest check, and the one that caught a real bug during development
(an early version used `non_drug_cost` instead of `total_cost` for C_treg,non-drug, silently
dropping the $57.90 administration fee that `treg_price_dependent_dose_cost()` always adds
regardless of price — every P\* was off by exactly that amount until the round-trip test exposed
it). `ejp_probabilistic()` solves the same closed form per PSA draw as pure post-hoc algebra over
`R/06_psa.R`'s existing output — no new Markov simulation — reporting the median/95% CI of the
P\*_k density; `ejp_p50()` independently computes the distinct P_50 quantity analysis_plan.md §9.1
asks for (the price at which 50% of draws are still cost-effective, from the empirical survival
function, not from `median()`), since the two are only guaranteed to coincide under symmetry.
`gross_margin_over_cogs()` implements the §9.1 companion figure — (P\* − COGS)/P\*, COGS read from
`model_tenham_derived_treg_dose_cost.csv`'s own line — carrying forward the plan's explicit scope
limit that this is not a manufacturer-profitability claim. 10 new tests (103 passing total).

**Market-comparator context added 2026-08-05** (not a model parameter, not part of any Gate):
`data/raw/market_comparator_cell_therapy_prices.csv` records two real allogeneic-cell-therapy
prices — Ryoncil ($194,000/infusion, Mesoblast's FDA-approved MSC therapy) and
tabelecleucel/Ebvallo (ICER-recommended $143,900–$273,700/cycle) — for Discussion-section
context on how this study's EJP compares to what the market/HTA bodies have actually priced
allogeneic cell therapies at. `docs/analysis_plan.md` §9.1 now also specifies a **gross margin
over COGS** output (P\* vs. the ten Ham-derived $4,979.19 manufacturing cost) for when
`R/08_ejp.R` (Gate 8, still a stub) is built — explicitly scoped as a COGS-margin figure, not a
manufacturer-profitability claim, which would need private company financials this project
doesn't have.

**A16 correction + biosimilar comparator re-pricing (2026-08-05)**, done together in one branch
per `docs/decision_resolutions_2026-08-05.md` §11 ("blocks everything downstream; every number
moves" — the top two priorities in that memo, ahead of everything else in its Bucket 3). A16
(`docs/model_audit_v6.md`): the deterministic Remission utility (and the whole chain built on it)
was running on 0.9554396356, a stray live PSA draw captured in the workbook snapshot
(`model_health_utilities.csv`), not the 0.82 literature-cited base value analysis_plan.md §7.1
item 19 and the PSA distribution file both name as the source — four independent lines of evidence
in `docs/model_audit_v6.md`'s A16 entry. `R/04_costs_utilities.R`'s `load_health_state_utilities()`
now derives the deterministic vector from `model_psa_parameter_distributions.csv`'s sourced values
instead of reading the retired snapshot file, so the base case is by construction the PSA's own
central draw. `tests/testthat/test-parameter-provenance.R` (new) makes this class of defect (A15,
A16 — a `data/processed/` snapshot value used deterministically with no `data/raw/`-traceable
source) fail loudly rather than needing manual rediscovery. Separately, A15 (the ~0.9494×
dose-cost discrepancy) closed with the same root cause identified (a shared multiplicative shock
across an unrelated TREG/IFX price pair, not independently reproducible from either build-up) —
no code change needed, `psa_base` was already the figure in use. **Comparator pricing** (UST/IFX/
ADA, analysis_plan.md §4.2): re-extracted at a single stated date against real, directly-queried
CMS data (Part B ASP biosimilar Q-codes for UST/IFX; NADAC biosimilar NDCs for ADA) — biosimilar-
inclusive pricing is now the base case (UST induction $4.5375/mg, UST maintenance $60.85/mg
DERIVED PROXY — see below, IFX $2.6803/mg, ADA $14.2482/mg), originator pricing retained verbatim
as the S8 comparability-with-Aliyev scenario (`load_drug_prices(pricing_basis =
"originator_pre_biosimilar")`). UST maintenance's HCPCS code (J3357) no longer carries a Part B
ASP payment limit at all in the current fee schedule — independent, unplanned confirmation that
self-administered SC dosing isn't a Part B ASP product (the reasoning already applied to ADA); no
NADAC ustekinumab biosimilar NDC exists yet either, so that one figure is a documented proxy, not
a directly observed price (full chain in `data/raw/cms_asp_and_hcup_cost_sources.csv`'s own row).
**Every deterministic and probabilistic result produced before this fix is superseded, not merely
revised** — QALYs fall, ICERs and the required cure fraction π\* rise, and the EJP falls
substantially on both changes; re-run `analysis/run_base_case.R` and friends before citing any
number from before 2026-08-05. Full test suite: 324 assertions, 0 failures.

**Lifetime horizon implemented (2026-08-05)**, closing the gap `R/05_deterministic_results.R`'s
own module header had flagged ("no life-table data sourced anywhere in this repository").
`data/raw/nchs_us_life_tables_2021.csv` (NCHS "United States Life Tables, 2021", age- and
sex-specific `qx` by single year of age 0–100) plus new `R/utils/life_table.R`
(`load_life_table()`, `death_prob_schedule()`) and `R/utils/transition_matrix.R`'s new
`age_adjust_matrix()` REPLACE (not add to) Aliyev's own flat, non-age-varying embedded trial
mortality with the life-table figure for the cohort's current attained age — the "no CD excess
mortality" assumption analysis_plan.md §7.1 item 7 already specified, now actually wired in.
Male/female sub-cohorts (50/50) run separately against their own curve and sum exactly (linearity,
not an approximation) — `R/02_markov_engine.R`'s `run_maintenance_arm_with_mortality()` and
`R/03_cure_fraction_module.R`'s `run_treg_arm_with_mortality()`. **Real bug found and fixed along
the way:** the Sustained Deep Remission (cured) pool had no death exit at all before this —
harmless at a 6-year horizon, silently wrong (immortal cured patients) at a lifetime one; SDR now
has a competing-risks death/relapse split each cycle. `R/05_deterministic_results.R`'s
`HORIZON_CYCLES_LIFETIME` runs the cohort from baseline age 35 (Aliyev's own cohort mean) to the
life table's terminal age 100; `run_base_case()`, `headroom_pi_star()`, `headroom_frontier()` and
`R/08_ejp.R`'s EJP functions all gained an opt-in `baseline_age`/`life_table` pair (`NULL` default
reproduces the exact pre-existing 6.15-year/10-year behaviour, byte-identical, fully backward
compatible). This resolves the structural-infeasibility problem the 6.15-year horizon had — no
cure fraction at any plausible price made Treg cost-effective at that horizon; at the lifetime
horizon the required cure fraction π\* at Treg's sourced acquisition price ($19,916.75) is
feasible at all three WTP thresholds (see `output/tables/headroom_at_sourced_price_lifetime.csv`
and analysis_plan.md §4.3 for the exact figures). PSA/EVPI/EVPPI/probabilistic EJP still run at
the 6.15-year horizon this pass — extending `age_adjust_matrix()` to a ~40,000-call PSA is an
unbenchmarked performance question, deliberately deferred. `analysis/run_full_analysis.R`'s
"6/6: Lifetime horizon" section now runs this end to end. Full test suite: 385 assertions,
0 failures.

**h sweep + π prior-sensitivity on EVPPI implemented (2026-08-05)**, closing the resolutions
memo's own priority-order item 3 (`docs/treg-cd_decision_resolutions_2026-08-05.md` §3.2/§3.3),
immediately after the lifetime-horizon work above. Two independent additions:

- **h sweep.** The memo flagged the prior hardcoded `relapse_hazard_annual = 0` default as
  *anti*-conservative, not conservative — zero relapse is the single most favourable assumption
  available to Treg, not a cautious one. `R/03_cure_fraction_module.R`'s new
  `duration_to_hazard()`/`hazard_to_duration()` re-parameterise the relapse hazard by its implied
  median Sustained Deep Remission duration T (h = ln 2 / T) — "a median 10 years of drug-free
  remission" is a quantity a clinician or reviewer can reason about directly, a bare hazard rate
  is not. `R/05_deterministic_results.R`'s new `RELAPSE_DURATION_GRID_YEARS` (2, 5, 10, 20, Inf
  years; base case T = 10) and `headroom_frontier_by_duration()` sweep this jointly with price at
  the **lifetime** horizon (the memo: "h interacts with the horizon... must be run jointly — at
  6.15 years h barely matters; over a lifetime it dominates"), stacking one full
  `headroom_frontier()` per duration into a single (π, T, price) surface with an exact,
  per-row-simulated `qaly_gain` column (the "expected discounted QALY gain per treated patient"
  the memo recommends as the primary reporting axis, computed directly rather than only via the
  π·g(h) approximation below). No existing function's default behaviour changed — every prior call
  site with `relapse_hazard_annual` left implicit still gets 0, unchanged; T-based sweeping is
  additive.
- **π·g(h) factorisation, verified exact.** The memo asked for a numeric check that Treg's
  incremental QALY over its π=0 track is *approximately* π × g(h) (g(h) = the discounted QALY
  value of one cure), hedging it as "very good but not exact" since relapsed SDR patients re-enter
  the ordinary Markov trace. `R/05_deterministic_results.R`'s new `verify_pi_factorization()`
  checked this directly: **the factorisation is exact to floating-point precision (~1e-13 relative
  error, ordinary numerical noise) at every price and both horizons tried, not merely close** —
  because π enters `run_treg_arm()` exactly once, at the landmark split, and every operation
  downstream (transitions, relapse redistribution, the cap sweep, cost/utility attachment) is
  linear, with no π-dependent nonlinearity anywhere in the recursion. A stronger, more useful
  finding than the memo anticipated, now documented with a reproducible check rather than asserted
  from the structural argument alone.
- **π prior-sensitivity on EVPPI.** U(0,1) (Decision 4's own fallback prior for π) maximises prior
  variance on π, which the memo warns makes π come out as the dominant EVPPI parameter "more or
  less by construction" — an artefact a competent reviewer would flag unless checked against
  alternative priors. `R/07_evpi_evppi.R`'s new `prior_reweight()` computes importance weights
  turning the *same* 10,000-draw PSA already on file into a sample representative of Beta(1,3)
  (mass toward low cure fractions, consistent with the Ovasave/CATS1 experience) or Beta(2,2)
  (symmetric, informative) — no re-simulation, exactly the memo's own "this is cheap" framing.
  `evpi_from_nb()` and `evppi_gam()` both gained an optional `weights` argument (`NULL` default:
  byte-identical to their pre-existing unweighted behaviour); `evppi_by_subset()` threads a
  `weights` argument through consistently to both the total-EVPI denominator and every subset's
  EVPPI numerator; `evppi_prior_sensitivity()` runs all three priors (U(0,1) reference, Beta(1,3),
  Beta(2,2)) and reports a `rank_within_prior` column so the memo's stated finding — "the RANKING
  of EVPPI by parameter subset is the finding, the LEVEL is prior-dependent" — is directly readable
  off the returned table rather than requiring a second pass.

`analysis/run_full_analysis.R` gained two new sections (7/8, 8/8) wiring both pieces in: the
(π, T, price) headroom surface at the lifetime horizon plus the factorisation check
(`output/tables/headroom_frontier_by_duration_lifetime.csv`,
`output/tables/pi_factorization_check.csv`), and the prior-sensitivity EVPPI table
(`output/tables/evppi_prior_sensitivity.csv`) computed from the PSA draws already produced in
section 3/8 — no new PSA run. Full test suite: 461 assertions, 0 failures.

**Half-cycle correction implemented (2026-08-05)**, closing A12 (`docs/model_audit_v6.md`: "Cohort
state occupancy is applied at cycle boundaries with no correction") — CHEERS 2022 item 17, and
analysis_plan.md §4.1's own "Apply" recommendation. `R/04_costs_utilities.R`'s new
`half_cycle_weights()` applies the standard trapezoidal-rule correction to every continuously
state-occupancy-driven cost/QALY calculation (`trace_qalys()`, `trace_costs()`,
`attach_sdr_costs_utilities()`'s SDR accrual) — both endpoints of a trace weighted at 0.5, every
interior timepoint at full weight — but deliberately NOT to one-time or discrete-event costs
(induction/dose acquisition, a maintenance dose on its specific dosing cycle), which aren't the
continuously-accruing quantity the correction models. `half_cycle_correction = TRUE` is now the
default everywhere from `R/04` up through `R/05`'s `run_comparator_arm_lifetime()`/
`run_treg_arm_lifetime()` — every downstream caller (`R/06` PSA, `R/07` EVPI/EVPPI, `R/08` EJP)
inherits the correction automatically with no code changes of its own, since none of them override
the new default. `FALSE` reproduces the exact pre-2026-08-05 figures for comparability. **Every
deterministic and probabilistic result produced before this fix is superseded** (same class of
change as A16 and the biosimilar re-pricing) — the shift is small in relative terms but is real
and one-directional at the per-arm level: since occupancy and utility/cost values are always
non-negative, halving both endpoints' weight can only ever reduce a trace's summed QALYs/cost
relative to the uncorrected figure, never raise it (full reasoning: A12's own entry in
`docs/model_audit_v6.md`). `docs/analysis_plan.md` §4.1 updated. No new test files added -- existing ones extended (new
`half_cycle_weights()`/corrected-vs-raw checks) and two exact-value tests repointed to
`half_cycle_correction = FALSE` to keep testing the raw per-cycle formula they were written to
check. Full test suite: 474 assertions (up from 461), 0 failures.

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
