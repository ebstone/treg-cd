# Project Design and Analysis Plan

## Allogeneic T-Regulatory Cell Therapy for Moderate-to-Severe Crohn's Disease: Economically Justifiable Price and Value of Information

**Prepared for:** D. Jadambaa, E. Stone, B. Abraham (Johns Hopkins Bloomberg School of Public Health)
**Target journal:** *PharmacoEconomics* (Springer/Adis, journal 40273) — Original Research Article
**Document version:** 1.1 | 4 August 2026
**Status:** Gate 0 closed 4 August 2026 (see Section 15) — Decisions 1, 4 and 5 recorded as final; Decisions 2, 3 and 6 recorded provisionally pending full co-author sign-off. Gate 1 transition probabilities closed 4 August 2026 (Aliyev et al. 2019 Appendix S2 sourced directly — see §14); repository foundations otherwise in progress.

---

## 0. Executive summary

### 0.1 What this redesign changes

| # | Change | Why |
|---|---|---|
| 1 | **Reframe the paper** from "here is the Treg ICER" to "here is what Treg must deliver, and at what price, to be worth it — and what we should measure next." Primary outputs become EJP, EVPI and EVPPI; ICER/NMB become supporting outputs. | An ICER computed from invented efficacy parameters is not a finding. A price threshold and a research-prioritisation result are defensible outputs for a technology with no clinical data. |
| 2 | **Replace the "20% better than ustekinumab" adjustment with an explicit mixture-cure structure**: a cure fraction π entering a Sustained Deep Remission state with its own relapse hazard, and non-cured patients modelled as efficacy-equivalent to ustekinumab. | The current parameterisation is structurally "a slightly better biologic," which does not represent the mechanism the paper argues for and is not how potentially curative cell therapies are modelled in HTA. |
| 3 | **Add a required-effectiveness (headroom) analysis as a co-primary result**: the minimum cure fraction, and the (π, price) frontier, at which Treg is cost-effective. | This converts the largest unknown from a fabricated input into the reported output. It is the single most important defensive move against the "you made the numbers up" critique. |
| 4 | **Extend the time horizon to lifetime**, retaining 6.15 years as a comparability scenario. | A 6-year horizon structurally cannot capture the value of a cure. Reviewers will raise this immediately. |
| 5 | **Rebuild in R with a full PSA** over transition probabilities, cure parameters, costs and utilities (the current PSA varies only 7 parameters), and deposit code and data in a public, Zenodo-archived repository. | EVPI/EVPPI are undefined over parameters that are not sampled. Also satisfies the journal's mandatory Data Availability Statement. |

### 0.2 What must be fixed regardless of the redesign

An audit of `IBD_CEA_v6_PSA.xlsm` and `IBD_CEA_v6_Univariate_Sensitivity_Analysis.xlsm` identified **fourteen defects or undocumented choices**, several of them material and at least three of which bias results in favour of Treg. Full detail with cell references is in **Appendix A**. The most consequential:

- **Costs are discounted against the wrong quantity.** In all three Markov sheets the cost discount factor is `=AQn*(1/(1+discount)^Nn)`, where column N is *cumulative conventional-therapy deaths*, not elapsed years (column A). QALYs are discounted correctly against column A. Costs are therefore materially under-discounted, which inflates the present value of the comparators' recurring drug costs relative to Treg's front-loaded costs.
- **The base-case Treg Surgery-row transition probabilities appear column-shifted** relative to ustekinumab and to the workbook's own Alpha/Beta scenario sheets.
- **The Treg arm is charged health-state costs in 2017 dollars** ($217 / $91 / $10) while ustekinumab and infliximab are charged 2025 dollars ($282.86 / $118.62 / $13.04).
- **The PSA varies only 5 utility parameters and 2 drug prices.** No transition probability and no efficacy parameter is sampled.
- **The Economically Justifiable Price worksheet is disconnected from the model**, using hardcoded QALYs (26.48 vs 25.73) that appear in no other sheet and returning $185,625 rather than the $95,382 reported in the manuscript.

None of these invalidate the project. All of them are reasons to rebuild the engine in code rather than patch the workbook.

### 0.3 What needs your sign-off

Six decisions (Section 15). The two that change the headline result most are **(1)** whether the 2-year biologic-maintenance cap is reinstated and **(4)** what anchors the base-case cure fraction. Nothing downstream should be coded until both are settled.

---

## 1. Title and elevator pitch

### 1.1 Proposed titles

**Recommended:** *Early Economic Evaluation of a Hypothetical Allogeneic T-Regulatory Cell Therapy for Moderate-to-Severe Crohn's Disease: Economically Justifiable Price and Value of Information*

**Alternative (higher-impact framing):** *What Would a Curative Cell Therapy for Crohn's Disease Have to Deliver? A Headroom, Pricing and Value-of-Information Analysis*

The first is safer for *PharmacoEconomics*' house style and signals the methods clearly. The second reads better and is more honest about what the paper actually does; it risks looking speculative to a conservative handling editor.

### 1.2 Elevator pitch

Allogeneic, "off-the-shelf" regulatory T-cell therapies are entering first-in-human testing in Crohn's disease, and they are being developed on a premise no existing therapy makes: that a single infusion could durably reset gut immune tolerance rather than suppress inflammation indefinitely. No efficacy data exist. Rather than assert an efficacy assumption and report an incremental cost-effectiveness ratio built on it, this study inverts the question. Using a hybrid decision-tree/Markov model extended with a mixture-cure structure — anchored on the ustekinumab and infliximab evidence base assembled by Aliyev et al. (2019) and the ten Ham et al. (2020) cell-therapy manufacturing costing framework — we ask three things a developer, payer or funder can act on today. What is the maximum price at which such a therapy would be cost-effective at conventional US willingness-to-pay thresholds? What fraction of patients would it have to durably cure, at a given price, to clear that bar? And where is the value of further evidence greatest — in pinning down the durability of remission, or in driving down manufacturing cost? By reporting the economically justifiable price alongside expected value of perfect and partial perfect information, the paper turns an unavoidable absence of data into an explicit research-prioritisation result, at the exact moment the field's first efficacy readouts are expected.

---

## 2. Structured abstract shell

Target 380–430 words. The journal specifies 150–250 words but explicitly permits up to 450 where needed for full reporting-guideline compliance, which CHEERS 2022 will require here. Results are bracketed placeholders pending re-analysis.

> **Background and Objective.** Allogeneic regulatory T-cell (Treg) therapies are in early clinical development for Crohn's disease on the premise of durable, potentially curative immune tolerance, but no efficacy data exist. We aimed to establish (i) the economically justifiable price (EJP) of such a therapy, (ii) the durable-remission ("cure") fraction it would need to achieve to be cost-effective at a given price, and (iii) where additional evidence would be most valuable, via expected value of perfect and partial perfect information (EVPI, EVPPI).
>
> **Methods.** We developed a hybrid decision-tree and cohort Markov model (8-week cycles; [lifetime] horizon; 3% annual discounting) comparing a hypothetical single-infusion allogeneic Treg therapy with ustekinumab, infliximab [and adalimumab] in [biologic-naïve] adults with moderate-to-severe Crohn's disease, from the US healthcare-sector perspective. Health states were Remission, Mild, Moderate-Severe, Moderate-Severe Responder, Surgery and Death, defined by Crohn's Disease Activity Index score. Biologic induction and maintenance transition probabilities were derived from Aliyev et al. (2019). The Treg arm was extended with a mixture-cure structure: a fraction of week-52 responders entered a Sustained Deep Remission state with a low annual relapse hazard; the remainder were assumed efficacy-equivalent to ustekinumab. Acquisition cost was derived from the ten Ham et al. (2020) manufacturing costing framework. Because no Treg efficacy data exist, the cure fraction and relapse hazard were treated as unknowns and varied across their full plausible range, anchored on three published analogs. Probabilistic analysis used [10,000] Monte Carlo draws; EVPPI was estimated by nonparametric regression.
>
> **Results.** [ICER, NMB and probability of cost-effectiveness by arm.] The EJP for Treg therapy was [$X] (95% credible interval [$X–$X]) at a $150,000/QALY threshold. At an acquisition price of [$X], Treg required a durable cure fraction of at least [X%] to be cost-effective. Population EVPI was [$X] per year, of which [X%] was attributable to the cure fraction and post-cure relapse hazard jointly and [X%] to acquisition cost.
>
> **Conclusions.** Cost-effectiveness of allogeneic Treg therapy in Crohn's disease is governed jointly by durability of remission and manufacturing cost, and the results define explicit, testable targets for both. [Directional conclusion.] The value of information attributable to durability greatly exceeds that attributable to price uncertainty, indicating that long-term follow-up of early-phase cohorts, rather than further cost engineering, is the higher-value next investment.

---

## 3. Research objectives

Stated as falsifiable aims. Each primary aim maps to one primary endpoint.

**Aim 1 (EJP).** Estimate the maximum per-dose acquisition price P\* at which a hypothetical allogeneic Treg therapy achieves an incremental net monetary benefit of zero versus the next-best non-dominated comparator, at willingness-to-pay thresholds of $50,000, $100,000 and $150,000/QALY. Report P\* both deterministically and as a full posterior distribution from the probabilistic analysis.
*Falsifiable as:* P\* is or is not greater than the $19,917/dose estimate derived from the ten Ham costing framework.

**Aim 2 (EVPI).** Estimate per-patient and population-level expected value of perfect information for the choice among comparators, as a function of Treg price and WTP threshold.
*Falsifiable as:* population EVPI at $150,000/QALY does or does not exceed the plausible cost of a confirmatory trial programme.

**Aim 3 (EVPPI).** Estimate EVPPI for pre-specified parameter subsets and rank them.
*Falsifiable as:* the durability block (cure fraction + relapse hazard) does or does not carry higher EVPPI than the cost block.

**Aim 4 (headroom — recommended addition as co-primary).** Identify the minimum durable cure fraction π\* required for cost-effectiveness across the plausible price range, and present the resulting (π, price) cost-effectiveness frontier.
*Rationale:* this is the one primary result that requires no invented efficacy input, and it directly answers the developer's question.

**Aim 5 (secondary — comparability).** Report ICERs, NMB and cost-effectiveness acceptability curves for all comparators under base-case assumptions, and demonstrate that the reconstructed biologic arms reproduce Aliyev et al. (2019) and are consistent with NICE TA456, as external validation.

**Aim 6 (optional extension).** Estimate expected value of sample information (EVSI) for a Treg trial of the size currently being run in the field (n ≈ 39), to quantify what a Phase 1/2a readout can and cannot resolve. Recommended only if Aims 1–5 land comfortably inside the word limit; EVSI is a strong differentiator for *PharmacoEconomics* but is a substantial additional coding effort.

---

## 4. Decision problem and perspective

### 4.1 Design parameters: current draft vs. recommendation

| Element | Current draft | Recommendation | Rationale |
|---|---|---|---|
| Population | Biologic-naïve, mean age 35, 71 kg, 50% male, moderate-severe CD (from Aliyev) | Keep as base case; add refractory scenario | See Section 5 and Decision 3 |
| Comparators | UST, IFX, (CT as post-failure pathway) | **Add ADA**; retain CT as the post-failure pathway, not a comparator arm | See Decision 2 |
| Intervention | Hypothetical allogeneic Treg, 2 doses | **Single infusion base case**; second dose as scenario | Aligns with the single-infusion design of the allogeneic Tr1 programme now in Phase 1/2a; also the more favourable and more likely commercial configuration |
| Perspective | US formal healthcare sector | Keep; **add societal scenario** | US Second Panel recommends reporting both. The Manceur et al. (2020) absenteeism/disability data are already extracted in the workbook's *Cell therapy costs* sheet, so the marginal effort is low, and productivity effects materially favour a durable therapy |
| Time horizon | 40 × 8-week cycles = 6.15 years | **Lifetime**, with 6.15-year and 10-year scenarios | A cure model with a 6-year horizon truncates precisely the benefit the paper is about. Requires US life tables and an explicit extrapolation assumption beyond trial-supported data — **implemented 2026-08-05, see below** |
| Cycle length | 8 weeks | Keep | Matches q8w maintenance dosing for both UST and IFX; changing it would require re-deriving every transition probability |
| Half-cycle correction | None | **Apply** (life-table or Simpson's 1/3) | CHEERS 2022 item 17; trivial to implement in code — **implemented 2026-08-05**, `R/04_costs_utilities.R`'s `half_cycle_weights()` (trapezoidal method); see A12, `docs/model_audit_v6.md`, for detail |
| Discount rate | 3% costs and QALYs | Keep 3%; scenarios at 0%, 1.5% and 5% | 1.5% is worth including explicitly given the live methodological debate about discounting one-time potentially curative therapies |
| WTP threshold | $150,000/QALY | Report $50k, $100k, $150k | $150,000 alone reads as the most permissive choice available |
| Currency / price year | 2025 USD (partially applied) | **2025 USD throughout**, index documented | See Section 8 |
| Cohort size | 1,000 | Report per-patient results | Aggregate-per-1,000 reporting in the current draft (e.g. NMB of $529,388) is a per-1,000 figure presented without units and will confuse reviewers |

### 4.2 A note on the comparator price base

Both comparators are now in a biosimilar market. Infliximab biosimilars have been available for several years; ustekinumab biosimilars entered the US market during 2025. The current model charges roughly $14,029 per 8-week ustekinumab maintenance cycle. **Action item:** re-extract ASP for all comparators at a single, stated pricing date and run a biosimilar-pricing scenario. If comparator prices have fallen materially, the EJP for Treg falls with them, and a reviewer will notice if this is not addressed. This is not currently reflected anywhere in the model or manuscript.

**Resolved 2026-08-05 (`docs/decision_resolutions_2026-08-05.md` §8; elevated to co-priority alongside A16, §0/§11 of that memo — "the EJP essentially *is* the net present value of the comparator drug spend a one-time therapy displaces," so comparator price is the single largest lever on every headline number).** This was under-scoped as a routine sourcing task; it is a base-case error. Re-extracted UST, IFX and ADA acquisition prices at a single stated date (2026-08-05), applying the same benefit-routing logic already used for ADA (Part B ASP where an ASP exists; NADAC where it doesn't) product-by-product rather than uniformly:

- **UST induction**: median CMS Part B ASP across the 6 ustekinumab biosimilar Q-codes now carrying one — **$4.5375/mg**, down from $12.808/mg originator.
- **UST maintenance**: J3357 (the code this project's pricing used) **no longer carries a Part B ASP payment limit at all** in the current fee schedule — direct, unplanned confirmation that self-administered SC maintenance dosing is not a Part B ASP product, resolving the nuance flagged when this item was scoped (the $4,695 Stelara IRA Maximum Fair Price is a Medicare Part D, 30-day-supply figure, not a $/mg number on the same benefit pathway this project prices UST from; it was not used for exactly that reason). NADAC does not yet list a ustekinumab biosimilar NDC either, so this figure is a **derived proxy** (the UST induction biosimilar:originator discount ratio applied to the prior SC originator ASP figure) rather than a directly observed price — **$60.85/mg**, down from $155.883/mg originator; flagged for revisiting once NADAC or a future ASP file actually lists an ustekinumab biosimilar SC NDC.
- **IFX**: median across 3 infliximab biosimilar Q-codes — **$2.6803/mg**, a modest ~14% cut from originator ($3.1041/mg the same quarter), consistent with a biosimilar market that converged years ago rather than one still in an early price war.
- **ADA**: median CMS NADAC across 4 adalimumab biosimilar 40mg-pen NDCs — **$14.2482/mg**, an ~83% cut from the brand-only NADAC figure ($84.1678/mg) this project sourced 2026-08-04.

**Biosimilar-inclusive pricing is now the base case; S8 (§10.3) inverts to make originator (pre-biosimilar) pricing the comparability-with-Aliyev scenario instead** — `R/04_costs_utilities.R`'s `load_drug_prices()` takes a `pricing_basis` argument (`"biosimilar_2026"` default, `"originator_pre_biosimilar"` for the scenario); full derivations are in `data/raw/cms_asp_and_hcup_cost_sources.csv`'s own rows. As anticipated, this already reshuffles which comparator is "next-best non-dominated" (§9.1) — IFX, not UST or ADA, is now typically the cheapest comparator at the base-case weight, changing which arm the EJP is solved against. All deterministic and probabilistic results generated before this fix are superseded, not merely revised, exactly as for A16 — the two were deliberately corrected together in one branch for that reason (`docs/decision_resolutions_2026-08-05.md` §11).

### 4.3 Lifetime horizon: implementation note (2026-08-05)

The lifetime horizon this section recommends is now implemented, closing the gap `R/05_deterministic_results.R`'s own module header had flagged ("no life-table data sourced anywhere in this repository"). Background mortality is US life-table, age- and sex-specific (item 7, §7.1), sourced from NCHS "United States Life Tables, 2021" (`data/raw/nchs_us_life_tables_2021.csv`) and applied via `R/utils/life_table.R`'s `death_prob_schedule()` and `R/utils/transition_matrix.R`'s `age_adjust_matrix()`, which **replaces** (not adds to) each maintenance matrix's Death-column entry with the life-table figure for the cohort's current attained age — Aliyev's own embedded trial mortality (~0.00006/2wk, flat and non-age-varying) is realistic for his short young-cohort follow-up but structurally wrong once the horizon runs into old age, which is the entire point of extending it (per this section's own "no CD excess mortality" framing, i.e. total mortality = general-population background only). Male and female sub-cohorts (50/50, per §5) are run separately against their own life-table curve and summed — exact by linearity, not an approximation. `HORIZON_CYCLES_LIFETIME` (`R/05_deterministic_results.R`) runs the cohort from baseline age 35 (§5) to the life table's terminal age 100. `run_base_case()`, `headroom_pi_star()`, `headroom_frontier()` and `R/08_ejp.R`'s EJP functions all take an opt-in `baseline_age`/`life_table` pair; `NULL` (the default) reproduces the exact pre-existing 6.15-year/10-year behaviour unchanged, so this is additive, not a breaking change to the earlier scenarios this section names.

This resolves the structural-infeasibility problem the 6.15-year horizon had: at that horizon, no cure fraction at any plausible Treg acquisition price could make Treg cost-effective against comparators, because six years is too short for even a full durable cure to recoup a five-figure one-time price against comparators' recurring drug spend. At the lifetime horizon, the required cure fraction π\* at Treg's sourced acquisition price ($19,916.75, `data/processed/model_tenham_derived_treg_dose_cost.csv`) is feasible at all three WTP thresholds — see `output/tables/headroom_at_sourced_price_lifetime.csv` for the exact figures; §10.1 covers how the (π, price) frontier is reported generally, §10.3's S5 is the lifetime/10-year/6.15-year structural scenario this closes. PSA, EVPI/EVPPI and probabilistic EJP still run at the 6.15-year horizon this pass — extending `age_adjust_matrix()` to a ~40,000-call PSA is an unbenchmarked performance question, deliberately deferred (`R/05_deterministic_results.R`'s own module header).

---

## 5. Population

### 5.1 The tension

The comparator evidence base is biologic-naïve. Aliyev et al. derived transition probabilities from UNITI/IM-UNITI, ACCENT I, CLASSIC-I and CHARM in biologic-naïve or predominantly biologic-naïve populations, and the current draft inherits a cohort of mean age 35, 71 kg, 50% male.

The Treg evidence base — such as it is — is not. The allogeneic Tr1 programme now in Phase 1/2a in Crohn's disease is enrolling adults aged 18–65 with **moderate-to-severe treatment-refractory** disease, having failed multiple prior advanced therapies. That is a population with lower spontaneous remission rates, higher surgery risk, higher baseline costs and a different comparator set (the relevant comparator is not "first biologic" but "next advanced therapy or surgery").

Borrowing structural plausibility from a refractory programme while modelling a biologic-naïve cohort is defensible only if stated plainly, and it will be a reviewer target if it is not.

### 5.2 Recommendation — Decision 3

**Run both (option c), with biologic-naïve as the base case.**

- *Base case: biologic-naïve.* Preserves direct comparability with Aliyev et al. and NICE TA456, preserves the validity of every comparator transition probability, and supports Aim 5.
- *Co-primary scenario: refractory.* Reported in full in the main text, not relegated to supplementary material, because it is the population in which the technology would actually first be used and priced.

For the refractory scenario, do **not** attempt to re-derive refractory transition probabilities from first principles. Instead apply an explicit, sourced set of multipliers to the biologic-naïve probabilities (reduced response and remission probabilities, elevated surgery hazard) drawn from published second-line and post-anti-TNF-failure CD cohorts, and vary those multipliers widely. State the approach as a scenario built on assumption, not as a separately estimated model.

**What would change this recommendation:** if a co-author judges that the paper's contribution is specifically to the developer's near-term decision rather than to the HTA literature, the refractory population becomes the better base case and biologic-naïve becomes the comparability scenario. That is a defensible inversion; it should be a deliberate choice, not a default.

### 5.3 Generalisability language for the manuscript

The manuscript must state explicitly that no Treg parameter is derived from the ongoing Crohn's trial, that the trial population differs from the modelled base-case population, and that the trial is cited only to establish that the modelled technology class is in clinical development. This belongs in Methods, not only in Limitations.

### 5.4 Refractory multipliers: implementation note (2026-08-05)

§5.2's recommended method — explicit, sourced multipliers rather than a re-derived refractory
model — is now implemented for response and remission (`R/utils/refractory_multipliers.R`,
`R/05_deterministic_results.R`'s `run_refractory_scenario()`). Source: UNITI-1 vs. UNITI-2
(Feagan et al. 2016, NEJM) — the same trial programme §5.1 already names as (part of) Aliyev's
own naive-population evidence base, which is exactly why this pair was chosen over an unrelated
second-line cohort: same drug, same doses, same endpoints, same investigators, differing only in
which patients were eligible (UNITI-1: primary/secondary anti-TNF nonresponse; UNITI-2:
conventional-therapy failure, 68.6% anti-TNF-naive). Two multipliers, both ustekinumab-specific
but applied to IFX/ADA's matrices too as a stated cross-drug proxy (no equivalently clean paired
refractory/naive trial design exists for either): induction response/remission (week 8, 6 mg/kg
IV) and maintenance remission (IM-UNITI week 44, 90 mg SC q8w, converted from its 36-week
cumulative figure to a per-2-week-cycle multiplier before use — see that file's own header for
why). Applied to UST/IFX/ADA's own induction and maintenance matrices only; CT's matrix and
Treg's own parameterisation are unaffected (Treg's reference clinical programme is already a
refractory population, per §5.1, so nothing about it needs adjusting for this scenario).

**Explicitly not included in this pass, a real flagged gap, not an oversight**: elevated surgery
hazard (the other half of §5.2's own recommendation). No equivalently clean paired
refractory-vs-naive surgery-rate data was found in the same trial programme in the time available;
a future pass could source one separately (e.g. from a dedicated post-anti-TNF-failure surgical
cohort) rather than reusing UNITI's own surgery data, which the induction trials weren't powered
or designed to report by prior-treatment-failure subgroup.

**Result, 6.15-year horizon** (`output/tables/refractory_scenario_results.csv`,
`analysis/run_full_analysis.R` step 9): QALYs fall for all three comparators under the refractory
adjustment, as expected (worse response/remission). Total cost also falls, for all three — a real,
non-obvious mechanical consequence, not a bug: patients who fail to achieve or hold remission
switch to the CT track more often, and CT is markedly cheaper than continuing branded-biologic
therapy, so the cost saved by more frequent CT-switching outweighs the extra cost of more relapse/
non-response cycles at this project's current CT/biologic price gap. Net effect on NMB is small
and mixed by comparator/WTP in this pass's numbers — worth a sentence in Results, not a headline
finding on its own.

---

## 6. Model structure

### 6.1 Overview

Two linked components, retained from the current design:

1. **Decision tree (induction, 8 weeks).** Each arm's cohort is partitioned at the end of induction into Remission, Mild, Moderate-Severe Responder (M-SR), or non-response. Non-responders enter the conventional-therapy (CT) Markov and remain on CT-derived transition probabilities.
2. **Cohort Markov (maintenance).** 8-week cycles. Living states: Remission, Mild, M-SR, Moderate-Severe (M-S), Surgery; absorbing state: Death. Patients on a biologic who deteriorate into M-S are switched to the CT track at the end of that cycle (this rule *is* correctly implemented in the current workbook).

The redesign adds one state and one transition rule, applicable only to the Treg arm.

### 6.2 The mixture-cure extension

**New state: Sustained Deep Remission (SDR).** Distinct from ordinary Remission. Characterised by:

- No ongoing drug acquisition cost.
- Reduced monitoring cost (recommend: remission-state monitoring cost, halved after year 2 — flag as assumption).
- Utility equal to the Remission utility in the base case (scenario: SDR utility set to general-population utility for age, on the argument that a durable cure removes the residual disutility of managed disease — this is a meaningful scenario because it is where much of the cure's QALY value lives).
- A low **annual relapse hazard** h, converted to an 8-week probability. Relapsed patients re-enter the ordinary Markov (base case: into Mild; scenario: into M-SR).
- **Cured patients never transition to the CT track.** Mechanistically, "flowing into conventional therapy" represents exhaustion of a maintenance drug; it has no meaning for a patient who is not on one. This resolves the second half of Decision 1.

**Where the cure split happens — recommended: a 52-week landmark.** Rather than splitting responders into cured/not-cured at the end of induction, apply the split at cycle 7 (week 56, the nearest cycle boundary to 52 weeks), conditional on the patient being in Remission at that point:

π = P(enter SDR | in Remission at week 56)

*Implementation note, 2026-08-04:* "cycle 7" above is the original 8-week-cycle design's figure; at this project's native 2-week cycle (§4.1's 8-week design was dropped — see `R/00_derive_transition_probs.R`'s header comment) week 56 is cycle 28. Also confirmed with E. Stone, not otherwise specified above: a Sustained Deep Remission relapse that occurs *after* the 2-year cap (§6.4) has already fired routes directly to CT rather than back to the (no-longer-available) biologic track — consistent with the cap's own "exhaustion of a maintenance drug" rationale, and with why the cap applies to non-cured Treg responders at all (Decision 1: avoiding an unearned advantage over UST/IFX/ADA, which would otherwise recur if only relapsers kept indefinite post-cap biologic-track access). Implemented in `R/03_cure_fraction_module.R`.

**Rationale.** (i) It matches how mixture-cure fractions are identified empirically — from a landmark of sustained response, not from initial response. (ii) It matches the pattern in the analogs, where early response and durable remission clearly dissociate. (iii) It keeps the cure fraction interpretable as a quantity a real trial could estimate, which is precisely what makes the EVPPI on it meaningful. (iv) It avoids the biologically implausible implication that curative status is determined at week 8.

**Non-cured Treg patients: recommend efficacy-equivalent to ustekinumab.** Drop the 20% adverse-transition reduction entirely. All incremental benefit in the Treg arm then flows through π and h, which makes the model interpretable, makes the headroom analysis clean, and prevents double-counting of benefit. The current approach — a cure narrative in the text and a 20% efficacy shave in the model — is the mismatch this redesign exists to fix. If co-authors wish to retain some non-cure advantage, model it explicitly as a hazard ratio on the M-S and Surgery transitions, sampled in the PSA, rather than as a deterministic 20% redistribution.

### 6.3 State-transition diagram

```
                       INDUCTION (decision tree, 8 weeks)
                       ─────────────────────────────────
                                  Cohort
                                    │
              ┌─────────────────────┼─────────────────────┐
              │                     │                     │
           UST arm               IFX arm            TREG arm   [+ ADA arm]
              │                     │                     │
        ┌─────┴─────┐         ┌─────┴─────┐         ┌─────┴─────┐
     Response   No response  Response  No resp.  Response   No response
        │            │           │         │         │            │
   (Rem/Mild/     CT track  (Rem/Mild/  CT track (Rem/Mild/    CT track
     M-SR)                     M-SR)               M-SR)


                    MAINTENANCE (Markov, 8-week cycles)
                    ───────────────────────────────────

   BIOLOGIC TRACK (UST / IFX / ADA)          CT TRACK
   ┌───────────────────────────┐             ┌───────────────────────────┐
   │  Remission ⇄ Mild ⇄ M-SR  │             │  Remission ⇄ Mild ⇄ M-SR  │
   │      │        │      │    │             │      │        │      │    │
   │      └────► Surgery ◄┘    │             │      └───► Surgery ◄┘     │
   │               │           │             │             │             │
   │            [M-S] ─────────┼────────────►│           [M-S]           │
   └───────────────────────────┘  switch to  └───────────────────────────┘
        │                         CT on loss        │
        │                         of response       │
        └──────────────► Death ◄────────────────────┘
                     (all states, background
                      mortality; no CD excess)

   TREG TRACK  (adds the cure structure)
   ┌──────────────────────────────────────────────────────────────┐
   │                                                              │
   │   Remission ⇄ Mild ⇄ M-SR ──► [M-S] ──► CT track             │
   │       │                                                      │
   │       │ at week-56 landmark, split with probability π        │
   │       ▼                                                      │
   │  ┌─────────────────────────────┐                             │
   │  │  SUSTAINED DEEP REMISSION   │  no drug cost               │
   │  │  (SDR)                      │  reduced monitoring         │
   │  └──────────┬──────────────────┘  never enters CT            │
   │             │ annual relapse hazard h                        │
   │             ▼                                                │
   │          Mild  (re-enters ordinary Markov)                   │
   └──────────────────────────────────────────────────────────────┘
```

For the manuscript figure, render this in a vector tool at 174 mm width per journal artwork specifications; the ASCII form above is the specification, not the deliverable.

### 6.4 The 2-year maintenance cap, by arm — Decision 1

Aliyev's source model caps biologic maintenance at 2 years, after which all remaining patients move to CT-derived transition probabilities. **This cap is absent from the current workbook.** Because the model runs 40 cycles rather than 13, patients currently remain on biologic-specific probabilities — and biologic-specific costs, at roughly $14,029 per cycle for ustekinumab — for the full horizon.

| Arm | Recommendation | Rationale |
|---|---|---|
| UST / IFX / ADA | **Reinstate the cap in the base case** (switch to CT-derived probabilities and CT costs after cycle 13); no-cap as a structural scenario | Restores fidelity to the source model whose parameters we are borrowing; without it we are using Aliyev's numbers in a structure Aliyev did not validate. Expect a material reduction in comparator costs and therefore in Treg's EJP |
| TREG — cured (SDR) | **No cap. Ever.** | A cured patient is not on maintenance therapy; there is nothing to discontinue |
| TREG — non-cured responders | Apply the same cap as the biologics in the base case | If non-cured patients are UST-equivalent, they should be treated identically; asymmetry here would be an unearned advantage |

**Be aware this decision moves the headline result substantially and in the direction unfavourable to Treg.** With the cap, comparator arms shed years of drug cost; Treg's cost advantage shrinks and its EJP falls. That is the correct answer if we are borrowing Aliyev's parameters, but co-authors should go into it with eyes open, and the no-cap scenario must be reported alongside.

**What would change this recommendation:** evidence that real-world persistence on ustekinumab or infliximab in responders substantially exceeds 2 years — which it plausibly does. A defensible alternative is to replace the arbitrary 2-year cap with a **sourced time-to-discontinuation curve** applied to all maintenance arms. That is better science and more work; if a co-author has access to a persistence estimate from claims data, it is worth the effort and would pre-empt the most likely reviewer objection to either choice.

---

## 7. Data sources and parameterisation plan

### 7.1 Parameter inventory

Distributions follow standard practice: Dirichlet for multinomial transitions out of a state (replacing Aliyev's independent Beta/Uniform treatment, which does not preserve row sums), Gamma for costs, Beta for probabilities and utilities on (0,1), log-normal for hazards.

| # | Parameter group | Base-case source | Distribution | Status |
|---|---|---|---|---|
| **Comparator efficacy** |
| 1 | UST induction transitions (M-S → Rem/Mild/M-SR/M-S) | Aliyev 2019 Suppl. Table 3, DEALE-converted to 8 weeks | Dirichlet | **Sourced** — re-derivation required (Appendix A, finding 5) |
| 2 | IFX induction transitions | Aliyev 2019 Suppl. Table 3 | Dirichlet | Sourced |
| 3 | ADA induction transitions | Aliyev 2019 Suppl. Table 3 | Dirichlet | Sourced — only if Decision 2 is yes |
| 4 | UST / IFX / ADA maintenance transitions | Aliyev 2019 Suppl. Table 4 | Dirichlet | Sourced |
| 5 | CT maintenance transitions | Aliyev 2019, IM-UNITI placebo arm | Dirichlet | Sourced |
| 6 | Surgery transitions (out of Surgery state) | Aliyev 2019 Suppl. Table 2 (8-week probabilities 0.060 / 0.075 / 0.34) | Dirichlet | Sourced |
| 7 | Background mortality | US life tables, age- and sex-specific; no CD excess mortality (per Aliyev / TA456) | Fixed | **Sourced and implemented (2026-08-05)** — NCHS 2021 US life tables (`data/raw/nchs_us_life_tables_2021.csv`), wired through the engine via `R/utils/life_table.R`; see §4.3 |
| **Treg efficacy — the unknowns** |
| 8 | Treg induction response and health-state split | Base case: **equal to UST** | Dirichlet | **Assumption** — recommended change from current 20%-shave |
| 9 | Treg cure fraction π (P(SDR \| Remission at week 56)) | See Section 7.2 | Beta, wide | **Assumption. Also a reported output (Aim 4)** |
| 10 | Post-cure annual relapse hazard h | See Section 7.2 | Log-normal, wide | **Assumption** |
| 11 | Non-cured Treg efficacy | Base case: **equal to UST**; scenario: HR on M-S and Surgery transitions | Log-normal on HR | **Assumption** |
| **Costs (2025 USD)** |
| 12 | Treg acquisition cost per dose | ten Ham et al. 2020 costing framework as operationalised in the workbook's *Treg Cost official* sheet: Case study 8 (facility D) materials cost inflated to 2025 USD, 24 batches/yr, 22 doses/batch, 50% materials economy of scale, 4× COGS-to-retail markup → **$19,917/dose**, range $9,958–$29,876 | Gamma (α ≈ 15.4, scale ≈ 1,297) | **Derived and wired in 2026-08-04** (`R/04_costs_utilities.R`'s `treg_dose_cost()`) — arithmetic is traceable in the workbook; the markup and scale assumptions are the authors'. **New finding while wiring this in:** `data/processed/model_dose_costs_and_psa_ranges.csv`'s "TREG dose cost" `value_usd` column ($18,908.74) disagrees with this same $19,917 figure by a consistent ~0.9494× factor — and the identical ~0.9494× gap appears independently for that file's IFX dose cost row too ($1,014.20 vs $1,068), ruling out coincidence. `treg_dose_cost()` uses the traceable $19,916.75 figure (from `model_tenham_derived_treg_dose_cost.csv` directly), not `value_usd`, whose derivation isn't reconstructable from any material in this repository — comparable in kind to A6/A7 (`docs/model_audit_v6.md`), not resolved, just not silently picked around. |
| 13 | Treg administration | Infusion administration + monitoring. **Revise:** the reference programme administers a ~30-minute infusion with possible overnight observation, and the trial protocol includes a cyclophosphamide preconditioning component | Gamma | **Partially resolved 2026-08-04.** Infusion administration: the ~30-minute-vs-1-hour distinction doesn't change the billing fee (CPT/HCPCS 96365 covers "initial, up to 1 hour" either way) — reuses the $57.90 rate already sourced for UST/IFX/ADA. Cyclophosphamide preconditioning: **price** sourced (HCPCS J9075, $0.0962/mg, same CMS ASP-methodology fee schedule as IFX) but the **dose** is not — the cited reference trial (NCT06721962, "RESTORE") confirms "low dose cyclophosphamide conditioning" in its public registry without publishing a specific mg figure; `treg_dose_cost()` defaults `cyclophosphamide_dose_mg` to 0, an explicit placeholder, not a claim it's unneeded. **Overnight observation-stay cost: PRICE now fully sourced (2026-08-05), still not a base-case default — see why below.** CMS's Comprehensive Observation Services payment (C-APC 8011) national base rate (wage-index unadjusted) history: $2,174.14 (CY2016) → $2,349.66 (CY2018) → $2,283.16 (CY2021) → $2,439.02 (CY2023), from independent sources citing the CMS OPPS final rule. **The CY2026 figure is now a primary-source-confirmed fact**: CMS's own "January 2026 Web Addendum A" gives APC 8011 directly — RW 29.2310, national unadjusted Payment Rate **$2,672.15** — which exactly reproduces (to the cent) against a $91.415 conversion factor independently back-calculated the same day from the companion "2026 January Web Addendum B" file (29.2310 × $91.415 = $2,672.15), a strong two-file cross-validation. This is a sourced number now, not a range. It still isn't wired into the base case for two reasons unrelated to pricing confidence: (a) $2,672.15 is the *national unadjusted* rate — a specific facility's payment is wage-index-adjusted; (b) C-APC 8011 pays per *qualifying stay* (≥8 hours of observation, no same/next-day surgical procedure), and whether a Treg infusion actually triggers that is a clinical assumption this project hasn't made. `treg_dose_cost()` takes this as an explicit opt-in `observation_stay_cost_usd` argument (default 0) — a scenario run can now cite $2,672.15 directly, with full confidence, once someone makes the clinical call on (b). |
| 14 | UST / IFX / ADA acquisition | CMS ASP/NADAC, biosimilar-era pricing, single stated date (2026-08-05) | Gamma | **Re-extracted and rebased 2026-08-05** — biosimilar-inclusive pricing is now the base case, originator pricing retained as the S8 comparability-with-Aliyev scenario; see §4.2 |
| 15 | Drug administration | CMS Physician Fee Schedule, HCPCS per product billing guides ($119.36 in current model) | Gamma | Sourced |
| 16 | Health-state monitoring/management | Malone et al. via Aliyev, inflated to 2025 USD — see Section 8 | Gamma | **Reconciled** — see Section 8 |
| 17 | Surgery state (per-cycle, not a one-time episode) | Aliyev $884 (2017) → $1,152.29 (2025), same treatment as items 16/18. Applied identically across TREG/UST/IFX/CT (Aliyev Costs Assumption #1). HCUP colectomy cost dropped from the base case — see Appendix A, finding A9 resolution | Gamma | **Reconciled 2026-08-04** — residual gap: $884 does not cleanly reconstruct from Suppl. Table 2's $1,475 PMPM; same open item as 16/18, not new |
| 18 | CT cycle cost | Aliyev $67 (2017) → ~$87.5 (2025) **plus** health-state costs | Gamma | **Fix required** (Appendix A, finding 9) |
| **Utilities** |
| 19 | Remission utility 0.82 | Aliyev / NICE TA352 (vedolizumab GEMINI) | Beta (re-parameterised from the 0.66–0.98 range) | Sourced — **base case rebuilt on this value 2026-08-05** (A16, `docs/model_audit_v6.md`); the deterministic chain is now derived in code from this item's own PSA file rather than read from a separately maintained snapshot that had drifted to a stray PSA draw (0.9554) |
| 20 | Mild:Remission 0.89; M-SR:Mild 0.89; M-S:M-SR 0.87; Surgery:M-SR 0.87 | Buxton 2007 CDAI→EQ-5D mapping via Aliyev | Beta, correlated | Sourced |
| 21 | SDR utility | Base case = Remission utility; scenario = age-matched general population | Beta | **Assumption** |
| **VOI inputs** |
| 22 | Eligible population per year | Section 9.4 | Fixed with scenario range | Derived |
| 23 | Decision-relevant horizon | 10 years (recommended) | Scenario 5 / 10 / 15 yrs | Assumption |

Two changes from Aliyev's approach are deliberate and should be stated in Methods: (i) Dirichlet rather than independent marginals for transition rows, so that sampled rows sum to one; (ii) Beta rather than Uniform for utilities, with the Uniform retained as a scenario for comparability.

### 7.2 The cure fraction and relapse hazard: how to parameterise honestly

This is the crux of the paper, and it is where a reviewer will look hardest. Three published analogs exist. None is Crohn's disease efficacy data, and none should be presented as such.

| Analog | What is actually reported | What it supports | Caveats |
|---|---|---|---|
| **Autologous polyclonal Treg therapy in type 1 diabetes (PolTREG PTG-007)** | 54 patients from Phase I/II trials followed 7–12 years. A *proportion* remained insulin-independent for 18–24 months; a *subset* remained in clinical remission (defined as low exogenous insulin requirement) at 7–12 years. Best results reported in patients co-treated with anti-CD20. The product has been qualified by EMA/CHMP for a centralised marketing-authorisation application (March 2026) | The **optimistic** bound: durable, multi-year benefit from a Treg product is not merely theoretical | Denominators for the durable subset are **not published in peer-reviewed form**; the evidence is company communications and conference presentation. Autologous, polyclonal, different disease, paediatric-onset population, co-administered rituximab in the best-responding subgroup |
| **Treg therapy in solid-organ transplant tolerance** | Liver-transplant protocols achieved immunosuppression withdrawal in 7/10 patients; kidney-transplant protocols achieved 0/16 full withdrawal, with rejection in 7/16 | The **moderate** bound, and the more important qualitative lesson: tolerance induction is strongly tissue- and context-dependent | Cross-indication borrowing without mechanistic justification; the gut is neither liver nor kidney |
| **Tr1 cell therapy in refractory Crohn's disease (Ovasave / CATS1)** | Effect was time-limited, on the order of ~5 weeks | The **pessimistic** bound, and the only one in the correct disease: Treg effects in Crohn's may not be durable at all | Different product, different manufacturing, autologous, and the programme did not proceed |

**Recommendation — Decision 4: do not let any single analog anchor the base case.** Instead:

1. **Make the base case explicitly agnostic.** Set the base-case cure fraction by structured expert elicitation among the clinical co-authors and, if feasible, two external IBD clinicians, using a documented protocol (SHELF or a simplified variant). Present the three analogs above to elicitees as anchoring evidence. Report the elicited distribution, not just its mean. This is a recognised, citable method for exactly this situation and converts "we assumed 25%" into "we elicited, and here is the protocol and the spread."
2. **Report the headroom result as co-primary (Aim 4).** The minimum π\* required for cost-effectiveness at each price is computed from the model, not assumed into it. If π\* comes out at, say, 8%, the paper can say something genuinely useful even though the true π is unknown: it is a target the field can be measured against.
3. **Carry all three analogs as named scenarios**, labelled by their source, so the reader can locate their own prior on the spectrum:

| Scenario | Cure fraction π | Annual relapse hazard h from SDR | Anchoring logic |
|---|---|---|---|
| Optimistic ("durable-remission analog") | High end of elicited range | Low (long plateau) | Pattern seen in the T1D long-term cohort |
| Moderate ("tolerance-induction analog") | Mid range | Moderate | Liver/kidney transplant split suggests partial, variable tolerance |
| Pessimistic ("time-limited-effect analog") | Near zero | High | Prior Tr1 experience in refractory Crohn's |
| Null | π = 0 | n/a | Treg collapses to a one-time-cost UST-equivalent. **Include this** — it is the honest floor and a useful reference case |

Deliberately, no numeric values are proposed here. Assigning "π = 0.30 because PolTREG" would be exactly the fabrication the design is meant to avoid: the published sources give a qualitative pattern and an unpublished denominator, not an estimable fraction. **The numbers should come from the elicitation in step 1 and be reported as elicited quantities.** If co-authors prefer to skip elicitation, the fallback is a uniform prior across the full 0–1 range for π, reported as such, with the headroom analysis carrying the paper. That is less informative but entirely defensible; an invented point estimate is not.

4. **Sensitivity to the relapse-timing shape.** A constant annual hazard from SDR is the simplest assumption and probably the right base case. Note in Limitations that mixture-cure models in oncology HTA often use a declining hazard, and that a constant hazard is conservative over a lifetime horizon.

**Status, 2026-08-05 (`docs/treg-cd_decision_resolutions_2026-08-05.md` §3):** Decision 4 is now **final, not "revisit at Gate 3"** — elicitation was not run for this submission (the case against it is strategic, not a feasibility gap: a three-person elicitation performed by the paper's own authors would be the easiest thing in the manuscript to attack); the uniform-prior fallback in the paragraph above is the committed approach. Point 4's own characterisation of a constant hazard as "conservative" needed correcting: **zero relapse hazard is the single MOST favourable assumption available to the intervention, not a conservative one** — a bug in the framing that would have shipped as a stated conservative assumption. `R/03_cure_fraction_module.R`'s `duration_to_hazard()`/`hazard_to_duration()` re-parameterise h by its implied median SDR duration T (h = ln2/T, base case T = 10 years, swept T ∈ {2, 5, 10, 20, ∞}) so this is reasoned about in units a reviewer can evaluate directly, and `R/05_deterministic_results.R`'s `headroom_frontier_by_duration()` sweeps T jointly with price and π at the **lifetime** horizon, not in isolation — h barely matters at Aliyev's 6.15-year horizon but dominates over a lifetime, so a one-way sweep of either alone would have understated both. Full detail: README.md's Status section.

### 7.3 Cycle-length conversion — reconstruction required

Aliyev's Appendix S2 documents the DEALE conversion (constant hazard, exponential survival) from trial endpoints to 2-week probabilities, with a worked example: ustekinumab week-6 remission 0.349 → annual rate 3.720 → 2-week probability 0.133. The current workbook contains 8-week probabilities as **hardcoded values in the Input Parameters sheet**, with no visible intermediate calculation. They cannot be independently verified from the file.

**Task:** re-derive every transition probability from the published trial endpoints, in code, using the same DEALE method, at 8-week cycle length; then diff against the hardcoded values currently in the workbook. Any discrepancy above a stated tolerance must be explained before it is adopted. The re-derivation script becomes `R/00_derive_transition_probs.R` and its output becomes a versioned CSV in `data/processed/`. This is the single highest-value reproducibility task in the project and should not be skipped on the grounds that the values "looked right."

---

## 8. Cost reconciliation — resolved

**2026-08-04 revision.** The cost-*cycle-length* question this section originally left open (were Aliyev's per-cycle costs meant for an 8-week cycle or something else?) is dissolved by the same decision that resolved it for transition probabilities: this project now runs Aliyev's native 2-week cycle throughout (`R/00_derive_transition_probs.R`'s third-revision header comment), and Appendix S2 states these are "costs per cycle" in the same document whose transition-probability section runs on that same 2-week cycle. There is no separate cost-cycle convention to reconcile — Aliyev's $217/$91/$10/$67/$884 apply directly, per 2-week cycle, no conversion needed. What remains is dollar-year inflation only (2017 → 2025 USD), addressed below.

Separately, re-verifying the *arithmetic* this section previously presented turned up an error, corrected here: the "Derivation in Aliyev" column below does **not** actually reproduce Aliyev's reported 2017 USD figures when computed — e.g. the stated M-S formula ($374 − $123 + $111) equals $362, not $217. This isn't a one-off transcription slip: every cost line shows the same ~0.6× gap between the raw PMPM-derived value and Aliyev's reported per-cycle figure, including Mild and Remission, which the appendix describes as a *direct* copy of the PMPM total with no subtraction formula at all ($152 reported as $91; $17 reported as $10). That's a systematic step in Aliyev's own derivation (possibly a Malone et al. PMPM-to-episode convention this project doesn't have access to) that Appendix S2 doesn't document and that **isn't independently reproducible from the materials in this repository** — the same situation as the IFX induction gap noted in Gate 1's `DERIVATION_NOTES.md`. This doesn't block using the figures: exactly as for the transition probabilities, this project takes Aliyev's *published* per-cycle costs as the primary source rather than re-deriving them, and the currency-year inflation step below checks out correctly against those published figures regardless of how he arrived at them.

| Health state | Aliyev S2 (2017 USD, per 2-week cycle) | PMPM-formula result (does not reproduce Aliyev's figure — see above) | Workbook / Supplement 1 (2025 USD) | Inflation factor (2017→2025) |
|---|---|---|---|---|
| Moderate-Severe | $217 | $374 − $123 + $111 = $362 | $282.86 (Supplement 1 rounds to $283) | 1.3035 |
| M-SR | $217 | assumed equal to M-S | $282.86 | 1.3035 |
| Mild | $91 | $152 (direct PMPM total, no formula) | $118.62 (Supplement 1: $119) | 1.3035 |
| Remission | $10 | $17 (direct PMPM total) | $13.04 (Supplement 1: $13) | 1.3035 |
| Conventional therapy | $67 | $111 (direct PMPM pharmacy) | $87.53 | ~1.306 |
| Surgery | $884 | $1,475 (direct PMPM total) | **$1,152.29** — Aliyev's own per-cycle figure, inflated by 1.3035 like every other row (was $35,518, an HCUP colectomy episode cost, replaced rather than inflated; resolved 2026-08-04, see Appendix A finding A9) | 1.3035 |

**The currency-year reconciliation holds**: the workbook/Supplement 1's 2025 USD figures are Aliyev's own 2017 USD figures (left column) inflated by ≈1.3035, consistently across every line — that part of the original reconciliation was correct and still stands.

**Three actions required:**

1. **Document the inflation index.** The factor of ~1.3035 must be named — which index, which base and end months. A 2017→2025 factor of 1.30 is plausible for a medical-care CPI or PCE health price index, but "plausible" is not a citation. This is CHEERS 2022 item 15 and reviewers check it.
2. **Fix the dollar-year inconsistency in the Treg arm.** The Treg arm currently uses **un-inflated 2017 values** — `costs_treg_MSR` = $217, `costs_treg_mild` = $91, `costs_treg_remission` = $10 — while UST and IFX use 2025 values. The Treg arm's non-drug maintenance costs are therefore understated by roughly 23%. (Note also `costs_treg_MS` = `=217+67` = $284, which is Aliyev's M-S plus CT cost in 2017 dollars and coincidentally lands near the correctly-inflated $282.86; it is not used by the Markov engine.)
3. **Adopt 2025 USD throughout** and state the price date once, in Methods. Recommend re-checking whether a 2026 price year is more appropriate given the likely publication date.

---

## 9. Outcome and analysis plan

### 9.1 Economically justifiable price

**Definition.** The EJP is the per-dose acquisition price P\* at which the incremental net monetary benefit of Treg versus the next-best non-dominated comparator equals zero at threshold λ:

INMB(P\*) = λ · (Q_treg − Q_comp) − (C_treg(P\*) − C_comp) = 0

Total Treg cost is linear in price: C_treg(P) = C_treg,non-drug + D̄ · P, where D̄ is the discounted expected number of doses per patient (D̄ = 1 under the recommended single-dose base case; D̄ < 2 under the two-dose scenario, because the second dose is given only to patients still on therapy and is discounted). Hence:

**P\* = [ C_comp + λ·(Q_treg − Q_comp) − C_treg,non-drug ] / D̄**

Three specification points that the current implementation gets wrong or leaves ambiguous:

- **Price, not total cost.** The current manuscript reports EJP as "a total cost of $95,382" and the workbook's EJP sheet treats the x-axis variable as total cost. Report the **per-dose price**, and report total cost separately if desired. A developer prices a dose.
- **The comparator must be the next-best non-dominated option**, identified on the efficiency frontier at the given λ — not fixed to infliximab by assumption.
- **Q_treg must be re-evaluated at P\***, which is trivial here because QALYs do not depend on price, but the identity should be stated so that a reader can verify the algebra.

**Probabilistic EJP (new).** For each PSA draw k, solve for P\*_k and report the median with a 95% credible interval, plus the full density. Also report the distinct quantity **P_50 = the price at which the probability of cost-effectiveness equals 50%**, and explain that these coincide only under symmetry. Reporting a single deterministic EJP for a technology whose efficacy is entirely hypothetical would be the most obviously criticisable choice in the paper.

**Presentation.** A price–probability curve (x = price, y = P(cost-effective)) at each λ, plus the (π, price) frontier from Aim 4, as two panels of one figure. This figure is the paper's centrepiece.

**Gross margin over COGS at P\* (added 2026-08-05; implemented the same day once `R/08_ejp.R` landed, `gross_margin_over_cogs()`).** Report gross margin over the ten Ham-derived manufacturing cost as a companion figure: (P\* − $4,979.19) / P\*, using `model_tenham_derived_treg_dose_cost.csv`'s "Total manufacturing (COGS) cost per treatment/dose" line as the cost basis (the same COGS `treg_dose_cost()` already marks up 4× to get the model's $19,916.75 acquisition-cost assumption — this is a different, EJP-anchored markup, not that one). **What this is:** an answer to "is the price this study says is cost-effective also a viable price for a manufacturer, at least at the COGS level" — a real, computable, useful number. **What this is NOT, and must not be presented as:** a claim about Tr1X's or any manufacturer's actual profitability. COGS margin excludes R&D recoupment (including the cost of failed programmes elsewhere in a manufacturer's portfolio), commercial infrastructure, cost of capital, and time-to-market risk — a full profitability assessment needs private company financials this project doesn't have and can't get. State this limitation explicitly alongside the figure, not just in a footnote.

For context when this figure is written up, `data/raw/market_comparator_cell_therapy_prices.csv` (added 2026-08-05) has two real allogeneic-cell-therapy pricing data points found while researching this: **Ryoncil** (remestemcel-L, Mesoblast's FDA-approved allogeneic MSC therapy for pediatric SR-aGvHD) has an actual US list price of $194,000/infusion (~$1.55M for the full 8-infusion course); **tabelecleucel/Ebvallo** (Atara's allogeneic EBV-specific T-cell therapy, EU/UK-approved, not FDA-approved) was rated cost-effective by ICER **if priced between $143,900 and $273,700 per treatment cycle** — the same kind of value-based pricing exercise as this section, for a different allogeneic cell product. Both are 8–14× the model's $19,917 acquisition-cost assumption. Neither is Treg/Crohn's-specific, and neither should be substituted for the ten Ham-derived cost parameter itself (one is a company list price, the other an HTA-recommended price — different quantities from a manufacturing-cost build-up) — but both are directly relevant Discussion-section context for how this study's own P\* and margin-over-COGS figures compare to what the market has actually shown allogeneic cell therapies can be priced at.

### 9.2 EVPI

**Definition.** Per-patient EVPI = E_θ[max_j NB_j(θ)] − max_j E_θ[NB_j(θ)], estimated from the PSA as the mean opportunity loss of the decision made under current information.

**Specification points:**

- EVPI here is **conditional on the Treg price**, because the decision set itself changes with price. Report EVPI as a *surface* over price × λ, not a single number. This is a genuinely novel and useful output for a pre-development technology and is worth a figure of its own.
- Frame the construct honestly in Methods: this is not a conventional payer-adoption EVPI, since the technology does not exist. It is the value of resolving uncertainty for a decision-maker choosing among these options *conditional on the technology being available at price P* — i.e. a prospective, developer- and funder-facing VOI. Say so explicitly rather than letting a methodologically literate reviewer say it first.

### 9.3 EVPPI

**Parameter subsets.** Minimum set, each chosen to map onto a real, fundable research action:

| Subset | Parameters | The research question it prices |
|---|---|---|
| A | Cure fraction π | Is it worth running a larger Phase 2 to estimate durable remission rate? |
| B | Post-cure relapse hazard h | Is it worth extending follow-up of existing early-phase cohorts? |
| A ∪ B | Durability block | Combined value of the whole durability question (report jointly — A and B are decision-relevant together) |
| C | Treg acquisition/manufacturing cost | Is further manufacturing cost-engineering worth funding? |
| D | UST/IFX/ADA maintenance transition probabilities | Would a comparator head-to-head trial change the decision? |
| E | Health-state utilities | Is a CD-specific utility study warranted? |
| F | Induction response probabilities (all arms) | Does short-term response need better estimation? |
| G | Surgery transitions and surgery cost | Does the surgery pathway need better data? |

**Method.** Estimate with the `voi` R package: nonparametric GAM regression (Strong, Oakley & Brennan 2014) for low-dimensional subsets (A, B, C, and the A∪B pair); the SPDE-INLA method (Heath, Manolopoulou & Baio) for the higher-dimensional blocks (D, E, F). Cross-check at least two subsets with `BCEA::evppi()` and report agreement. Where subsets are dimension-heavy, consider a principal-components pre-reduction and report it.

**Sample size.** Use **10,000 PSA draws minimum**, not the current 1,000. EVPPI estimates from 1,000 draws are unstable, particularly for regression-based estimators, and a reviewer familiar with the method will ask. Report a convergence check (EVPPI estimate vs. number of draws).

**Statistical hygiene note.** The current workbook computes an ICER per PSA iteration and appears to summarise across them. Ratios must not be averaged. All probabilistic summaries should be computed on the **incremental net monetary benefit** scale.

**π prior-sensitivity, added 2026-08-05 (`docs/treg-cd_decision_resolutions_2026-08-05.md` §3.2).** Subset A's own PSA prior, π ~ U(0,1) (Decision 4's fallback), maximises prior variance on π — under it, π is liable to come out as the dominant EVPPI parameter "more or less by construction," which a reviewer could (correctly) read as an artefact of the prior rather than a finding. `R/07_evpi_evppi.R`'s `evppi_prior_sensitivity()` re-estimates the full subset table under two named alternative priors, Beta(1,3) (mass toward low cure fractions, per the Ovasave/CATS1 analog in §7.2) and Beta(2,2) (symmetric, informative), via `prior_reweight()`'s importance-weighting of the *same* PSA draws — no re-simulation. The result to report in Methods: the **ranking** of subsets by EVPPI is the finding; the **level** is prior-dependent, stated as such rather than presented as a single number.

### 9.4 Population-level VOI — Decision 6

Population EVPI = per-patient EVPI × effective population size, where:

Effective population = Σ over the decision horizon of (annual eligible incident + initially eligible prevalent patients), discounted at 3%.

**Building blocks and their status:**

| Input | Value | Source | Status |
|---|---|---|---|
| US Crohn's disease prevalence | 305 per 100,000; ≈1.011 million Americans | Lewis et al., *Gastroenterology* 2023 | **Sourced** |
| US IBD incidence | 10.9 per 100,000 person-years (high-probability algorithm; 15.9 including lower-probability) | Lewis et al. 2023 | **Sourced** — CD-specific incidence must be extracted from the paper's age/sex tables rather than assumed as a share of IBD |
| Proportion with moderate-to-severe disease eligible for advanced therapy | — | To be sourced | **Open — must not be assumed.** Sensitivity across a wide range |
| Proportion biologic-naïve (base case) vs. refractory (scenario) | — | To be sourced | **Open** |
| Decision horizon | **10 years recommended** (scenarios 5 and 15) | Assumption | Rationale: the reference Phase 1/2a completes around 2027; a Phase 3 programme plus review plausibly places a real adoption decision in the early-to-mid 2030s, after which competing modalities would likely supersede the specific decision modelled |
| Discount rate on future research value | 3% | Consistent with the model | — |

**Recommendation.** Report population VOI for the *incident + prevalent eligible* population under the base-case (biologic-naïve) definition, and separately for the refractory-eligible population, since the second is much smaller and gives a much more conservative — and for a first-indication technology, arguably more decision-relevant — figure. Present both; do not choose one silently. Show the arithmetic in a supplementary table so a reader can substitute their own denominators.

---

## 10. Sensitivity and scenario analysis plan

### 10.1 Deterministic

One-way analysis across all parameters in Section 7.1 over their credible ranges, presented as a tornado on incremental NMB (Treg vs. next-best comparator). Port the existing 9-parameter Excel tornado to R and extend it to the full parameter set. **Also run two-way analysis on (π, price)** — this generates the Aim 4 frontier and is the more informative display for this paper than the tornado.

### 10.2 Probabilistic

10,000 draws. Outputs: cost-effectiveness plane, CEACs across λ ∈ [$0, $300,000], cost-effectiveness acceptability frontier, and the EJP posterior. Correlations must be preserved where they exist — in particular, the utility ratios are structurally correlated (each state's utility is derived multiplicatively from the Remission utility) and must be sampled as a chain, as the current workbook correctly does.

### 10.3 Structural scenarios

| # | Scenario | Why it must be run |
|---|---|---|
| S1 | 2-year maintenance cap: on vs. off | Decision 1; materially moves the result |
| S2 | ADA included vs. excluded | Decision 2 |
| S3 | Biologic-naïve vs. refractory population | Decision 3 — **implemented 2026-08-05, §5.4**; response/remission multipliers only, run at the 6.15-year horizon (`run_refractory_scenario()`, `output/tables/refractory_scenario_results.csv`); surgery-hazard elevation not yet sourced |
| S4 | Cure-fraction anchor: optimistic / moderate / pessimistic / null | Decision 4 |
| S5 | Time horizon: lifetime / 10-year / 6.15-year | Tests the truncation critique directly and reproduces the current draft — **lifetime arm implemented 2026-08-05, §4.3**; deterministic base case and headroom frontier only (PSA/EVPI/EVPPI/probabilistic EJP still 6.15-year) |
| S6 | Treg dosing: 1 dose vs. 2 doses | Current draft assumes 2; the reference clinical programme is a single infusion |
| S7 | SDR utility = Remission vs. general-population | Locates how much of the cure's value is utility-driven |
| S8 | Comparator pricing: **biosimilar-era pricing (base case, as of 2026-08-05) vs. originator/pre-biosimilar pricing (now the scenario)** | Section 4.2 — inverted from the original S8 framing; `load_drug_prices(pricing_basis = "originator_pre_biosimilar")` runs this scenario |
| S9 | Perspective: healthcare sector vs. societal | Second Panel; favours durable therapy |
| S10 | Discount rate: 3% / 0% / 1.5% / 5% | Discounting of one-time curative therapies is contested |
| S11 | Relapse from SDR re-enters at Mild vs. M-SR | Minor; cheap to run |
| S12 | Non-cured Treg = UST-equivalent vs. HR-advantaged | Tests the Section 6.2 recommendation |

Twelve structural scenarios is more than can be reported in 6,000 words. Recommend S1, S3, S4, S5 and S8 in the main text; the remainder in Electronic Supplementary Material with a summary table in the main text.

---

## 11. CHEERS 2022 compliance plan

The journal requires CHEERS 2022 for economic evaluations (Husereau et al., *PharmacoEconomics* 2022, doi:10.1007/s40273-021-01112-8). The completed checklist should be submitted as an Online Resource.

The full 28-item mapping, and the two additional journal requirements (LLM-use disclosure; the Declarations section), have been moved to **`docs/CHEERS_2022_checklist.md`** so it can be maintained as a fillable submission artifact independent of this plan. Update the "Manuscript location" column there as the manuscript is drafted.

---

## 12. Software, tools and reproducibility

### 12.1 Recommended stack

| Purpose | Package | Note |
|---|---|---|
| Markov engine | **Custom, matrix-based** | Recommended over `heemod`. The week-56 landmark cure split and the arm-specific CT-switch rules are awkward to express in `heemod`'s formula interface, and a hand-written engine (a list of transition matrices multiplied through a state-occupancy vector) is ~150 lines, fully auditable and much easier to embed in a PSA loop. Use `heemod` only to cross-validate a simplified version |
| VOI | `voi` | Primary. `evppi()` with `method = "gam"` and `method = "inla"` |
| VOI cross-check | `BCEA` | `evppi()`; also gives publication-ready CEAC/CEAF plots |
| CEA visualisation | `dampack` | CE plane, CEAC, tornado, one-way and two-way analysis |
| Data handling | `tidyverse` | |
| Distribution fitting | `fitdistrplus`, `SHELF` | `SHELF` if expert elicitation is run — it also handles the elicitation session itself |
| Environment locking | `renv` | Required for the DAS |
| Reporting | `quarto` | Generate all tables/figures from code; no hand-copied numbers into the manuscript |
| Testing | `testthat` | At minimum: transition rows sum to 1; cohort size conserved; QALY/cost discounting matches a closed-form check; known-input regression tests |

### 12.2 Repository structure

```
treg-crohns-cea/
├── README.md                     # project summary, how to reproduce, citation, licence
├── LICENSE                       # MIT for code; CC-BY-4.0 for data/docs
├── CITATION.cff
├── renv.lock
├── data/
│   ├── raw/                      # verbatim source extracts, one file per source
│   │   ├── aliyev2019_tableS2_parameters.csv
│   │   ├── aliyev2019_tableS3_induction.csv
│   │   ├── aliyev2019_tableS4_maintenance.csv
│   │   ├── tenham2020_manufacturing_costs.csv
│   │   ├── cms_asp_<date>.csv
│   │   ├── hcup_surgery_costs.csv
│   │   └── cure_fraction_elicitation.csv
│   ├── processed/                # model-ready, generated by scripts — never hand-edited
│   └── data_dictionary.md        # every parameter: definition, unit, distribution, source, dollar year
├── R/
│   ├── 00_derive_transition_probs.R   # DEALE re-derivation (Section 7.3) — run first
│   ├── 01_decision_tree.R
│   ├── 02_markov_engine.R
│   ├── 03_cure_fraction_module.R
│   ├── 04_costs_utilities.R
│   ├── 05_deterministic_results.R
│   ├── 06_psa.R
│   ├── 07_evpi_evppi.R
│   ├── 08_ejp.R                  # deterministic + probabilistic EJP, headroom frontier
│   └── utils/
├── analysis/
│   ├── run_base_case.R
│   ├── run_scenario_analyses.R
│   └── run_full_analysis.R       # single entry point, reproduces every number in the paper
├── tests/testthat/
├── output/
│   ├── figures/
│   └── tables/
└── docs/
    ├── analysis_plan.md          # this document, versioned
    ├── model_structure.md
    ├── model_audit_v6.md         # Appendix A, preserved as the record of what changed and why
    └── CHEERS_2022_checklist.md
```

**Two conventions worth enforcing from day one:** (i) `data/processed/` is generated, never edited — every number in it traces to a script and a file in `data/raw/`; (ii) `analysis/run_full_analysis.R` regenerates every figure and table in the manuscript from a clean checkout. If a number appears in the paper that cannot be produced by that script, it does not go in the paper. That rule is what prevents a recurrence of the disconnected EJP worksheet described in Appendix A.

**Note on the tree above:** the `data/raw/` and `data/processed/` filenames shown are illustrative of the target naming scheme. The files actually populated in this repository during the initial data pull use the names and split documented in `data/data_dictionary.md`, which is authoritative — e.g. `aliyev2019_appendixS1_table1_assumptions.csv` and `aliyev2019_appendixS1_table2_parameters.csv` rather than the `tableS2`/`tableS3`/`tableS4` names shown here, and no `cure_fraction_elicitation.csv` yet, since the elicitation in Section 7.2 has not been run. `R/`, `analysis/`, `tests/testthat/` and `output/` contain scaffolding stubs only — no analysis code has been written yet (Gate 0/1 in Section 14 have not closed).

### 12.3 Archiving and the Data Availability Statement

The journal mandates a Data Availability Statement for original research and strongly encourages repository deposition. A GitHub URL alone is not a persistent identifier. **Recommendation:** link the repository to Zenodo, tag a release at submission (`v1.0-submission`), and cite the resulting DOI in the DAS and reference list, per the journal's data-citation guidance (Creator, Title, Publisher, Year, DOI). Re-tag at acceptance if the code changes during review.

Draft DAS:

> All model code, input parameter files and analysis scripts required to reproduce every result in this article are openly available at https://github.com/<org>/treg-crohns-cea and archived at Zenodo (DOI: 10.5281/zenodo.XXXXXXX). No individual participant data were used; all inputs are published aggregate estimates, cited in Table X and the Electronic Supplementary Material.

Also submit the parameter tables as `.csv` or `.xlsx` Online Resources, per the journal's supplementary-file format requirements.

---

## 13. Anticipated limitations and likely reviewer critiques

Stated proactively; each should have a corresponding sentence in the manuscript's Discussion, and several should be pre-empted in Methods.

| Critique | Expected form | Our response |
|---|---|---|
| **"You invented the efficacy."** | The most likely rejection reason. A reviewer sees a CEA of a technology with zero clinical data | The paper's primary outputs are the price and cure fraction *required* for cost-effectiveness, plus VOI — none of which requires assuming efficacy. Structure the abstract and Results so this is unmissable in the first screen. The null scenario (π = 0) is reported |
| **Cross-disease/indication borrowing** | "Type 1 diabetes and kidney transplant tell you nothing about Crohn's" | Correct, and we say so. The analogs anchor the *range* and the elicitation, not the point estimate. The one Crohn's-specific analog is the pessimistic one, and it is reported |
| **Population mismatch** | The reference clinical programme is refractory; the base case is biologic-naïve | Both populations run; the mismatch is stated in Methods, not buried in Limitations |
| **Comparator prices were stale / pre-biosimilar** | Ustekinumab pricing has moved | **Resolved 2026-08-05** — biosimilar-era pricing re-extracted at a stated date and made the base case; originator pricing retained as scenario S8 (§4.2) |
| **Structural uncertainty is not captured by PSA** | Whether a cure exists at all is a structural question, and PSA cannot price it | Reported as a distinct uncertainty category (CHEERS item 20), with π = 0 as the structural bound |
| **Constant relapse hazard from SDR** | Oncology mixture-cure models often use declining hazards | Acknowledged; constant hazard is conservative over a lifetime horizon; scenario available |
| **Lifetime extrapolation from ≤1 year of comparator trial data** | Standard critique of any long-horizon CD model | Report 6.15-year and 10-year scenarios; state the extrapolation assumption explicitly |
| **DEALE assumes constant hazards** | Inherited from Aliyev | Acknowledge; note it applies equally to all arms |
| **Utilities are mapped, not measured** | Buxton 2007 CDAI→EQ-5D algorithm | Inherited from Aliyev and NICE TA352/TA456; scenario on utility source |
| **VOI for a non-existent technology is unconventional** | A methodologically sophisticated reviewer may challenge the construct | Frame it up front as conditional/prospective VOI for developers and funders, not payer-adoption EVPI |
| **No adverse events modelled for a cell therapy** | Allogeneic cell products carry immunogenicity, infusion-reaction and (where preconditioning is used) cytopenia risks | **Currently a real gap.** At minimum, model administration-related costs including any preconditioning and observation stay, and run a scenario with a disutility and cost for serious adverse events. Reviewers of a cell-therapy CEA will expect this |
| **Manufacturing cost is a single-source academic costing** | ten Ham et al. is a small-scale academic costing, and the 4× markup is the authors' assumption | Already varied ±50%; make the markup assumption explicit in the main text and report EJP, which sidesteps the cost assumption entirely for the paper's primary conclusion |

The last row points at the strongest defensive framing available: **the EJP result does not depend on the manufacturing cost estimate at all.** Lead with it.

---

## 14. Timeline and task allocation

Sequenced so that structural decisions precede parameterisation, and parameterisation precedes VOI coding. Durations are working-week estimates for part-time academic effort; role assignments are proposals for the co-authors to adjust.

| Phase | Weeks | Tasks | Lead | Gate |
|---|---|---|---|---|
| **0. Sign-off** | 1 | Review this document; resolve Decisions 1–6 in a single meeting; record decisions in `docs/analysis_plan.md` | All three | **Gate 0: no coding before this closes — closed 2026-08-04 (Decisions 1, 4, 5 final; 2, 3, 6 provisional, see §15)** |
| **1. Foundations** | 2–3 | Set up repository, `renv`, testing harness. Re-derive all transition probabilities from published endpoints via DEALE (§7.3) and diff against v6. Resolve and document the inflation index (§8). Re-extract comparator ASP at a stated date | Model lead (Jadambaa) | **Gate 1: re-derived probabilities reconcile with v6, or discrepancies are explained in writing** — **scope narrowed 2026-08-04 (E. Stone), then closed 2026-08-04 (second revision):** an initial from-scratch DEALE re-derivation attempt was retired after Aliyev et al. 2019 Appendix S2 (the paper's own transition-probability methodology appendix) became available — it showed Aliyev's method for IFX/ADA requires each trial's separate placebo-arm endpoint, not present in `data/raw/aliyev2019_appendixS1_table2_parameters.csv`, so from-scratch re-derivation from that file alone was not achievable. Appendix S2's own Supplementary Tables 3 (induction) and 4 (maintenance) were transcribed and verified by direct visual comparison against the appendix's table images, added to `data/raw/`, and are now used directly (validated for internal consistency; maintenance cycle-converted 2-week → 8-week via exact Markov-chain matrix power, not DEALE). All four therapies (UST/IFX/ADA/CT), both phases — no remaining data gaps. See `R/00_derive_transition_probs.R` and `data/processed/DERIVATION_NOTES.md`. |
| **2. Engine** | 3–5 | Build decision tree, Markov engine, cure module, cost/utility module. Reproduce the v6 base case as a regression test before adding the cure structure — then reproduce Aliyev's published results as external validation | Model lead | **Gate 2: v6 and Aliyev both reproduced** — **note 2026-08-04:** "reproduce v6" is superseded by the Gate 1 decision to source Aliyev's own published matrices directly rather than the v6 workbook's derivation; Gate 2's live criterion is reproducing Aliyev's own published results (Aim 5), not the workbook. `R/02_markov_engine.R` (cohort Markov core: `simulate_cohort()`, the two-track `run_maintenance_arm()` implementing the M-S-to-CT switch and §6.4's 2-year cap, `discount_factor()`), `R/01_decision_tree.R` (induction split into the biologic/CT initial occupancy vectors `run_maintenance_arm()` needs), and `R/03_cure_fraction_module.R` (§6.2's mixture-cure extension: week-56 landmark — cycle 28 at the native 2-week cycle — split into Sustained Deep Remission, relapse hazard, and cap-aware relapse re-entry into CT once the 2-year cap has fired, confirmed with E. Stone 2026-08-04 and not otherwise specified in this plan text) are built, wired together, and tested end to end for UST; `R/04` (costs/utilities, all four arms) and `R/05` (deterministic results, first pass — base case plus the Aim 4 headroom frontier) are also now built (README.md's Status section is the authoritative, current record — this cell is not kept in sync line-by-line). **Aim 5 external validation is a real, open scope gap, not just an unstarted task:** this project transcribed Aliyev's transition-probability appendix only, not his own cost/utility parameters or a results table, and deliberately re-sources costs from current CMS pricing rather than his 2017 figures (R/04's module header) — so there is nothing of his to numerically diff the full costed model against from data already in this repo, only his input matrices (already confirmed used unmodified). Closing Gate 2 on Aim 5's own terms would need either transcribing more of his paper or restating the criterion as validating the epidemiological (undiscounted, cost-free) trace instead — an open decision, not resolved here. |
| **3. Parameterisation** | 4–6 (parallel) | Cure-fraction elicitation protocol, session, and fitting. Source refractory-population multipliers. Source eligible-population denominators for VOI. Adverse-event and administration cost inputs | Clinical lead (Abraham) with senior author (Stone) | **Gate 3: parameter table frozen and reviewed** |
| **4. Deterministic results** | 6–7 | Base case, all structural scenarios, tornado, two-way (π, price) frontier, EJP | Model lead | `R/05_deterministic_results.R` first pass landed 2026-08-05 (base case, Aim 4 headroom frontier); `R/08_ejp.R` first pass landed the same day (deterministic + probabilistic EJP, gross margin over COGS) — tornado and the other structural scenarios still open, see README.md's Status section |
| **5. PSA and VOI** | 7–9 | 10,000-draw PSA; CEACs; probabilistic EJP; EVPI surface; EVPPI by subset with convergence checks and `BCEA` cross-check | Model lead | **Gate 4: EVPPI stable across draw counts** — `R/06_psa.R` first pass landed 2026-08-05 (10,000-draw PSA over the utility chain, Treg price, and π; CEAC and cost-effectiveness-plane outputs); `R/07_evpi_evppi.R` first pass landed the same day (EVPI price×λ surface, EVPPI by subset A/C/E with a convergence check) — probabilistic EJP still needs `R/08`; the `voi`/`BCEA` cross-check this row asks for could not be run (environment lacks a Fortran toolchain — see README.md's Status section for the substitute method used and why it's methodologically equivalent for the subsets actually computable so far) |
| **6. Writing** | 9–12 | Draft to journal structure with decimal headings; CHEERS checklist; figures at journal specification; Declarations; DAS | Senior author (Stone) drafting, all revising | |
| **7. Submission prep** | 12–13 | Zenodo release and DOI; ESM assembly; internal read-through against CHEERS line by line; cover letter | All three | |

**Timing consideration worth acting on.** Initial safety and efficacy results from the ongoing allogeneic Tr1 Crohn's programme have been publicly anticipated for late 2026. If they emerge during preparation, they will not supply usable efficacy parameters for a biologic-naïve cohort — but they will change what a reviewer expects the paper to cite, and a paper submitted without acknowledging a public readout will look dated. Build in a check at Gate 4 and again at submission. This argues for a compressed timeline rather than a leisurely one.

---

## 15. Open decisions requiring co-author sign-off

Each requires an explicit yes/no or selection. Record the outcome, with date, in `docs/analysis_plan.md`.

**☑ Decision 1 — Reinstate the 2-year biologic-maintenance-to-CT cap?**
*Recommendation:* **Yes for UST/IFX/ADA in the base case; never for cured (SDR) Treg patients; yes for non-cured Treg responders. Report the no-cap structural scenario alongside.**
*Rationale:* we are borrowing Aliyev's transition probabilities, and they were estimated for a structure containing that cap; using them without it is an unvalidated extrapolation. Be aware this moves the result against Treg. *Would change if:* a sourced real-world persistence curve is available, which would be better than either binary choice.
*Recorded 2026-08-04 (E. Stone):* **Yes — adopted as recommended.** Cap applies to UST/IFX/ADA and non-cured Treg responders after cycle 13; never applies to cured (SDR) Treg patients. No-cap scenario to be retained alongside the base case. Final.

**☑ Decision 2 — Add adalimumab as a third biologic comparator?**
*Recommendation:* **Yes.**
*Rationale:* Aliyev supplies the parameters at no derivation cost, it restores the original three-way comparison, and omitting the cheapest widely-used biologic from a cost-effectiveness comparison invites the criticism that the comparator set was chosen to flatter the intervention. *Would change if:* word-count pressure becomes acute — but ADA belongs in the model even if it is reported only in supplementary material.
*Recorded 2026-08-04 (E. Stone), provisional pending full co-author sign-off:* **Yes — adopted as recommended.** ADA included as a comparator arm, reported in supplementary material if main-text word count is tight.

**☑ Decision 3 — Base-case population: biologic-naïve, refractory, or both?**
*Recommendation:* **Both, with biologic-naïve as base case and refractory as a co-primary scenario reported in the main text.**
*Rationale:* biologic-naïve preserves comparability with Aliyev and NICE TA456 and keeps every comparator parameter valid; refractory is where the technology would actually first be used. *Would change if:* the co-authors decide the paper's primary audience is the developer rather than the HTA community, in which case invert the two.
*Recorded 2026-08-04 (E. Stone), provisional pending full co-author sign-off:* **Both — adopted as recommended.** Biologic-naïve is the base case; refractory reported as a co-primary scenario in the main text.

**☑ Decision 4 — What anchors the base-case cure fraction and relapse hazard?**
*Recommendation:* **None of the three analogs. Anchor the base case on documented structured expert elicitation, present all three analogs as named scenarios plus a null (π = 0) case, and report the required cure fraction (headroom) as a co-primary result.**
*Rationale:* the published analogs give a qualitative pattern and, in the optimistic case, an unpublished denominator — not an estimable fraction. Reporting the required π converts the largest unknown from an input into an output. *Would change if:* elicitation is judged infeasible, in which case fall back to an explicit uniform prior over π and let the headroom analysis carry the paper. Do not fall back to a point estimate.
*Recorded 2026-08-04 (E. Stone):* **Elicitation deferred for now — fall back to the documented alternative.** π is swept as an explicit uniform variable over the full 0–100% range (not a point estimate), with the headroom analysis (minimum π\* for cost-effectiveness at a given price) as the co-primary result per Aim 4. In addition to the continuous sweep, report named scenario points at π = 0% (null/floor case), 50%, 75%, and 90%, each carrying its own relapse hazard h consistent with §7.2's optimistic/moderate/pessimistic framing. Structured expert elicitation remains open to revisit if it becomes feasible before Gate 3. Final for now, revisit at Gate 3.
*Amended 2026-08-05 (`docs/treg-cd_decision_resolutions_2026-08-05.md` §3.1):* **Now final, not "revisit at Gate 3."** A three-person elicitation performed by this paper's own authors would be the single easiest thing in the manuscript to attack — the uniform-prior-plus-headroom-frontier position is stronger precisely because it concedes the unknown rather than dressing it. Structured elicitation is instead named in Discussion as the natural next step once first efficacy readouts exist. §3.2/§3.3's own amendments (prior-sensitivity reweighting on EVPPI; the h sweep parameterised by median SDR duration T rather than a bare hazard) are implemented — see §7.2 and §9.3's own status notes and README.md's Status section.

**☑ Decision 5 — Cost reconciliation and dollar year.**
*Recommendation:* **Adopt 2025 USD throughout** (or 2026, if publication timing warrants). The reconciliation is resolved in §8: Aliyev's 2017 per-cycle costs ($217 M-S and M-SR, $91 Mild, $10 Remission, $67 CT) inflated by ≈1.3035 give the workbook's 2025 figures ($282.86 / $118.62 / $13.04 / $87.53); surgery was replaced with an HCUP episode cost rather than inflated. **Two fixes required:** name and cite the inflation index, and correct the Treg arm, which is currently charged un-inflated 2017 health-state costs.
*Would change if:* a co-author can identify a different index actually used, in which case document that one.
*Recorded 2026-08-04 (E. Stone):* **2025 USD throughout — adopted as recommended**, per the reconciliation already resolved in §8. The two required fixes (name/cite the inflation index; correct the Treg arm's un-inflated 2017 health-state costs) carry forward as implementation tasks for the costs/utilities module (`R/04_costs_utilities.R`), not open decisions. Final.

*Addendum, 2026-08-04:* the surgery/HCUP substitution named in this decision's own rationale line ("surgery was replaced with an HCUP episode cost rather than inflated") was never actually revisited — it was recorded as accepted fact, not evaluated. It should not have been: Appendix A finding A9 already flagged it as needing an explicit decision, and it does not survive review. Aliyev's Surgery state is not a colectomy proxy (see A9 resolution in `docs/model_audit_v6.md`); the HCUP colectomy episode cost ($30,389–$35,518) is dropped from the base case and replaced with Aliyev's own per-cycle Surgery figure, inflated on the same basis as items 16/18 (→ $1,152.29, 2025 USD), applied identically across all four arms including Treg. Item 17 above and `data/processed/model_health_state_costs.csv` updated accordingly. This closes A9's "explicit decision" requirement as part of Decision 5 rather than as a separate decision, since it is the same cost-reconciliation exercise, not a new modelling choice.

*Addendum, 2026-08-05:* two further items from the same cost/utility-reconciliation family, resolved together in one branch per `docs/decision_resolutions_2026-08-05.md` §11 ("blocks everything downstream; every number moves" — the two highest-priority items in that memo, above everything else in Bucket 3):

- **A16 (Remission utility, Appendix A / `docs/model_audit_v6.md`) closed.** The deterministic utility base case was running on 0.9554396356 (Remission), traced to a stray live PSA draw captured in the workbook snapshot, not the 0.82 literature-cited base value this section's own item 19 names as the source (see `docs/model_audit_v6.md`'s A16 entry for the four independent lines of evidence). `R/04_costs_utilities.R`'s `load_health_state_utilities()` now derives the deterministic vector from `data/processed/model_psa_parameter_distributions.csv`'s sourced 0.82-based values instead of reading the retired snapshot file, so the deterministic base case is by construction the central draw of the PSA. **Every deterministic result produced before this fix (§10.1's base case, the EJP in §9.1) is superseded, not revised** — QALYs fall ~8–15% across every arm, ICERs and the required cure fraction π\* rise, and the EJP falls. `tests/testthat/test-parameter-provenance.R` now makes this class of defect (a `data/processed/` snapshot value used deterministically without a `data/raw/`-traceable source) fail loudly rather than requiring manual rediscovery.
- **A15 (dose-cost discrepancy, Appendix A) closed**, no fix required: the ~0.9494× factor between `model_dose_costs_and_psa_ranges.csv`'s `value_usd` and `psa_base` columns is the same class of snapshot artefact as A16 (a shared multiplicative shock, identical to six figures across two independently-derived drugs, that cannot have originated inside either build-up). `psa_base` — already the figure `R/04_costs_utilities.R` uses — is independently traceable end to end; `value_usd` is not and was never wired in. Downgraded from "unresolved" to "resolved" in `docs/model_audit_v6.md`; dropped from the manuscript limitations list entirely, since it was never a modelling choice to defend.
- **§4.2's comparator-pricing action item resolved the same day, in the same branch**, on the argument that it belongs beside A16 as a co-priority rather than under ordinary prioritisation: at π=0 the EJP collapses to essentially the discounted comparator drug spend displaced, so comparator acquisition price is the single largest lever on every headline number. See §4.2 for the full re-pricing (UST/IFX/ADA re-extracted at a single stated date, biosimilar-inclusive pricing now the base case, originator pricing the S8 scenario) and the resolution of the Part B ASP vs. NADAC benefit-routing question for UST maintenance specifically.

*Addendum, 2026-08-05:* item 4 of the executive summary's own redesign list (§0.1) — extending the time horizon to lifetime — is now implemented, not just recommended. See §4.3 for the full detail (life-table sourcing, why background mortality replaces rather than adds to Aliyev's own embedded trial mortality, the male/female sub-cohort mechanics) and §7.1 item 7's status update. This closes the "structurally cannot capture the value of a cure" problem the 6-year horizon named as the reason for the redesign item in the first place: at the 6.15-year horizon, no cure fraction at any plausible price made Treg cost-effective; at the lifetime horizon, the required cure fraction π\* at Treg's sourced price is feasible at all three WTP thresholds (§10.1/§10.3 S5, `output/tables/headroom_at_sourced_price_lifetime.csv`). PSA/EVPI/EVPPI/probabilistic EJP remain at the 6.15-year horizon this pass (§4.3's closing paragraph) — not yet a decision requiring sign-off, a deferred performance question.

**☑ Decision 6 — Population and horizon for monetised EVPI/EVPPI.**
*Recommendation:* **Incident plus prevalent eligible US moderate-to-severe CD population, over a 10-year decision horizon, discounted at 3%; reported separately for the biologic-naïve-eligible and refractory-eligible denominators.**
*Rationale:* prevalence (≈1.011 million with CD; 305 per 100,000) and IBD incidence (10.9 per 100,000 person-years) are sourced from Lewis et al. 2023, but the moderate-to-severe and treatment-line fractions are not yet sourced and must not be assumed. A 10-year horizon reflects the plausible interval to a real adoption decision. *Would change if:* a co-author prefers a shorter horizon on the grounds that the modality will be superseded faster — run it as a scenario either way.
*Recorded 2026-08-04 (E. Stone), provisional pending full co-author sign-off:* **10-year horizon at 3% discount, incident + prevalent eligible population — adopted as recommended.** 5- and 15-year horizons to be run as scenarios per the original recommendation. Moderate-to-severe and treatment-line fractions still require sourcing before this can be monetised (not resolved by this decision).

**One further item for the record, not a decision but a disclosure:** agree now how LLM assistance in model development, code or drafting will be described in the Methods section, per the journal's authorship policy.


---

## Appendix A. Audit of the IBD CEA v6 workbooks

Findings from a direct read of `IBD_CEA_v6_PSA.xlsm` and `IBD_CEA_v6_Univariate_Sensitivity_Analysis.xlsm` (formulas and cached values), summarised in Section 0.2 above. Fourteen defects or undocumented choices were identified (A1–A14), several material and at least three biasing results in favour of Treg (A1, A3, A4), one against Treg (A2), and one favouring comparators (A8).

The full findings, with cell references and effect direction for each, have been moved to **`docs/model_audit_v6.md`** so they can be preserved as a standalone record of what changed and why, independent of this plan.

---

## Appendix B. Provenance of external claims used in this design

Explicitly separating what is sourced from what is assumed, per the working instructions.

### B.1 Verified by search, August 2026

| Claim | Status |
|---|---|
| *PharmacoEconomics* Original Research Article: up to 6,000 words plus unlimited ESM; no limit on tables, figures or references; word count excludes abstract, references, figure legends and table captions | Confirmed on the journal's submission-guidelines page |
| Economic evaluations must follow the CHEERS checklist (linked to doi:10.1007/s40273-021-01112-8) | Confirmed |
| Abstract 150–250 words, extendable to 450 for reporting-guideline compliance | Confirmed |
| Decimal headings, maximum three levels | Confirmed |
| "Declarations" section before the reference list (Funding, Competing interests, Ethics approval, Consent, Data/Material/Code availability, Author contributions) | Confirmed |
| Data Availability Statement mandatory for original research; repository deposition strongly encouraged; data citation with persistent identifier recommended | Confirmed |
| Manuscript in .docx or LaTeX; supplementary spreadsheets as .csv or .xlsx; supplementary text as PDF; supplements referred to as "Online Resource" | Confirmed |
| LLM use must be documented in Methods; LLMs cannot be authors | Confirmed — newly relevant, not in the original brief |
| An allogeneic engineered Tr1 cell therapy is in an ongoing Phase 1/2a open-label dose-escalation study in moderate-to-severe **treatment-refractory** Crohn's disease (NCT06721962, "RESTORE"), recruiting, adults 18–65, ≥40 kg; administered as a single ~30-minute IV infusion with possible overnight observation; the registered study also involves cyclophosphamide | Confirmed. **No efficacy results published.** Company communications have anticipated initial safety and efficacy data in late 2026 — treat as a timeline signal only |
| The same product is in a parallel Phase 1 trial for GvHD prevention (NCT06462365); a related CAR-Tr1 product for progressive MS has received IND clearance | Confirmed — useful only as platform-plausibility context |
| Autologous polyclonal Treg therapy in type 1 diabetes (PTG-007): 54 patients from Phase I/II trials followed 7–12 years; a proportion insulin-independent for 18–24 months; a subset in clinical remission at 7–12 years; EMA/CHMP qualified the product for a centralised marketing-authorisation application in March 2026 | Confirmed **as company/press communications and conference presentation**. The company stated an intention to publish; no peer-reviewed long-term publication was located. Cite accordingly and do not present denominators that are not published |
| US Crohn's disease prevalence 305 per 100,000; ≈1.011 million Americans; IBD prevalence 721 per 100,000 ≈ 2.39 million; IBD incidence 10.9 per 100,000 person-years (15.9 including lower-probability algorithm) | Confirmed (Lewis et al., *Gastroenterology* 2023;165(5):1197–1205.e2) |
| Ryoncil (remestemcel-L), an FDA-approved allogeneic MSC cell therapy for pediatric SR-aGvHD, has a US list price of $194,000/infusion (~$1.55M for the full 8-infusion course) | Confirmed — company-set WAC, reported by BioWorld and Managed Healthcare Executive (`data/raw/market_comparator_cell_therapy_prices.csv`) |
| Tabelecleucel/Ebvallo, an allogeneic EBV-specific T-cell therapy (EU/UK-approved, not FDA-approved), was rated cost-effective by ICER if priced $143,900–$273,700 per treatment cycle | Confirmed — an HTA value-based price recommendation, not a company list price (the product isn't commercially priced in the US); reported by Managed Healthcare Executive (`data/raw/market_comparator_cell_therapy_prices.csv`) |
| CMS Comprehensive Observation Services (C-APC 8011), CY2026: national unadjusted Payment Rate $2,672.15 (RW 29.2310); CY2026 OPPS conversion factor $91.415 | Confirmed directly from primary sources — CMS's own "January 2026 Web Addendum A" and "2026 January Web Addendum B" (OPPS APCs / OPPS Payment by HCPCS Code for CY 2026), user-supplied 2026-08-05 and independently re-verified the same day by reading the primary CSV extractions of both zipped files directly (APC 8011 row read verbatim: SI J2, RW 29.2310, Payment Rate $2,672.15); the two files cross-validate to the cent (29.2310 × 91.415 = 2672.15) |
| Nine ustekinumab biosimilars launched in the US through 2025–2026 at reported 80–99% list-price discounts (Wezlana, Selarsdi, Pyzchiva, Yesintek, Otulfi, Steqeyma, Imuldosa, Starjemza, plus an unbranded ustekinumab-aekn at 95% off from 2026-01-01); ustekinumab biosimilar market share reached 42% in Q1 2026; Stelara carries a 2026 IRA Maximum Fair Price of $4,695 (~66% below list) | Confirmed by independent web search (Pharmacy Times, GaBI, BioSpace/MedImpact, Pearce IP on Samsung Bioepis' own Q1/Q3 2026 market reports, Healio) — `docs/decision_resolutions_2026-08-05.md` §12a. **Not used for the actual re-pricing** (§4.2): the $4,695 MFP figure is a Medicare Part D, 30-day-supply quantity on a different benefit pathway than this project prices UST maintenance from, and converting it would require guessing CMS's dosing-equivalence methodology; re-priced instead from real CMS Part B ASP biosimilar Q-code data (UST/IFX) and CMS NADAC biosimilar NDC data (ADA), both queried directly, 2026-08-05 |

### B.2 Taken from the project files, not independently verified

Aliyev et al. transition probabilities, adjustment ratios, PMPM cost derivations and utility mappings; the ten Ham et al. manufacturing cost figures and the workbook's adaptation of them ($19,917/dose retail from a 4× markup on a 50%-economies-of-scale COGS); CMS ASP and HCUP figures as entered in the workbook. All require source re-verification during Phase 1.

### B.3 Assumptions introduced by this design, flagged as such

The mixture-cure structure itself; the week-56 landmark for the cure split; SDR utility equal to Remission utility; a constant post-cure relapse hazard; cured patients never entering conventional therapy; non-cured Treg patients being efficacy-equivalent to ustekinumab; the 10-year decision horizon for population VOI; and the recommendation of a single-dose base case. None of these is a finding. Each is a modelling choice that should be named as such in the manuscript and tested in Section 10's scenarios.

### B.4 Deliberately not supplied

No numeric value is proposed anywhere in this document for the Treg cure fraction or the post-cure relapse hazard. The published analogs support a qualitative range and, at best, an unpublished denominator. Assigning a point estimate to them and presenting it as parameterisation "informed by" those programmes would misrepresent the evidence. Section 7.2 sets out the two defensible alternatives: documented expert elicitation, or an explicit uniform prior with the headroom analysis carrying the result.
