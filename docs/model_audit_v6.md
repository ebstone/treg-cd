# Audit of the IBD CEA v6 workbooks

Findings from a direct read of `IBD_CEA_v6_PSA.xlsm` and `IBD_CEA_v6_Univariate_Sensitivity_Analysis.xlsm` (formulas and cached values). Cell references are given so each can be checked independently. Severity reflects likely effect on reported results. Direction of bias is stated where it can be determined.

Extracted from `docs/analysis_plan.md` v1.0, Appendix A, so it stands as its own record of what changed and why as the rebuild in `R/` proceeds. Cross-referenced from the analysis plan's Section 0.2 and Section 15 decisions.

## A.1 Errors

**A1 — Costs are discounted against cumulative deaths, not elapsed time. [High severity]**
In all three Markov sheets, the discounted-cost column uses `=AQn*(1/(1+discount)^Nn)`. Column N holds *cumulative conventional-therapy deaths*; elapsed time in years is column A (`=Bn*(8/52)`). The QALY column is correct (`=ADn*(1/(1+discount)^An)`). At cycle 7 the exponent applied to costs is 0.343 (deaths) instead of 1.077 (years); at the final cycle the implied discount factor is ≈0.899 rather than the correct ≈0.833.
*Effect:* costs are systematically under-discounted in every arm. Because the comparators' costs are recurring and Treg's are front-loaded, this inflates comparator costs relative to Treg — **biases in favour of Treg**.

**A2 — Base-case Treg "Surgery" transition row appears column-shifted. [High severity]**
`Input Parameters` C34:C39 give Surgery → M-S 0.042, M-SR 0.0276, Mild 0.228, Remission 0, Surgery 0.7024. Ustekinumab's corresponding row is 0.0405 / 0.0267 / 0.0267 / 0.6776 / 0.2283, and the workbook's own Alpha (10%) and Beta (30%) Treg sheets give Remission 0.706 / 0.763 and Surgery-self 0.205 / 0.160. Under the documented 20% rule the base case should be approximately Remission 0.73, Surgery-self 0.18, Mild 0.03. The values for Mild, Remission and Surgery appear to have been misplaced.
*Effect:* Treg patients entering the Surgery state largely remain there, accruing the surgery state cost ($34,679/cycle) and the lowest utility. **Biases against Treg.** Correcting it will improve Treg's reported cost-effectiveness; the magnitude must be quantified.

**A3 — Second Treg dose is charged to the wrong cohort. [Moderate]**
`TREG Markov` AJ10 = `SUM(AF10:AI10)+IF(treg_total_doses>1, H16*cost_treg_dose, 0)`. The cost is charged in row 10 (cycle 7, ≈1.08 years — correct timing) but multiplied by `H16`, the number alive in the Treg track at row 16 (cycle 13, ≈2 years). It should be `H10`.
*Effect:* undercounts the patients receiving dose 2. **Biases in favour of Treg.** Moot if the single-dose base case is adopted.

**A4 — Dollar years are mixed across arms. [Moderate]**
Treg health-state costs are un-inflated 2017 values (`costs_treg_MSR` 217, `costs_treg_mild` 91, `costs_treg_remission` 10) while ustekinumab and infliximab use 2025 values ($282.86, $118.62, $13.04). See analysis plan §8.
*Effect:* Treg's non-drug maintenance costs understated by ≈23%. **Biases in favour of Treg** (small in absolute terms relative to drug costs).

**A5 — Infliximab M-SR state cost appears to double-count administration. [Low–moderate]**
`Input Parameters`: `costs_ifx_MS` = dose + 402, `costs_ifx_mild` = dose + 238, `costs_ifx_remission` = dose + 133 — each equal to the health-state cost plus one $119.36 administration fee. But `costs_ifx_MSR` = dose + 522 ≈ 402.22 + 119.36, i.e. the M-S figure *plus a second* administration fee.
*Effect:* overstates infliximab costs in the M-SR state. Biases against infliximab.

**A6 — Infliximab induction cost differs between the two workbooks. [High severity for the PSA]**
Univariate workbook: `cost_inx_MS` = $3,563/patient (hardcoded; consistent with three induction infusions at weeks 0, 2 and 6). PSA workbook: `=IF(psa_switch=0,1068,psa_costs_ifx_dose)+…` = $1,372/patient — a single infusion.
*Effect:* the PSA understates infliximab induction cost by ≈$2,191/patient. The probabilistic results are therefore **not** the stochastic counterpart of the deterministic results. Any PSA-based claim in the current manuscript is affected.

**A7 — The Economically Justifiable Price worksheet is disconnected from the model. [High severity]**
`Economically Justifiable Price` D9 = 26.48 (QALY Treg), D10 = 25.73 (QALY IFX), D15 = 73,125.2 (Costs IFX) — all hardcoded. The QALY values appear nowhere else in the workbook (the model's discounted per-patient QALYs are ≈4.07 and ≈3.95). The sheet returns EJP = $185,625, while the manuscript reports $95,382. The "graph" column additionally computes `(price − Costs_IFX)/ΔQALY`, treating the x-axis variable as *total cost*, not price.
*Effect:* the reported EJP cannot be reproduced from the workbook. This number must be recomputed from scratch.

**A8 — The conventional-therapy track omits health-state costs. [Moderate]**
`costs_ct_MS` = `costs_ct_MSR` = `costs_ct_mild` = `costs_ct_remission` = 88, i.e. a flat per-cycle cost in every non-surgical state. In Aliyev, the CT cycle cost ($67, 2017 USD) is the CT *drug* cost, applied **in addition to** health-state costs.
*Effect:* understates CT-track costs in all arms. Because the share of the cohort on the CT track differs sharply by arm (≈79% for UST, ≈76% for IFX, ≈34% for Treg at induction), the bias is not symmetric — it **favours the comparators**, partially offsetting A1 and A4.

**A9 — Surgery is charged as a full episode cost in every cycle of occupancy. [Moderate]**
`costs_*_sur` ≈ $35,518, an HCUP colectomy episode cost, is multiplied by the Surgery-state occupancy each cycle. With an 8-week probability of remaining in Surgery of 0.228 (and 0.702 in the Treg base case, see A2), a meaningful share of the cohort is charged a full surgical episode repeatedly. Aliyev used a per-cycle state cost derived from PMPM data ($884, 2017 USD).
*Effect:* likely overstates surgical costs in all arms; interacts badly with A2. Needs an explicit decision: tunnel state with a one-time episode cost, or a per-cycle state cost.

*Resolved 2026-08-04:* **per-cycle state cost, not a one-time episode cost** — matching how every other health state in Aliyev's design is costed, and matching Aliyev's Costs Assumption #1 (Suppl. Table 1) that health-state costs are constant across all biologics and CT. Aliyev's Surgery state is not a colectomy proxy: it is explicitly documented as covering "all surgeries and procedures" (Lichtenstein et al. 2005), most performed outpatient, costed via the Malone et al. severe-fulminant PMPM claims category — chosen only because it was the sole category with reported inpatient costs, not because Surgery-state occupants are assumed to be undergoing major inpatient operations. The native 2-week Aliyev matrix (Suppl. Table 4) confirms the state is short-staying by design (Surgery→Remission = 0.868, Surgery→Surgery = 0.098 per 2-week cycle; expected occupancy ≈1.1 cycles), consistent with a recurring claims-average cost rather than a one-time DRG charge.

Base case now uses Aliyev's own per-cycle Surgery figure ($884, 2017 USD) inflated by the same ≈1.3035 factor already applied to Aliyev's other per-cycle costs (Decision 5, analysis_plan.md §8) → **$1,152.29 (2025 USD)**, applied identically across TREG/UST/IFX/CT in `data/processed/model_health_state_costs.csv`. The HCUP colectomy figure is retained in `data/raw/cms_asp_and_hcup_cost_sources.csv` for provenance and as a candidate input for a future episode-cost structural scenario, but is no longer used in the base case.

**Open item carried forward:** Aliyev's own materials do not cleanly reconcile on this figure — the manuscript appendix worked example states $884 (2017 USD) while Suppl. Table 2 lists the underlying Severe-Fulminant PMPM at $1,475 (nominally "2008 USD," though the table's footnote states all costs were converted to 2017 USD). $884 was used here for consistency with how the model already treats Aliyev's other four per-cycle costs (M-S/Mild/Remission/CT), which come from the same appendix worked-example section rather than Table 2's raw PMPM inputs — but the derivation from $1,475 (or from the underlying $374/$123/$111 building blocks) to $884 is not independently reconstructable from the materials on hand. Verify against the original Aliyev appendix table image, or with the Aliyev co-authors if reachable, before this goes to submission. This is the same unresolved reconciliation already noted for M-S/Mild/Remission/CT (data_dictionary.md, analysis_plan.md §8) — Surgery does not add a new gap, it inherits the existing one.

## A.2 Reproducibility and documentation gaps

**A10 — Treg transition probabilities are hardcoded and the documented method does not reproduce them.**
All Treg values in `Input Parameters` C10:C39 (and the Alpha/Beta blocks) are literal numbers, not formulas referencing ustekinumab. The **maintenance** rows are approximately consistent with the documented 20% rule. The **induction** row is not: base-case M-S → M-S is 0.3368 against ustekinumab's 0.791, a ≈57% relative reduction rather than 20%. The Alpha/Beta/base induction values are linear in the stated advantage parameter but imply a zero-advantage intercept of ≈0.421, not 0.791. Since the induction split determines what fraction of the cohort ever reaches the biologic track at all, this is the single most influential set of numbers in the model, and it is currently unreproducible from the manuscript's stated method.
*Action:* re-derive in code, or document the actual derivation. Do not carry these values forward unexamined.

**A11 — Cycle-length conversion cannot be verified.** All 8-week transition probabilities are hardcoded; the DEALE intermediate calculations are not in the file. See analysis plan §7.3.

**A12 — No half-cycle correction.** Cohort state occupancy is applied at cycle boundaries with no correction.

**A13 — Cycle count and cycles-per-year are misstated in the manuscript.** The workbook runs 40 cycles (rows 4–43); the manuscript states 39. Forty 8-week cycles = 6.15 years, so the workbook is internally consistent and the manuscript text is wrong. Separately, Supplement 1 states "# cycles per year 6.15"; the correct figure, used by the workbook, is 6.5 (52/8). The 6.15 is the horizon in years, not a rate.

**A14 — The PSA samples only seven parameters.** `PSA` rows 3–12: five utility parameters (Uniform over their CI bounds) and two drug unit costs (Gamma). No transition probability, no efficacy parameter, no health-state cost and no utility-for-Treg-specific-state is sampled. ICERs are computed per iteration (column R) rather than incremental net monetary benefit.
*Effect:* the reported PSA characterises uncertainty in utilities and drug prices only. It cannot support EVPI or EVPPI, which are undefined over parameters held fixed. Rebuilding the PSA is a prerequisite for the primary endpoints of the redesign, not an optional improvement.

## A.3 What the audit implies for the reported results

The manuscript's base-case numbers (IFX $73,124 / 3.951 QALYs; Treg $80,855 / 4.068; UST $118,672 / 3.904; ICER $65,918; NMB $519,527 / $529,389 / $466,926 per 1,000 patients) are exactly reproducible from the univariate workbook's cached values, so the reported figures are internally faithful to that file. The issue is not arithmetic slippage between model and manuscript; it is that the file itself contains the defects above.

Note also that the NMB column in `Results` carries both a live formula and a hardcoded "Static" column holding precisely the manuscript's reported values, while the live formula in the PSA workbook returns different numbers. Whatever the history, the practical consequence is that the manuscript's results are a frozen snapshot that the live workbook no longer reproduces. This, more than any individual defect, is the argument for rebuilding in scripted, version-controlled code.
