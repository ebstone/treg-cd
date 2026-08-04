# Data Dictionary

Source data for the Treg-vs-biologics Crohn's disease CEA (Jadambaa, Stone,
Abraham), pulled together from the manuscript, the Aliyev et al. 2019
(*Pharmacotherapy*) supplementary appendices, and the live `IBD_CEA_v6_PSA.xlsm`
model workbook.

## `data/raw/` — external, published source extracts

| File | Contents | Source | Caveats |
|---|---|---|---|
| `aliyev2019_appendixS1_table1_assumptions.csv` | 29 model assumptions + justifications (Induction, Maintenance, Quality of Life, Costs) | Aliyev et al. 2019, Appendix S1, Supplementary Table 1 | The source PDF's text layer extracts in scrambled (non-row) order. Rows marked `verified=yes` were confidently re-paired by content; 6 rows marked with a `LOW CONFIDENCE` note have justification text that is ambiguous or recurs verbatim across multiple assumptions — **check these against the original PDF table image before relying on them.** |
| `aliyev2019_appendixS1_table2_parameters.csv` | 59 source parameters (trial endpoints, adjustment ratios, costs, utilities) with Beta/Gamma/Uniform distributions | Aliyev et al. 2019, Appendix S1, Supplementary Table 2 | Extracted cleanly (row order intact); high confidence. Only the active-drug-arm trial endpoints — does not include each trial's separate placebo-arm endpoint, which Appendix S2's own adjustment-ratio method requires; not sufficient on its own to independently re-derive `aliyev2019_appendixS2_table3_...`/`..._table4_...` below (`R/00_derive_transition_probs.R` header comment). |
| `aliyev2019_appendixS2_table3_induction_transition_probabilities.csv` | Full induction-phase transition matrix (4x4, UST/IFX/ADA), long format | Aliyev et al. 2019, Appendix S2, Supplementary Table 3 | Transcribed 2026-08-04 and verified by direct visual comparison against the table image in `phar2208-sup-0002-appendixs2.docx`/`.pdf` (Wiley/*Pharmacotherapy* supplementary material — not committed to this repository; cite via the journal DOI rather than redistributing the file). ADA and IFX rows are numerically identical in the published table (base-case IFX:ADA efficacy ratio = 1.00), confirmed against the image, not a transcription error. Used directly by `R/00_derive_transition_probs.R`, no conversion. |
| `aliyev2019_appendixS2_table4_maintenance_transition_probabilities.csv` | Full maintenance-phase transition matrix (UST/IFX/ADA/CT), 2-week cycle, long format | Aliyev et al. 2019, Appendix S2, Supplementary Table 4 | Same transcription/verification note as Table 3 above. Biologic arms (UST/IFX/ADA) have no `Moderate-Severe` row by design (non-responders exit to the CT track rather than continuing on this matrix). `R/00_derive_transition_probs.R` converts this 2-week matrix to the study's 8-week cycle via exact Markov-chain matrix power (`M_2wk^4`), not DEALE. |
| `manuscript_supplement1_health_states_and_utilities.csv` | CDAI-defined health states and their EQ-5D utilities as used in *this* manuscript | Final manuscript, Supplement 1 | |
| `manuscript_supplement1_induction_transition_probabilities.csv` | UST/IFX induction transition probabilities (= Aliyev Supplementary Table 3) plus the authors' derived TREG induction row | Final manuscript, Supplement 1 | TREG row is the authors' own construction (20% reduction in adverse transitions vs. UST), not a published trial result. |
| `manuscript_supplement1_costs_and_discounting.csv` | Sample induction cost calc, per-cycle monitoring/management costs, discount rate/cycle parameters | Final manuscript, Supplement 1 | |
| `manuscript_table1_base_case_results.csv` | Base-case QALYs, costs, ICERs | Final manuscript, Table 1 | |
| `manuscript_table2_nmb_results.csv` | Net monetary benefit at WTP=$150,000, under Base/Alpha(10%)/Beta(30%) Treg-advantage scenarios | Final manuscript, Table 2 + Supplement 1 sensitivity tables | |
| `tenham2020_manufacturing_cost_casestudies.csv` | 8 published academic cell-therapy manufacturing cost case studies, 2018 EUR | ten Ham et al. 2020, *Cytotherapy* | Raw published figures; currency/inflation conversions are in `data/processed/model_tenham_derived_treg_dose_cost.csv`. |
| `cms_asp_and_hcup_cost_sources.csv` | CMS ASP drug unit costs, CMS Physician Fee Schedule administration cost, HCUP colectomy surgery cost, FSS unit costs | CMS, HCUPnet, Aliyev et al. 2019 Table 2 | URLs are general program pages, not deep links to the specific pricing file used (the workbook cites monthly pricing files by name but the underlying files were not included in the project). |

## `data/processed/` — extracted directly from `IBD_CEA_v6_PSA.xlsm` (`Input Parameters`, `Sensitivity Analysis`, `PSA`, `Treg Cost official` sheets)

| File | Contents | Notes |
|---|---|---|
| `model_maintenance_transition_probabilities.csv` | Full 8-week-cycle maintenance transition matrices for TREG (Base/Alpha=10%/Beta=30% scenarios), UST, IFX, and CT (179 rows, long format) | This is what the model actually runs on. Values for IFX/UST/CT are the DEALE-converted 8-week versions of the Aliyev 2-week probabilities; the conversion formulas themselves are not present in the workbook (hardcoded), so they could not be independently re-derived here — see project notes. |
| `model_health_state_costs.csv` | Per-8-week-cycle cost by health state, for each of TREG/UST/IFX/CT | |
| `model_dose_costs_and_psa_ranges.csv` | TREG and IFX per-dose drug costs plus their PSA (Gamma) ranges | |
| `model_health_utilities.csv` | Health-state utility values as calculated in the workbook | |
| `model_utility_ratios_and_psa.csv` | The 4 utility ratios (Mild:Remission, M-SR:Mild, M-S:M-SR, Surgery:M-SR) plus the Remission base utility, each with its PSA (Uniform) range | |
| `model_key_parameters.csv` | Study population (1,000) and annual discount rate (3%) | |
| `model_univariate_sensitivity_tornado.csv` | The 9 parameters varied in the manuscript's tornado plot (Figure 3), with NMB for IFX/TREG/UST at base/low/high for each | Reproduces the manuscript's stated finding that Treg efficacy advantage and Treg cost are the top two NMB drivers. |
| `model_psa_parameter_distributions.csv` | The 7 parameters actually sampled in the workbook's probabilistic sensitivity analysis (5 utility ratios, Uniform; 2 drug costs, Gamma) | **This PSA does not sample the Treg efficacy-advantage parameter or any transition probability**, despite the efficacy advantage being the top tornado driver — a known limitation flagged in prior model review, not fixed here. |
| `model_tenham_derived_treg_dose_cost.csv` | The workbook's derivation from ten Ham's Case Study 8 (materials/personnel/facility costs) through the modeled scale-up assumptions (44→22 doses/batch, 50% materials cost reduction, 4x retail markup) to the final $4,979 manufacturing cost / $19,917 retail price per dose | Reconciles exactly with the manuscript's stated figures. |

## Known open issues (not resolved by this data pull — see prior project notes)

- The `Economically Justifiable Price` sheet in the workbook currently computes **$185,625**, which does not match the manuscript's stated EJP of **$95,382**. Not corrected here; flagged for the team.
- The maintenance Markov sheets do not implement Aliyev's 2-year biologic-duration cap (patients stay on biologic-specific transition probabilities for the full ~6.15-year horizon rather than switching to CT-derived probabilities after cycle 13).
- The Treg base-case surgery-state transition row does not follow the same redistribution rule as the Alpha/Beta scenario rows (see `model_maintenance_transition_probabilities.csv`, TREG base scenario, `from_state=Surgery`).
