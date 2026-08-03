# CHEERS 2022 Checklist

Working checklist for submission as an Online Resource, per *PharmacoEconomics*' requirement that economic evaluations follow CHEERS 2022 (Husereau et al., *PharmacoEconomics* 2022, doi:10.1007/s40273-021-01112-8).

This is a draft mapping carried over from `docs/analysis_plan.md` Section 11. Verify item numbering and wording against the published checklist before submission, and fill in the "Manuscript location" column once the manuscript is drafted — the "Where addressed" column below currently points at analysis-plan section numbers, not manuscript sections.

| # | CHEERS 2022 item | Where addressed (analysis plan) | Manuscript location | Status |
|---|---|---|---|---|
| 1 | Title | §1.1 | | State that it is an economic evaluation and name the interventions |
| 2 | Abstract | §2 | | Structured; may run to 450 words |
| 3 | Background and objectives | Manuscript intro; §3 | | |
| 4 | Health economic analysis plan | This document (`docs/analysis_plan.md`) | | Cite it: state that an analysis plan was developed a priori and deposited in the repository |
| 5 | Study population | §5 | | Must state both populations and which is base case |
| 6 | Setting and location | §4.1 | | US |
| 7 | Comparators | §4.1, Decision 2 | | Justify inclusion/exclusion of ADA explicitly |
| 8 | Perspective | §4.1 | | Healthcare sector base case; societal scenario |
| 9 | Time horizon | §4.1 | | Justify lifetime and report the 6.15-year comparability scenario |
| 10 | Discount rate | §4.1 | | 3%, with scenarios |
| 11 | Selection of outcomes | §3 | | QALYs; plus EJP/EVPI/EVPPI as decision outcomes |
| 12 | Measurement of outcomes | §7.1 | | CDAI-defined states |
| 13 | Valuation of outcomes | §7.1 | | Buxton 2007 mapping; note it is a mapping, not directly elicited utilities |
| 14 | Measurement/valuation of resources and costs | §7.1, §8 | | |
| 15 | Currency, price date, conversion | §8 | | **Currently non-compliant** — index unnamed, dollar years mixed |
| 16 | Rationale and description of model | §6 | | Diagram required; cure structure must be justified against precedent in oncology cell-therapy HTA |
| 17 | Analytics and assumptions | §6, §7 | | Half-cycle correction; DEALE re-derivation; Dirichlet sampling |
| 18 | Characterising heterogeneity | §5 | | Population scenarios serve this |
| 19 | Characterising distributional effects | — | | **Gap.** Not addressed and probably out of scope; state so explicitly rather than omitting |
| 20 | Characterising uncertainty | §9, §10 | | **New obligation:** report structural uncertainty around the cure assumption as a category distinct from parameter uncertainty; π = 0 is the honest structural bound |
| 21 | Engagement with patients and others | §7.2 | | If expert elicitation is run, describe the protocol here; if patients are not engaged, say so |
| 22 | Study parameters | §7.1 + ESM | | Full parameter table with distributions and sources, as an Online Resource |
| 23 | Summary of main results | Results | | Disaggregated costs and QALYs by arm |
| 24 | Effect of uncertainty | §9, §10 | | CEACs, EVPI/EVPPI, EJP credible interval |
| 25 | Effect of engagement | §7.2 | | How elicitation altered the analysis, if applicable |
| 26 | Findings, limitations, generalisability | §13 | | |
| 27 | Source of funding | Declarations | | |
| 28 | Conflicts of interest | Declarations | | Declare any relationship with cell-therapy developers, including none |

## Additional journal requirements

- Use of a large language model must be documented in the Methods section. LLMs cannot be listed as authors. Agree the wording early.
- Declarations must appear as a "Declarations" section before the reference list, covering Funding, Competing interests, Ethics approval, Consent, Data and code availability, and Author contributions. This is a modelling study using published aggregate data and requires no IRB review; state the exemption and its basis explicitly.
