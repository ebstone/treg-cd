# data/

Source data for the Treg-vs-biologics Crohn's disease cost-effectiveness
analysis, pulled together from the manuscript, the Aliyev et al. 2019
supplementary appendices, and the `IBD_CEA_v6_PSA.xlsm` model workbook, per
the project's `data/raw/` + `data/processed/` convention (Section 12.2 of the
research project design).

- **`raw/`** — verbatim extracts from published sources (Aliyev 2019, ten Ham
  2020, CMS, HCUP) and from the manuscript's own tables/supplement. Nothing
  in this folder is computed; it is transcription only.
- **`processed/`** — parameter values as the live Excel model actually uses
  them, extracted directly from `IBD_CEA_v6_PSA.xlsm` via `openpyxl` (cached
  cell values, not formulas). Treat this as a snapshot of the current model,
  not as independently re-derived or verified numbers.
- **`data_dictionary.md`** — one row per file: what it contains, where it
  came from, and any caveats (a few are worth reading before you use the
  file — see especially the Aliyev Table 1 assumptions file and the PSA
  parameter list).

Both `.xlsm` workbooks themselves (the frozen reference implementation)
should also be committed unchanged, per the plan — they are the ground truth
these CSVs were extracted from.

Original workbook: `IBD_CEA_v6_PSA.xlsm` (also `IBD_CEA_v6_Univariate_Sensitivity_Analysis.xlsm`,
same `Input Parameters`/model structure, different PSA sample size).
