# Cost and utility attachment by health state and arm (analysis_plan.md §7.1, §8). Combines a
# per-cycle occupancy trace (R/02_markov_engine.R / R/03_cure_fraction_module.R's on_biologic/
# on_ct/on_sdr matrices) with per-cycle costs and utilities to produce discounted QALYs and
# non-drug costs. 2025 USD throughout; Aliyev's native 2-week cycle throughout (no cycle-length
# conversion anywhere in this pipeline -- R/00_derive_transition_probs.R's third-revision header
# comment).
#
# FIRST PASS, 2026-08-04. Implements utility attachment and the arm-independent per-cycle
# health-state monitoring cost in full. UST/IFX drug acquisition + administration costs added
# 2026-08-04 (second pass). ADA dosing sourced 2026-08-04 (third pass) and its drug-cost layer
# wired in the same day (fourth pass) -- see "UST/IFX/ADA drug acquisition/administration costs"
# below. Treg's own dose cost sourced and wired in the same day (fifth pass) -- see "Treg's own
# dose cost" below; NOT fully resolved, see that section for what's still open (cyclophosphamide
# dose; observation-stay cost PRICE now fully confirmed at $2,672.15, CY2026, directly from CMS's
# own Addendum A -- wired in as an explicit opt-in `treg_dose_cost()` argument defaulting to 0,
# since whether a Treg infusion actually qualifies for this billing category is still an open
# clinical/modelling question, not a pricing one). Non-drug-cost outputs from the first pass
# ("non_drug_cost_by_cycle" etc.) are kept as separately named fields even now that drug costs
# exist elsewhere in this module -- don't conflate the two.
#
# ---- Why this doesn't read data/processed/model_health_state_costs.csv ------------------------
#
# That file's non-Surgery rows (Moderate-Severe/M-SR/Mild/Remission) are extracted as-is from the
# v6 Excel workbook (data_dictionary.md: "extracted directly from IBD_CEA_v6_PSA.xlsm") and are
# NOT simply "Aliyev's per-cycle health-state cost" -- they are workbook-era composites that bake
# in an implicit drug-acquisition assumption for the biologic arms. E.g. UST's $14,312 M-SR
# figure decomposes as $282.86 (Aliyev's own state cost, see below) + $14,029.14 -- and
# $14,029.14 = $155.883/mg (data/raw/cms_asp_and_hcup_cost_sources.csv) x 90mg, i.e. a full q8w
# UST maintenance dose, matching the "$14,029 per cycle for ustekinumab" analysis_plan.md §6.4
# cites. That composite was computed under the workbook's original 8-week-cycle design, where one
# cycle *was* one dosing interval. This project now runs a native 2-week cycle: naively applying
# that composite every 2-week cycle would charge a full q8w dose four times as often as it's
# actually given -- a ~4x overstatement of UST/IFX drug cost, not a rounding issue. Same class of
# bug as A9 (Surgery) and the stale cost_per_8wk_cycle_usd_2025 label (both resolved 2026-08-04,
# data/data_dictionary.md), just not caught until now because nothing had consumed this file yet.
#
# The CT track has the mirror-image problem: A8 (docs/model_audit_v6.md) already documents that
# the workbook's flat $88 CT-track figure is the CT *drug* cost only ($67 2017 USD x ~1.306 =
# $87.53), not a health-state cost -- the workbook omits the state-specific monitoring cost A8
# says should be added on top. Reusing $88 as if it were a state cost would carry A8's bug
# forward, not fix it.
#
# So this module treats "cost" as two independently-sourced layers, matching how Aliyev actually
# costed his own model (Appendix S2; analysis_plan.md §8's reconciliation table):
#   1. A monitoring/management cost that varies by health state but is constant across every arm
#      (Aliyev Costs Assumption #1 -- already established for Surgery in the A9 fix). Sourced
#      directly from Aliyev's own reported 2017 USD per-cycle figures below, not from the
#      workbook's arm-specific composites.
#   2. An arm-specific drug-acquisition/administration cost, layered on top only for arms and
#      cycles where a dose is actually given. Implemented for UST/IFX/ADA -- see below -- and,
#      partially, Treg (see "Treg's own dose cost" below for what's still open).
#
# ---- UST/IFX/ADA drug acquisition/administration costs (UST/IFX added 2026-08-04; ADA added ---
# ---- the same day, in a later pass) -------------------------------------------------------------
#
# Sourced from data/raw/biologic_dosing_schedule.csv (dose size, route, and frequency, cited to
# the STELARA/REMICADE/HUMIRA FDA-approved Prescribing Information) and
# data/raw/cms_asp_and_hcup_cost_sources.csv ($/mg prices + the flat IV administration fee).
# UST/IFX prices are CMS ASP-methodology Medicare Part B payment limits; ADA's is NOT -- see
# point 3 below. Three things this pairing needed that neither file resolves on its own:
#
# 1. **Patient weight**, for UST's induction tier lookup and IFX's mg/kg dosing (both weight-
#    based; ADA's dosing is flat/fixed, no weight dependency at all). Uses
#    ASSUMED_PATIENT_WEIGHT_KG below, Aliyev's own cohort mean (analysis_plan.md §5: "mean age
#    35, 71 kg, 50% male") -- not manuscript_supplement1_costs_and_discounting.csv's 70kg worked
#    example, which is just a rounded illustration of that same figure, not an independent
#    parameter.
# 2. **Vial rounding for IFX**: 71kg x 5mg/kg = 355mg doesn't divide evenly into IFX's 100mg
#    vials (data/raw/aliyev2019_appendixS1_table2_parameters.csv). Rounded UP to the nearest whole
#    vial (400mg), per Aliyev's own Induction Assumption #8 ("No vial sharing... Excess
#    medication is assumed to be wasted," data/raw/aliyev2019_appendixS1_table1_assumptions.csv)
#    -- applied to every IFX dosing event (induction AND maintenance), not just induction, since
#    Remicade vials are single-use/non-preserved and every administration independently wastes
#    any partial vial. UST's induction doses (260/390/520mg) and ADA's doses (160/80/40mg) are
#    already exact vial/pen multiples (130mg and 40mg respectively) by design -- no rounding
#    needed for either.
# 3. **A current ADA price that isn't a CMS ASP figure.** ADA is entirely self-administered SC --
#    no dose (including induction) has an applicable Medicare Part-B ASP price, which is why
#    Aliyev et al. 2019 itself prices ADA via an FSS unit cost rather than ASP. Sourced instead
#    from CMS Medicaid's NADAC (National Average Drug Acquisition Cost) -- $84.1678/mg, from the
#    40mg pen (matching the same 40mg unit Aliyev prices ADA against), effective 2026-07-22,
#    queried directly from the CMS/Medicaid open-data API, cross-validated against the 80mg pen's
#    independently-computed $/mg rate (within 0.03%) -- data/raw/cms_asp_and_hcup_cost_sources.csv.
#
# Cost is split the same way the module header already explains: a one-time INDUCTION cost
# (induction_drug_cost(), decision-tree-level, not part of any per-cycle trace -- analysis_plan.md
# §6.1's 8-week induction window) and a recurring MAINTENANCE cost
# (maintenance_drug_cost_by_cycle(), charged only on cycles a real dose is actually given, and
# only against mass still on that arm's own biologic track, so patients who've switched to CT or
# hit the 2-year cap stop being charged automatically, with no separate cap-awareness needed
# here). Dosing cadence read from `schedule`, not hardcoded: UST/IFX maintenance is every 4th
# native 2-week cycle (q8w); **ADA maintenance is every SINGLE cycle (EOW, interval=1) -- a
# materially different cadence, not a typo.** ADA's induction is two unequal sequential doses
# (160mg then 80mg), summed from two separate schedule rows rather than one row x n_doses like
# IFX's three identical doses.
#
# ---- Biosimilar re-pricing (2026-08-05) -- base case moved off originator pricing --------------
#
# The three points above describe the ORIGINAL 2026-08-04 sourcing pass; as of 2026-08-05
# (docs/model_audit_v6.md; docs/analysis_plan.md §15/§4.2; S8, §10.3) that pricing is no longer the
# base case. At pi=0 the Treg arm collapses to a one-time-cost UST-equivalent product, so the EJP
# is essentially the discounted comparator drug spend displaced -- comparator acquisition price is
# therefore the single largest lever on every headline number in this project, and originator
# (pre-biosimilar) ASP/NADAC pricing was the wrong base case for a therapy that would reach an
# adoption decision well after ustekinumab, infliximab and adalimumab biosimilars matured (nine
# ustekinumab biosimilars launched in the US through 2025-2026 at reported 80-99% list-price
# discounts; infliximab's biosimilar market has been established for years).
#
# `load_drug_prices()` below now takes a `pricing_basis` argument (default `"biosimilar_2026"`),
# and `data/raw/cms_asp_and_hcup_cost_sources.csv` carries two rows per drug accordingly:
#   - **UST induction**: median CMS Part B ASP across 6 ustekinumab biosimilar Q-codes (a real
#     current figure, not an estimate) -- $4.5375/mg, vs. $12.808/mg originator.
#   - **UST maintenance**: J3357 (the SC/maintenance HCPCS code the original sourcing pass used) no
#     longer carries a Part B ASP payment limit AT ALL in the current fee schedule -- independent,
#     unplanned confirmation that self-administered SC dosing isn't really a Part B ASP product,
#     the same reasoning point 3 above already applied to ADA. Unlike ADA, NADAC does not yet list
#     any ustekinumab biosimilar NDC either, so this row is a DERIVED PROXY (the UST induction
#     biosimilar:originator discount ratio applied to the previous SC originator ASP figure), not a
#     directly observed price -- see the CSV row's own note for the full chain and its limitation.
#   - **IFX**: median across 3 infliximab biosimilar Q-codes -- $2.6803/mg, a modest ~14% cut from
#     originator ($3.1041/mg the same quarter), consistent with a market that converged years ago
#     rather than one still in an early price war (contrast UST's much larger cut above).
#   - **ADA**: median CMS NADAC across 4 adalimumab biosimilar 40mg-pen NDCs -- $14.2482/mg, an
#     ~83% cut from the brand-only NADAC figure ($84.1678/mg) point 3 above sourced; same NADAC
#     methodology as before, only the product mix (brand-only -> biosimilar-inclusive) changes.
# All four biosimilar-era figures and their full derivations are in the CSV's own rows/notes, not
# repeated here. The pre-fix figures are retained verbatim under `pricing_basis =
# "originator_pre_biosimilar"` -- S8 (analysis_plan.md §10.3) inverts to that as the
# comparability-with-Aliyev scenario; pass it to `load_drug_prices()` and thread the resulting
# `prices` list into `run_comparator_arm_lifetime()`/`run_treg_arm_lifetime()`'s own `prices`
# override argument (R/05_deterministic_results.R) to run it -- no separate scenario-running
# harness exists yet for S1-S12 collectively (R/05's own header), and this fix doesn't add one.
#
# **Not applied to Treg's non-cured track.** attach_treg_costs_utilities() reuses UST's own
# maintenance matrix for Treg's non-cured patients (efficacy-equivalent to UST, analysis_plan.md
# §6.2) but that is an efficacy-only equivalence, not a claim that non-cured Treg patients are
# actually receiving ustekinumab. Confirmed with E. Stone, 2026-08-04: non-cured Treg patients
# are NOT charged UST's drug cost -- they received the Treg product (its own one-time/two-dose
# cost, still not wired in, see below), and giving them UST's ongoing drug cost on top would
# double up an acquisition cost they never actually incurred. This is why
# attach_maintenance_costs_utilities() takes an explicit `therapy` argument that must be supplied
# for the drug-cost layer to activate at all -- attach_treg_costs_utilities() deliberately never
# passes one for its Markov portion.
#
# ---- Treg's own dose cost (added 2026-08-04) ---------------------------------------------------
#
# treg_dose_cost() below is a one-time, decision-tree-level cost (R/00's third-revision note),
# same architectural pattern as induction_drug_cost() above -- not a per-cycle trace function,
# but still housed in this module rather than R/01_decision_tree.R, for consistency with that
# precedent. Three components, sourced to very different standards:
#
# 1. **Acquisition ($19,916.75/dose): fully traceable.** analysis_plan.md §7.1 item 12 already
#    carries "Derived" status (not "Action required") -- the ten Ham et al. 2020 manufacturing
#    costing framework's arithmetic is reproduced end to end in
#    data/processed/model_tenham_derived_treg_dose_cost.csv. treg_dose_cost() reads that file's
#    own final "Retail price per treatment/dose" line directly, NOT
#    data/processed/model_dose_costs_and_psa_ranges.csv's "TREG dose cost" `value_usd` column
#    ($18,908.74) -- found while wiring this in: those two figures disagree by a consistent
#    ~0.9494x factor, and the SAME ~0.9494x gap appears independently for that file's IFX dose
#    cost row too ($1,014.20 vs its own $1,068 psa_base), which rules out coincidence and points
#    to an unflagged inconsistency in the source workbook (comparable to A6/A7,
#    docs/model_audit_v6.md) -- not resolved here, just not silently picked around.
# 2. **Infusion administration: reuses the existing $57.90 CMS PFS rate** (prices$
#    iv_administration_usd, HCPCS 96365), same as UST induction/IFX/ADA's IV doses. Item 13
#    (§7.1) flags that "the reference programme administers a ~30-minute infusion" versus the
#    old model's 1-hour assumption -- but CPT/HCPCS 96365 covers "initial, up to 1 hour"
#    regardless of the exact duration within that window, so this doesn't actually change the
#    fee; noted here so it isn't mistaken for an unresolved gap.
# 3. **Cyclophosphamide preconditioning: PRICE sourced, DOSE is not.** Item 13 also flags that
#    "the trial protocol includes a cyclophosphamide preconditioning component," currently
#    modelled as zero. $/mg is now sourced (HCPCS J9075, data/raw/cms_asp_and_hcup_cost_sources.csv,
#    same WV Bureau for Medical Services ASP-methodology fee schedule already used for IFX). The
#    actual DOSE is NOT sourced: the cited reference trial (NCT06721962, "RESTORE") confirms "low
#    dose cyclophosphamide conditioning" in its public registry but does not publish a specific
#    mg or mg/kg figure. `cyclophosphamide_dose_mg` defaults to 0 in treg_dose_cost() -- a
#    known-incomplete placeholder, not a clinical claim that preconditioning is unnecessary or
#    free. Do not fill this in with a dose borrowed from a different program's regimen (e.g. a
#    CAR-T lymphodepletion protocol) without flagging it as a proxy, not a citation -- ask
#    E. Stone/the clinical co-authors first.
#
# **Overnight observation-stay cost: PRICE now fully sourced (2026-08-05, against CMS's own**
# **Addendum A), still NOT wired into the base case -- see why below.** CMS's Comprehensive
# Observation Services payment (C-APC 8011 under OPPS) is a bundled/packaged rate with specific
# qualifying criteria (>=8 hours of G0378, no status-T procedure same/next day), not a simple
# per-night figure. Historical national, wage-index-*unadjusted* base rates: $2,174.14 (CY2016),
# $2,349.66 (CY2018), $2,283.16 (CY2021), $2,439.02 (CY2023) -- all from independent secondary
# sources citing the CMS OPPS final rule. **The CY2026 figure is now a primary-source-confirmed
# fact, not an estimate**: CMS's own "January 2026 Web Addendum A" (OPPS APCs for CY 2026),
# user-supplied 2026-08-05, gives APC 8011 directly -- Relative Weight 29.2310, national
# unadjusted Payment Rate **$2,672.15** (SI J2). This exactly reproduces (to the cent) against
# the $91.415 conversion factor independently back-calculated from the companion Addendum B file
# supplied the same day (29.2310 x $91.415 = $2,672.15), which is about as strong a
# cross-validation as two CMS files can give each other. **$2,672.15 (CY2026, national,
# wage-index unadjusted) is therefore the sourced figure**, not a range, superseding the
# $2,439.02-ish planning range this header carried a few hours earlier. Two caveats remain, and
# are the reason this still isn't a base-case default: (a) $2,672.15 is the *national unadjusted*
# rate -- a specific facility's actual payment is wage-index-adjusted, higher or lower depending
# on region, a decision this project hasn't made; (b) C-APC 8011 pays per *qualifying stay*, not
# per Treg dose -- whether a Treg infusion actually triggers 8+ qualifying observation hours (vs.
# same-day discharge) is a clinical assumption this project hasn't made either, and is not a
# pricing question at all. Given that, treg_dose_cost() below still takes
# `observation_stay_cost_usd` as an explicit opt-in argument defaulting to 0 (same
# "sourced-price, opt-in-only" pattern already used for cyclophosphamide) -- the base case is
# unchanged, but a scenario run can now supply $2,672.15 directly, cited with confidence, without
# editing this function.
#
# **Base case is a single dose** (n_doses = 1), analysis_plan.md §4.1's recorded recommendation
# ("Single infusion base case; second dose as scenario") and §9.1's EJP formula (D-bar = 1 base
# case). The 2-dose STRUCTURAL SCENARIO (S6, §10.3) is deliberately not implemented here: the
# plan's own D-bar formula treats the second dose as conditional on the patient still being on
# therapy at that timepoint, and discounted accordingly -- that needs the per-cycle trace, not a
# one-time decision-tree cost, so belongs with the scenario-running logic (R/05/§10.3), not this
# function. treg_dose_cost() enforces n_doses == 1 rather than silently accepting a larger value
# and getting the 2-dose scenario's actual mechanics wrong.
#
# **Update, 2026-08-05:** S6 was subsequently formally DROPPED (R/09_scenarios.R's own module
# header, docs/analysis_plan.md §10.3/§10.4), not merely deferred as "not implemented here" above
# implied at the time it was written -- no dosing-interval assumption for a second dose exists
# anywhere in this project's sourcing, and single-dose was already the deliberately-chosen base
# case regardless. `n_doses == 1` is therefore this project's ONLY Treg dosing assumption now, not
# a base case awaiting a pending scenario alternative -- the enforcement above is permanent, not
# a placeholder.

if (!exists("MAINTENANCE_STATES")) source("R/utils/transition_matrix.R")
if (!exists("discount_factor")) source("R/02_markov_engine.R")
if (!exists("load_population_utility_norms")) source("R/utils/population_utility.R")

# ---- Health-state monitoring/management cost (Aliyev Appendix S2; analysis_plan.md §8) --------

#' Aliyev's own per-cycle health-state costs, 2017 USD, at his native 2-week cycle -- constant
#' across every arm (Costs Assumption #1, Aliyev Suppl. Table 1; already established for Surgery
#' in the A9 fix, data/data_dictionary.md, resolved 2026-08-04). M-SR uses the same figure as
#' Moderate-Severe ("assumed equal to M-S", analysis_plan.md §8's reconciliation table) -- Aliyev
#' does not report M-SR separately. Death is always 0.
ALIYEV_MONITORING_COST_2017_USD <- c(
  "Moderate-Severe" = 217, "Moderate-Severe Responder" = 217, "Mild" = 91,
  "Remission" = 10, "Surgery" = 884, "Death" = 0
)

#' 2017 -> 2025 USD inflation factor. Reproduces the workbook's published M-S/Mild/Remission
#' figures exactly (analysis_plan.md §8) but the underlying index is still uncited -- Decision
#' 5's first required fix (analysis_plan.md §15), not resolved by this module.
INFLATION_2017_TO_2025 <- 1.3035

#' Named vector of per-cycle monitoring cost (2025 USD) by health state, applied identically to
#' every arm. Using this uniformly for TREG as well as UST/IFX/CT also resolves A4/Decision 5's
#' second required fix (Treg's health-state costs were un-inflated 2017 dollars) as a side
#' effect of sourcing from Aliyev directly rather than from the arm-specific workbook rows.
health_state_monitoring_costs <- function() {
  ALIYEV_MONITORING_COST_2017_USD * INFLATION_2017_TO_2025
}

# ---- Societal-perspective productivity cost (scenario S9, analysis_plan.md §4.1/§10.3) ---------
#
# Manceur et al. 2020 (J Med Econ, data/raw/manceur2020_productivity_cost.csv): $2,168 incremental
# annual work-loss cost per Crohn's disease patient vs. matched non-IBD controls (claims-cohort
# analysis, n=6,715 CD / 33,575 controls, ~5-year average follow-up). Applied as a FLAT per-cycle
# addition to every living (non-Death) health state's monitoring cost, identically across every
# arm -- the source figure is a population average over a real CD cohort's mixed active-disease/
# remission experience, not stratified by momentary disease activity, so stratifying it by health
# state here would fabricate a precision the source doesn't support (this project's whole
# methodology exists to avoid exactly that). Whether Treg's own Sustained Deep Remission (cured)
# state should be exempted from this cost -- paralleling the SDR-utility argument, S7, that a
# durable cure removes the residual burden of managed disease -- is a stated, deliberately
# undecided modelling choice: societal_monitoring_costs() below applies the flat addition to SDR
# exactly like every other living state, the more conservative (larger societal-cost) choice,
# and does not zero it out.
ANNUAL_INDIRECT_COST_USD <- 2168

#' Per-cycle productivity/indirect cost (2025-vintage, no separate inflation applied: Manceur et
#' al.'s own study period runs through March 2017, comparable in currency-year terms to this
#' project's other 2017-sourced Aliyev figures before the same INFLATION_2017_TO_2025 factor would
#' apply -- but that factor is NOT applied here, a stated simplification: Manceur's own cost year
#' is not cleanly identified in the source the way Aliyev's is, so inflating it would imply a
#' precision this project doesn't actually have. Flagged, not silently assumed away.
productivity_cost_per_cycle <- function(cycle_weeks = 2) {
  ANNUAL_INDIRECT_COST_USD * cycle_weeks / 52
}

#' The societal-perspective variant of health_state_monitoring_costs() -- adds
#' productivity_cost_per_cycle() to every living state, Death excluded (Death has no ongoing
#' work-loss cost by construction). Scenario S9 (perspective: healthcare-sector vs. societal)
#' is exactly the choice between passing this function's output or health_state_monitoring_costs()'s
#' own output as the `monitoring_costs` argument everywhere that argument is already threaded
#' (trace_costs(), attach_maintenance_costs_utilities(), attach_sdr_costs_utilities()) -- no
#' changes needed to any of those functions themselves.
societal_monitoring_costs <- function(cycle_weeks = 2) {
  base <- health_state_monitoring_costs()
  addition <- ifelse(names(base) == "Death", 0, productivity_cost_per_cycle(cycle_weeks))
  base + addition
}

# ---- CT-track drug cost (Aliyev Appendix S2 item 18; A8 fix) ----------------------------------

#' Flat per-cycle CT (conventional therapy) drug cost, added on top of the state-specific
#' monitoring cost for every living CT-track state (A8, docs/model_audit_v6.md: the workbook's
#' $88 flat figure IS this drug cost, mislabeled as the whole health-state cost, with the
#' monitoring component missing entirely). $67 (2017 USD) inflated by the same factor as every
#' other line in this module, for internal consistency -- analysis_plan.md §8 reports this one
#' line's factor as "~1.306" (giving $87.53) rather than 1.3035 ($87.33); both are the workbook's
#' own published rounding of the same underlying number, not a discrepancy this module needs to
#' resolve.
CT_DRUG_COST_2017_USD <- 67
ct_drug_cost <- function() CT_DRUG_COST_2017_USD * INFLATION_2017_TO_2025

# ---- UST/IFX/ADA drug acquisition + administration cost (data/raw/biologic_dosing_schedule.csv) --

#' Assumed patient body weight for weight-based dosing (UST induction tier selection, IFX mg/kg
#' conversion). Aliyev's own cohort mean (analysis_plan.md §5: "mean age 35, 71 kg, 50% male"),
#' not the 70kg figure in manuscript_supplement1_costs_and_discounting.csv's worked example,
#' which is a rounded illustration of the same underlying population, not an independent source.
ASSUMED_PATIENT_WEIGHT_KG <- 71

#' Assumed patient starting age, same source and cohort as ASSUMED_PATIENT_WEIGHT_KG above
#' (analysis_plan.md §5: "mean age 35"). Used by the lifetime-horizon path
#' (R/05_deterministic_results.R's HORIZON_CYCLES_LIFETIME, R/utils/life_table.R's
#' death_prob_schedule()) as the cohort's cycle-0 attained age -- the point background mortality
#' starts accruing from.
ASSUMED_PATIENT_AGE_YEARS <- 35

#' IFX vial size (data/raw/aliyev2019_appendixS1_table2_parameters.csv: "FSS IFX Unit Cost (100
#' mg Vial)"). Used for vial-rounding -- see module header.
IFX_VIAL_SIZE_MG <- 100

#' UST/IFX/ADA dosing schedule (data/raw/biologic_dosing_schedule.csv, sourced 2026-08-04:
#' STELARA/REMICADE/HUMIRA Prescribing Information).
load_dosing_schedule <- function(raw_dir = "data/raw") {
  utils::read.csv(file.path(raw_dir, "biologic_dosing_schedule.csv"), stringsAsFactors = FALSE)
}

#' $/mg (or $/mg-equivalent) drug prices and the flat IV administration fee
#' (data/raw/cms_asp_and_hcup_cost_sources.csv).
#'
#' **Biosimilar re-pricing (2026-08-05, docs/model_audit_v6.md; docs/analysis_plan.md §15/§4.2;
#' S8, §10.3).** UST/IFX/ADA each now carry TWO rows in the source file, distinguished by a
#' `component`/`pricing_basis` pair rather than a single fixed HCPCS code: `"biosimilar_2026"` (the
#' default, and now the base case -- current biosimilar-inclusive pricing, re-extracted at a single
#' stated date) and `"originator_pre_biosimilar"` (this project's exact pre-fix figures, retained
#' verbatim, unchanged -- the comparability-with-Aliyev scenario S8 inverts to). Pass
#' `pricing_basis = "originator_pre_biosimilar"` to run that scenario; the two-argument
#' `price_for_component()` closure below is what actually branches on it. Cyclophosphamide and the
#' IV administration fee are unrelated products (not repriced by this fix) and keep the original
#' fixed-HCPCS-code lookup (`price_for_code()`) -- a code typo there still fails loudly rather than
#' silently returning NA, same discipline as before.
load_drug_prices <- function(raw_dir = "data/raw", pricing_basis = "biosimilar_2026") {
  df <- utils::read.csv(file.path(raw_dir, "cms_asp_and_hcup_cost_sources.csv"), stringsAsFactors = FALSE)

  price_for_code <- function(code) {
    val <- df$payment_limit_usd_per_unit[df$hcpcs_code == code]
    stopifnot(length(val) == 1, !is.na(val))
    val
  }
  price_for_component <- function(component_name) {
    rows <- df[df$component == component_name & df$pricing_basis == pricing_basis, ]
    stopifnot(nrow(rows) == 1, !is.na(rows$payment_limit_usd_per_unit))
    rows$payment_limit_usd_per_unit
  }

  list(
    ust_induction_usd_per_mg = price_for_component("UST induction"),
    ust_maintenance_usd_per_mg = price_for_component("UST maintenance"),
    ifx_usd_per_mg = price_for_component("IFX"),
    ada_usd_per_mg = price_for_component("ADA"),
    cyclophosphamide_usd_per_mg = price_for_code("J9075"),
    iv_administration_usd = price_for_code("96365")
  )
}

#' $/mg price for one therapy, from the `prices` list load_drug_prices() returns. UST has two
#' different prices (induction vs. maintenance -- different HCPCS codes/formulations); IFX and
#' ADA each have one price covering both phases.
drug_price_per_mg <- function(therapy, phase, prices) {
  switch(therapy,
    UST = if (phase == "Induction") prices$ust_induction_usd_per_mg else prices$ust_maintenance_usd_per_mg,
    IFX = prices$ifx_usd_per_mg,
    ADA = prices$ada_usd_per_mg,
    stop("Unknown therapy: ", therapy)
  )
}

#' UST's induction dose is a discrete weight-tiered label dose, already vial-exact (260/390/520mg
#' are exact 2x/3x/4x multiples of the 130mg vial Aliyev prices induction against) -- a tier
#' lookup against the sourced schedule, not a formula. Tier boundaries read from `schedule`
#' rather than hardcoded here, so drift between this function and the sourced CSV is impossible.
ust_induction_dose_mg <- function(weight_kg, schedule) {
  rows <- schedule[schedule$therapy == "UST" & schedule$phase == "Induction", ]
  lower_ok <- is.na(rows$weight_tier_min_kg) | weight_kg > rows$weight_tier_min_kg
  upper_ok <- is.na(rows$weight_tier_max_kg) | weight_kg <= rows$weight_tier_max_kg
  match <- rows[lower_ok & upper_ok, ]
  stopifnot(nrow(match) == 1)
  match$dose_amount
}

#' IFX's dose is mg/kg, not vial-exact -- round UP to the nearest whole vial (see module header:
#' Aliyev's own "no vial sharing" assumption). Applied identically to induction and maintenance
#' doses (both are single administration events).
ifx_dose_mg <- function(weight_kg, schedule, vial_size_mg = IFX_VIAL_SIZE_MG) {
  mg_per_kg <- unique(schedule$dose_amount[schedule$therapy == "IFX"])
  stopifnot(length(mg_per_kg) == 1)
  raw_mg <- weight_kg * mg_per_kg
  ceiling(raw_mg / vial_size_mg) * vial_size_mg
}

#' Dose size in mg for one schedule row, resolving IFX's weight-based mg/kg dosing to an actual
#' vial-rounded mg amount; every other therapy's dose_amount is already a flat/fixed mg figure
#' (UST's induction tiers, UST's 90mg maintenance, ADA's 160/80/40mg doses), read as-is.
dose_mg_for_row <- function(row, weight_kg, schedule) {
  if (row$therapy == "IFX") ifx_dose_mg(weight_kg, schedule) else row$dose_amount
}

#' One-time induction-phase drug cost (analysis_plan.md §6.1's 8-week decision-tree induction) --
#' NOT part of any per-cycle Markov trace; the caller (eventually R/05_deterministic_results.R)
#' adds this once per patient/arm on top of the per-cycle results below.
#'
#' UST's induction is a single weight-tiered dose -- one row selected by
#' ust_induction_dose_mg()'s tier lookup, not summed with the schedule's other UST induction rows
#' (those are alternative weight tiers, not sequential doses). IFX and ADA sum EVERY induction row
#' for that therapy: IFX has one row with n_doses=3 (three identical vial-rounded doses at weeks
#' 0/2/6); ADA has two rows with n_doses=1 each (160mg then 80mg, unequal sequential doses -- see
#' module header on why these are separate rows rather than one row x n_doses like IFX). Each
#' dose's administration fee is added only if that row's route is IV (UST induction, IFX); ADA is
#' always SC, so no administration fee at any point.
induction_drug_cost <- function(therapy, weight_kg, schedule, prices) {
  therapy <- match.arg(therapy, c("UST", "IFX", "ADA"))

  if (therapy == "UST") {
    dose_mg <- ust_induction_dose_mg(weight_kg, schedule)
    return(dose_mg * drug_price_per_mg("UST", "Induction", prices) + prices$iv_administration_usd)
  }

  rows <- schedule[schedule$therapy == therapy & schedule$phase == "Induction", ]
  stopifnot(nrow(rows) >= 1)
  price_per_mg <- drug_price_per_mg(therapy, "Induction", prices)
  total <- 0
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, ]
    admin_fee <- if (row$route == "IV") prices$iv_administration_usd else 0
    total <- total + row$n_doses * (dose_mg_for_row(row, weight_kg, schedule) * price_per_mg + admin_fee)
  }
  total
}

#' Per-cycle maintenance drug cost, added on top of on_biologic occupancy at the real dosing
#' cadence read from `schedule` (data/raw/biologic_dosing_schedule.csv) -- every 4th native
#' 2-week cycle for UST/IFX (q8w), every SINGLE cycle for ADA (EOW -- a materially different
#' cadence, not a typo, see module header) -- for as long as a patient remains on that arm's own
#' biologic track. `on_biologic_trace` mass naturally goes to 0 for patients who've switched to
#' CT or hit the 2-year cap (R/02_markov_engine.R), so this needs no separate cap-awareness.
#'
#' UST and ADA maintenance are SC (self-injection): drug cost only, no IV administration fee --
#' Aliyev's own Costs Assumption #3 (data/raw/aliyev2019_appendixS1_table1_assumptions.csv)
#' covers SC patients' visit costs via the monitoring/management cost, not a separate procedure
#' fee. IFX maintenance is IV: drug cost + administration fee every dose, same as induction.
maintenance_drug_cost_by_cycle <- function(on_biologic_trace, therapy, weight_kg, schedule, prices) {
  therapy <- match.arg(therapy, c("UST", "IFX", "ADA"))
  row <- schedule[schedule$therapy == therapy & schedule$phase == "Maintenance", ]
  stopifnot(nrow(row) == 1)

  n_cycles <- nrow(on_biologic_trace) - 1
  living <- setdiff(colnames(on_biologic_trace), "Death")
  mass_on_biologic <- rowSums(on_biologic_trace[, living, drop = FALSE])

  is_dose_cycle <- rep(FALSE, n_cycles + 1)
  dose_cycles <- seq(row$first_maintenance_cycle_native_2wk, n_cycles,
                      by = row$maintenance_interval_cycles_native_2wk)
  is_dose_cycle[dose_cycles + 1] <- TRUE  # +1: trace row 1 is cycle 0

  admin_fee <- if (row$route == "IV") prices$iv_administration_usd else 0
  per_dose_cost <- dose_mg_for_row(row, weight_kg, schedule) * drug_price_per_mg(therapy, "Maintenance", prices) + admin_fee

  mass_on_biologic * is_dose_cycle * per_dose_cost
}

# ---- Treg's own dose cost (data/processed/model_tenham_derived_treg_dose_cost.csv) -------------

#' Treg's own one-time dose acquisition cost (ten Ham et al. 2020 manufacturing costing
#' framework; analysis_plan.md §7.1 item 12, "Derived"). Reads the derivation file's own final
#' retail-price line directly -- see module header on why this deliberately does NOT read
#' model_dose_costs_and_psa_ranges.csv's "TREG dose cost" `value_usd` column instead (a real,
#' previously-unflagged ~0.9494x discrepancy found while wiring this in).
load_treg_dose_acquisition_cost <- function(proc_dir = "data/processed") {
  lines <- readLines(file.path(proc_dir, "model_tenham_derived_treg_dose_cost.csv"))
  blank <- which(lines == "")
  stopifnot(length(blank) >= 1)
  second_block <- utils::read.csv(
    text = paste(lines[(blank[1] + 1):length(lines)], collapse = "\n"), stringsAsFactors = FALSE
  )
  val <- second_block$value_usd[
    second_block$final_buildup_item == "Retail price per treatment/dose (used as cost_treg_dose base case)"
  ]
  stopifnot(length(val) == 1)
  as.numeric(val)
}

#' Treg's one-time per-dose cost: acquisition (ten Ham-derived, fully traceable) + infusion
#' administration (same $57.90 CMS PFS rate as UST/IFX/ADA's IV doses -- a 30-minute infusion
#' doesn't change this fee, see module header) + cyclophosphamide preconditioning (price sourced,
#' dose is NOT -- `cyclophosphamide_dose_mg` defaults to 0, a known-incomplete placeholder, not a
#' clinical claim that preconditioning is free or unnecessary; module header explains why a real
#' number isn't filled in here) + an optional overnight observation-stay cost
#' (`observation_stay_cost_usd`, defaults to 0 -- module header's 2026-08-04 sourcing pass gives
#' $2,672.15 (CY2026, confirmed directly from CMS's own Addendum A -- module header) as the
#' sourced figure to supply explicitly for a scenario run; the base case stays at 0, unchanged,
#' since whether a Treg infusion actually qualifies for this billing category at all is a
#' separate, still-open clinical/modelling question, not a pricing one).
#'
#' `n_doses` must be 1 (the recorded single-infusion base case, analysis_plan.md §4.1/§9.1) -- the
#' 2-dose structural scenario (S6, §10.3) needs per-cycle-trace-level logic this one-time function
#' doesn't have (see module header), so is deliberately rejected here rather than silently
#' computed wrong as `n_doses x this function's per-dose output`.
treg_dose_cost <- function(n_doses = 1, cyclophosphamide_dose_mg = 0, observation_stay_cost_usd = 0,
                            proc_dir = "data/processed", raw_dir = "data/raw", prices = NULL) {
  stopifnot(n_doses == 1)
  if (is.null(prices)) prices <- load_drug_prices(raw_dir)
  acquisition <- load_treg_dose_acquisition_cost(proc_dir)
  cyclophosphamide_cost <- cyclophosphamide_dose_mg * prices$cyclophosphamide_usd_per_mg
  n_doses * (acquisition + prices$iv_administration_usd + cyclophosphamide_cost
             + observation_stay_cost_usd)
}

# ---- Utilities (data/processed/model_psa_parameter_distributions.csv) -------------------------
#
# A16 fix (2026-08-05, docs/model_audit_v6.md): the deterministic utility vector used to be read
# straight from model_health_utilities.csv, whose Remission figure (0.9554396356) turned out to be
# a stray live PSA draw captured in the workbook snapshot, not a parameter -- four independent
# lines of evidence (full-double-precision reproduction of the chain below off that one number;
# every one of the five chained inputs sitting at a scattered, non-base percentile of its own
# stated Uniform range; ~16 significant figures on all five, the signature of volatile
# `=lo+RAND()*(hi-lo)` cells; and the fact that only 0.82, not 0.9554, reproduces the manuscript's
# own published utility table -- see the resolution memo folded into A16 below and
# docs/analysis_plan.md §15). model_health_utilities.csv is retired as a base-case input as a
# result -- kept in the repo, unused, with a provenance note in its own header row and in
# data/data_dictionary.md, not deleted (Data Availability Statement expectations: what a number
# WAS, not just what it now is, should stay reconstructable).
#
# The deterministic vector below is now DERIVED, not read: Remission's base value (0.82) and the
# four multiplicative ratios come from model_psa_parameter_distributions.csv's own `value` column
# -- the same file, and the same Beta-re-parameterised-from-Uniform base values (analysis_plan.md
# §7.1 items 19-20), that R/06_psa.R's sample_utility_chain_draws() already draws its Uniform
# bounds from. Chaining here uses the identical arithmetic (Mild = Remission x ratio, M-SR = Mild
# x ratio, M-S = M-SR x ratio, Surgery = M-SR x ratio) so the deterministic base case is, by
# construction, the central draw of the PSA rather than a separately maintained figure that can
# drift from it -- closing the actual root cause, not just this one instance of it (see
# test-parameter-provenance.R).

#' Named vector of per-state utility values (deterministic base case), derived from
#' model_psa_parameter_distributions.csv's sourced Remission base value and the four chained
#' ratios -- see module section header above for why this is a derivation, not a file read, as of
#' the A16 fix. Death is always 0 (not a row in the PSA file; no health-related utility applies).
load_health_state_utilities <- function(proc_dir = "data/processed") {
  df <- utils::read.csv(file.path(proc_dir, "model_psa_parameter_distributions.csv"), stringsAsFactors = FALSE)

  #' One row's base `value`, by its `variable` label -- a code typo fails loudly (stopifnot) rather
  #' than silently returning NA, same discipline as R/04's own load_drug_prices()/price_for().
  value_for <- function(variable_name) {
    val <- df$value[df$variable == variable_name]
    stopifnot(length(val) == 1, !is.na(val))
    val
  }

  remission <- value_for("Remission")
  mild <- remission * value_for("Mild:Remission Ratio")
  msr <- mild * value_for("M-SR:Mild Ratio")
  ms <- msr * value_for("M-S:M-SR Ratio")
  surgery <- msr * value_for("Surgery:M-SR Ratio")

  c(
    "Moderate-Severe" = ms, "Moderate-Severe Responder" = msr, "Mild" = mild,
    "Remission" = remission, "Surgery" = surgery, "Death" = 0
  )
}

# ---- Discounted per-trace QALYs and non-drug costs ---------------------------------------------

#' Per-cycle discount factors for a trace of length n_cycles+1 (cycles 0..n_cycles), reusing
#' R/02_markov_engine.R's discount_factor() rather than recomputing the formula here.
discount_factors_for_trace <- function(n_cycles, cycle_weeks = 2, annual_rate = 0.03) {
  vapply(0:n_cycles, discount_factor, numeric(1), cycle_weeks = cycle_weeks, annual_rate = annual_rate)
}

#' Half-cycle correction weights (analysis_plan.md §4.1: "Apply (life-table or Simpson's 1/3)" --
#' CHEERS 2022 item 17; closes A12, docs/model_audit_v6.md: "Cohort state occupancy is applied at
#' cycle boundaries with no correction"). The standard trapezoidal-rule correction for a trace of
#' `n_cycles` + 1 timepoints (rows 0..n_cycles, i.e. `n_cycles` actual cycle-length intervals
#' between them): true within-cycle occupancy is continuous and unobserved, so each interval's
#' contribution is approximated by the straight-line average of its two bounding timepoints,
#' (occ[i-1] + occ[i]) / 2 for i = 1..n_cycles. Summed algebraically that is
#' 0.5*occ[0] + occ[1] + ... + occ[n_cycles-1] + 0.5*occ[n_cycles] -- i.e. every INTERIOR
#' timepoint at full weight, both ENDPOINTS (baseline and terminal occupancy) at half weight.
#' `n_cycles == 0` (a single row, no interval elapsed at all) returns weight 1 -- there is nothing
#' to correct.
half_cycle_weights <- function(n_cycles) {
  stopifnot(n_cycles >= 0)
  if (n_cycles == 0) return(1)
  c(0.5, rep(1, n_cycles - 1), 0.5)
}

#' Discounted per-cycle QALYs for a full occupancy trace (rows = cycles 0..n, columns = states).
#' Undiscounted per-cycle utility mass is trace %*% utilities[states] (matches R/02's own
#' trace[t, ] %*% m style), multiplied by cycle length in years, the half-cycle correction weight
#' for that row, and the cycle's discount factor.
#'
#' `half_cycle_correction = TRUE` (default, closing A12 as of 2026-08-05 -- see
#' half_cycle_weights()): every result this function feeds now applies it by default, matching
#' analysis_plan.md §4.1's "Apply" recommendation. Pass `FALSE` to reproduce the exact
#' pre-2026-08-05 uncorrected behaviour (e.g. for a comparability scenario, or a hand-computed
#' test that wants the raw per-row mass) -- every prior deterministic/probabilistic result in this
#' repository was produced with the equivalent of `FALSE`, and is superseded by the `TRUE` default
#' now that it exists (README.md's Status section).
trace_qalys <- function(trace, utilities, cycle_weeks = 2, annual_rate = 0.03,
                         half_cycle_correction = TRUE) {
  states <- colnames(trace)
  stopifnot(!is.null(states), all(states %in% names(utilities)))
  n_cycles <- nrow(trace) - 1
  utility_mass <- as.numeric(trace %*% utilities[states])
  if (half_cycle_correction) utility_mass <- utility_mass * half_cycle_weights(n_cycles)
  utility_mass * cycle_weeks / 52 * discount_factors_for_trace(n_cycles, cycle_weeks, annual_rate)
}

#' Discounted per-cycle non-drug cost for a full occupancy trace. `add_ct_drug_cost = TRUE` adds
#' the flat CT-track drug cost (see module header) to every living state's occupancy each cycle --
#' pass this for a CT-track trace (e.g. run_maintenance_arm()'s on_ct), not a biologic-track one.
#'
#' `half_cycle_correction` -- see trace_qalys()'s docstring; same default, same weights, same
#' supersession note. Applied to the COMBINED cost mass (monitoring cost plus, if requested, the
#' CT drug cost) so both components of a cycle's cost get the same within-cycle treatment, not
#' just the state-occupancy-driven part.
trace_costs <- function(trace, monitoring_costs, cycle_weeks = 2, annual_rate = 0.03,
                         add_ct_drug_cost = FALSE, half_cycle_correction = TRUE) {
  states <- colnames(trace)
  stopifnot(!is.null(states), all(states %in% names(monitoring_costs)))
  n_cycles <- nrow(trace) - 1
  cost_mass <- as.numeric(trace %*% monitoring_costs[states])
  if (add_ct_drug_cost) {
    living <- setdiff(states, "Death")
    cost_mass <- cost_mass + rowSums(trace[, living, drop = FALSE]) * ct_drug_cost()
  }
  if (half_cycle_correction) cost_mass <- cost_mass * half_cycle_weights(n_cycles)
  cost_mass * discount_factors_for_trace(n_cycles, cycle_weeks, annual_rate)
}

# ---- Whole-arm attachment ------------------------------------------------------------------

#' Attach costs and QALYs to one run_maintenance_arm()-style result (on_biologic, on_ct, total;
#' R/02_markov_engine.R). Non-drug costs: on_ct gets the CT drug cost added, on_biologic doesn't
#' (its own drug cost, if any, is the separate `drug_cost_by_cycle` field below). QALYs use
#' on_biologic + on_ct (equivalent to `total` for a plain maintenance arm) since utility depends
#' only on health state, not which track a patient is on.
#'
#' `therapy` activates the UST/IFX/ADA drug-cost layer (module header) -- pass "UST", "IFX", or
#' "ADA" plus `weight_kg`/`schedule`/`prices` to charge on_biologic's own maintenance drug cost;
#' `drug_cost_by_cycle` stays all-zero if `therapy` is NULL (the default), which is deliberate for
#' CT-only arms (CT's drug cost is already the separate ct_drug_cost() layer on on_ct, not this
#' one) and for Treg's non-cured track (see attach_treg_costs_utilities() -- never passes
#' `therapy` here, by design, not by omission).
#'
#' `half_cycle_correction` (default TRUE, A12/analysis_plan.md §4.1) -- passed through to the
#' non-drug (state-occupancy-driven) cost and QALY calls only, NOT to `drug_cost_by_cycle` below:
#' a maintenance dose is a discrete, punctate event that either happens in a given cycle or
#' doesn't (`maintenance_drug_cost_by_cycle()`'s own `is_dose_cycle` mask), not a continuously
#' accruing quantity like time spent occupying a health state -- the correction exists to
#' approximate the latter's within-cycle path, and would misapply to the former (e.g. halving the
#' cost of a dose actually given at cycle 0 or the final cycle, which isn't the mechanism A12
#' describes at all).
attach_maintenance_costs_utilities <- function(arm_result, utilities, monitoring_costs,
                                                cycle_weeks = 2, annual_rate = 0.03,
                                                therapy = NULL, weight_kg = NULL,
                                                schedule = NULL, prices = NULL,
                                                half_cycle_correction = TRUE) {
  states_trace <- arm_result$on_biologic + arm_result$on_ct
  non_drug_cost_by_cycle <- (
    trace_costs(arm_result$on_biologic, monitoring_costs, cycle_weeks, annual_rate,
                add_ct_drug_cost = FALSE, half_cycle_correction = half_cycle_correction) +
    trace_costs(arm_result$on_ct, monitoring_costs, cycle_weeks, annual_rate,
                add_ct_drug_cost = TRUE, half_cycle_correction = half_cycle_correction)
  )

  n_cycles <- nrow(arm_result$on_biologic) - 1
  if (!is.null(therapy)) {
    stopifnot(!is.null(weight_kg), !is.null(schedule), !is.null(prices))
    drug_cost_by_cycle <- (
      maintenance_drug_cost_by_cycle(arm_result$on_biologic, therapy, weight_kg, schedule, prices) *
      discount_factors_for_trace(n_cycles, cycle_weeks, annual_rate)
    )
  } else {
    drug_cost_by_cycle <- rep(0, n_cycles + 1)
  }

  list(
    qalys_by_cycle = trace_qalys(states_trace, utilities, cycle_weeks, annual_rate,
                                  half_cycle_correction = half_cycle_correction),
    non_drug_cost_by_cycle = non_drug_cost_by_cycle,
    drug_cost_by_cycle = drug_cost_by_cycle,
    total_cost_by_cycle = non_drug_cost_by_cycle + drug_cost_by_cycle
  )
}

#' SDR-specific cost and QALYs (analysis_plan.md §6.2). `on_sdr` is a plain numeric vector
#' (cycles 0..n), not a states matrix -- R/03_cure_fraction_module.R's own representation (SDR
#' has no internal state structure). Utility = Remission's utility by default (the base case);
#' Cost = Remission's monitoring cost, halved after `halve_after_cycle` (default the 2-year cap
#' boundary, cycle 52 at this project's native 2-week cycle) -- explicitly "recommend...flag as
#' assumption" in the plan text (§6.2), not an independently sourced figure. `drug_cost_by_cycle`
#' is always all-zero: SDR patients are off therapy by definition, not an omission to fill in
#' later.
#'
#' `sdr_utility = NULL` (default): use `utilities[["Remission"]]`, exactly as before this
#' parameter existed -- fully backward compatible. Scenario S7 (§7.1 item 21, §10.3) passes a
#' general-population utility value/schedule instead (R/utils/population_utility.R) -- either a
#' single scalar (recycled across every cycle, R's own vector recycling, the right behaviour for
#' the 6.15-year/10-year horizons which don't track attained age) or a length-(n_cycles+1) vector
#' (a genuinely age-varying schedule, for the lifetime horizon). No shape-checking needed here:
#' ordinary R recycling of a length-1 or length-(n_cycles+1) vector against `on_sdr` does the
#' right thing either way; anything else would already be a caller error the multiplication itself
#' would surface via a recycling-length warning.
#'
#' `half_cycle_correction` (default TRUE) -- SDR occupancy is exactly the same kind of
#' continuously-accruing quantity as any other health state's (patients spend time IN it, same as
#' Remission or Mild), so it gets the identical trapezoidal treatment as
#' trace_qalys()/trace_costs() apply elsewhere, via the same half_cycle_weights() helper -- not a
#' separate rule invented for this function.
attach_sdr_costs_utilities <- function(on_sdr, utilities, monitoring_costs, cycle_weeks = 2,
                                        annual_rate = 0.03, halve_after_cycle = 52,
                                        half_cycle_correction = TRUE, sdr_utility = NULL) {
  n_cycles <- length(on_sdr) - 1
  cycles <- 0:n_cycles
  discount <- discount_factors_for_trace(n_cycles, cycle_weeks, annual_rate)
  weights <- if (half_cycle_correction) half_cycle_weights(n_cycles) else 1
  if (is.null(sdr_utility)) sdr_utility <- utilities[["Remission"]]
  remission_cost <- monitoring_costs[["Remission"]]
  cost_rate <- ifelse(cycles > halve_after_cycle, remission_cost / 2, remission_cost)
  non_drug_cost_by_cycle <- on_sdr * cost_rate * weights * discount

  list(
    qalys_by_cycle = on_sdr * sdr_utility * cycle_weeks / 52 * weights * discount,
    non_drug_cost_by_cycle = non_drug_cost_by_cycle,
    drug_cost_by_cycle = rep(0, n_cycles + 1),
    total_cost_by_cycle = non_drug_cost_by_cycle
  )
}

#' Attach costs and QALYs to a run_treg_arm()-style result (on_biologic, on_ct, on_sdr, total;
#' R/03_cure_fraction_module.R). Combines attach_maintenance_costs_utilities()'s logic for the
#' Markov portion with attach_sdr_costs_utilities()'s rule for the SDR portion -- arm_result$total
#' isn't used directly here since it's already collapsed to a headcount scalar per cycle
#' (rowSums(on_biologic) + rowSums(on_ct) + on_sdr) and can't be re-split by state.
#'
#' Deliberately never passes `therapy` to the Markov portion below -- see module header ("Not
#' applied to Treg's non-cured track"): non-cured Treg patients are efficacy-equivalent to UST
#' but were never actually given ustekinumab, so they are not charged its drug cost.
#'
#' `half_cycle_correction` (default TRUE) -- forwarded unchanged to both the Markov portion and
#' the SDR portion, so the whole arm gets one consistent within-cycle treatment.
#'
#' `sdr_utility = NULL` (default) -- forwarded unchanged to attach_sdr_costs_utilities() (see its
#' own docstring for scenario S7); the Markov (non-SDR) portion is unaffected either way, since
#' the general-population-utility scenario is specifically about the SDR state's own utility.
attach_treg_costs_utilities <- function(arm_result, utilities, monitoring_costs, cycle_weeks = 2,
                                         annual_rate = 0.03, halve_after_cycle = 52,
                                         half_cycle_correction = TRUE, sdr_utility = NULL) {
  markov <- attach_maintenance_costs_utilities(
    list(on_biologic = arm_result$on_biologic, on_ct = arm_result$on_ct),
    utilities, monitoring_costs, cycle_weeks, annual_rate,
    half_cycle_correction = half_cycle_correction
  )
  sdr <- attach_sdr_costs_utilities(
    arm_result$on_sdr, utilities, monitoring_costs, cycle_weeks, annual_rate, halve_after_cycle,
    half_cycle_correction = half_cycle_correction, sdr_utility = sdr_utility
  )
  list(
    qalys_by_cycle = markov$qalys_by_cycle + sdr$qalys_by_cycle,
    non_drug_cost_by_cycle = markov$non_drug_cost_by_cycle + sdr$non_drug_cost_by_cycle,
    drug_cost_by_cycle = markov$drug_cost_by_cycle + sdr$drug_cost_by_cycle,
    total_cost_by_cycle = markov$total_cost_by_cycle + sdr$total_cost_by_cycle
  )
}

#' Lifetime (sum of discounted cycles) QALYs and cost for one arm's attach_*() output.
#' `induction_cost` (default 0) is added in undiscounted -- it's already a cycle-0 dollar amount
#' from induction_drug_cost(), and discount_factor(0) = 1 regardless -- for whichever arm's
#' one-time induction drug cost the caller has separately computed; this function only sums
#' per-cycle series, it doesn't call induction_drug_cost() itself.
summarise_arm <- function(attached, induction_cost = 0) {
  list(
    qalys = sum(attached$qalys_by_cycle),
    non_drug_cost = sum(attached$non_drug_cost_by_cycle),
    drug_cost = sum(attached$drug_cost_by_cycle) + induction_cost,
    total_cost = sum(attached$total_cost_by_cycle) + induction_cost
  )
}
