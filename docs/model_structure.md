# Model Structure

Technical specification of the decision-tree + Markov + mixture-cure engine **as actually
implemented** in `R/00_derive_transition_probs.R` through `R/09_scenarios.R`, kept in sync with
the code rather than with the design rationale — that's `docs/analysis_plan.md` §6's job (the
proposal, the alternatives considered, and why each choice was made). This file exists to answer
"what does the code do," not "why," and is written last (2026-08-05), after the engine was
essentially complete, specifically so it reflects the code rather than something that had to be
rewritten every time the engine changed. If this file and the code ever disagree, the code is
right and this file is stale — treat that as a documentation bug, not a design question.

---

## 1. Health states

Six states, shared by every arm's maintenance-phase matrix (`MAINTENANCE_STATES`,
`R/00_derive_transition_probs.R`):

| State | Meaning | Notes |
|---|---|---|
| Moderate-Severe (M-S) | Active, non-responding disease | For UST/IFX/ADA, this state has no live transition row of its own — see §3's CT-switch rule |
| Moderate-Severe Responder (M-SR) | Responding but not in Mild/Remission | Aliyev's own intermediate response category |
| Mild | Mild disease activity | |
| Remission | CDAI-defined remission | The best non-cure state; utility/cost anchor for several downstream calculations |
| Surgery | Post-surgical state | Its own row has an unusually high one-step return probability to Remission in Aliyev's published matrix — see §8 |
| Death | Absorbing | |

Treg's arm adds a seventh state that exists **only** in `R/03_cure_fraction_module.R`'s own
representation, not in `MAINTENANCE_STATES`:

| State | Meaning | Notes |
|---|---|---|
| Sustained Deep Remission (SDR) | Durable, drug-free cure | A scalar per cycle (no internal sub-structure — a patient is either in SDR or not), not a column of the shared state matrix. See §6. |

---

## 2. Two linked components

```
INDUCTION (decision tree, one-time split)          MAINTENANCE (Markov, native 2-week cycle)
────────────────────────────────────────           ──────────────────────────────────────────
   Cohort, arm by arm                                 Biologic-arm cohort  <---CT-switch--->  CT-arm cohort
        │                                                     │                                    │
  Aliyev Table 3 (per therapy)                          Aliyev Table 4 (per therapy)          Aliyev Table 4 (CT)
        │                                                     │                                    │
  ┌─────┴─────┐                                     M-S / M-SR / Mild / Remission / Surgery / Death, every cycle
Responder    Non-responder
(M-SR/Mild/    (enters CT
 Remission)     directly)

TREG ARM ONLY, grafted onto the UST pipeline above:
  week-56 landmark (cycle 28) --split responders in Remission--> Sustained Deep Remission (cured)
                                                              \-> continue on UST-equivalent Markov (non-cured)
```

Every comparator arm (UST, IFX, ADA) runs induction once, then maintenance for the rest of the
horizon. CT is not an independent comparator arm — it is the fallback track every biologic arm's
non-responders and treatment failures are routed into (§3). Treg reuses UST's own induction split
and maintenance matrix for everything except the cure mechanic itself (§6) — there is no
independently-sourced Treg transition matrix anywhere in this codebase.

---

## 3. Induction (decision tree)

**Implementation:** `R/01_decision_tree.R`'s `run_decision_tree()`, fed by
`R/00_derive_transition_probs.R`'s `load_published_induction()`.

**Source:** Aliyev et al. 2019 Appendix S2, Supplementary Table 3, transcribed verbatim
(`data/raw/aliyev2019_appendixS2_table3_induction_transition_probabilities.csv`). Table 3's
Moderate-Severe row is a **per-2-week-cycle** transition probability derived from a week-6 (UST)
or week-4 (IFX/ADA) trial endpoint (Appendix S2's own worked example) — **not** the terminal
end-of-induction split itself (B1 fix, 2026-08-06, `docs/model_audit_v6.md` A17; earlier
revisions of this project read it as a one-step split, which was wrong and understated every
comparator's true responder pool by ~1.75–2.4×). `load_published_induction()`
(`R/00_derive_transition_probs.R`) now runs the full induction matrix through `simulate_cohort()`
for `INDUCTION_CYCLES[[therapy]]` cycles (UST = 3, IFX/ADA = 2) and reads the terminal row, which
is what feeds `run_decision_tree()` — that function's own signature/mechanics are unchanged by
the fix. Each therapy's corrected terminal row gives four terminal shares: `to_moderate_severe`,
`to_moderate_severe_responder`, `to_mild`, `to_remission`. UST/IFX/ADA each have their own row;
**CT has no induction row at all** — patients only ever reach CT by failing a biologic first
(§3), never by starting there.

**Mechanics:** `to_moderate_severe_responder`/`to_mild`/`to_remission` become the initial
occupancy vector for the biologic-arm maintenance trace; `to_moderate_severe` becomes the initial
occupancy for the CT-arm trace (patients who never responded start directly on CT). Both vectors
are normalised so their combined total is exactly 1 — a mechanical correction for Aliyev's own
publication-rounding (row sums off by up to ~0.1%, e.g. Table 3's induction rows sum to 0.9994,
not 1.0000 exactly), not a modelling assumption.

Scenario S3 (refractory population, §9) and S12 (non-cured HR-advantage, §6) both operate on
copies of this induction row / the maintenance matrix it feeds — never on the row Aliyev actually
published, which stays available unmodified for every other caller.

---

## 4. Maintenance (cohort Markov)

**Implementation:** `R/02_markov_engine.R`'s `run_maintenance_arm()` /
`run_maintenance_arm_with_mortality()`, fed by `load_published_maintenance()`
(`R/00_derive_transition_probs.R`) and `build_transition_matrix()`
(`R/utils/transition_matrix.R`).

**Cycle length:** native 2 weeks, matching Aliyev's own Supplementary Table 4 exactly — no
cycle-length conversion anywhere in this codebase (an earlier 8-week design, chosen to align with
UST/IFX's q8w dosing cadence, was dropped: drug costs attach at whatever cycle real dosing
happens regardless of the transition-probability cycle length, so the alignment bought nothing,
and running Aliyev's matrices unmodified at their native cycle directly strengthens external
validation against his own published results).

**Matrix construction:** `build_transition_matrix()` builds each therapy's square
state-by-state matrix from Aliyev's long-format table, padding any state with no reported row as
an absorbing self-loop (probability 1) purely so the matrix is square and every state has *some*
row — not a modelling claim. UST/IFX/ADA's own matrices have no live `Moderate-Severe` row for
exactly this reason (see below); rows are then renormalised to sum to exactly 1 (correcting the
same ~0.1%-level publication rounding noted in §3, which would otherwise compound
multiplicatively every cycle over a long horizon).

**The CT-switch rule** (`R/02_markov_engine.R`'s internal step function, `ms_idx`/
`apply_cap_now` logic): every cycle, whatever occupancy a biologic arm's own matrix would have
placed in Moderate-Severe is instead swept entirely into the CT-arm's own occupancy at that same
cycle (`ct_next[ms_idx] <- ct_next[ms_idx] + bio_next[ms_idx]; bio_next[ms_idx] <- 0`) — "flowing
into conventional therapy" represents exhausting a maintenance drug, mechanistically meaningless
for a patient already on CT, so this switch only ever moves mass ONE way, biologic → CT, never
back. This is the live implementation of the padded-self-loop matrices noted above: the
Moderate-Severe row those matrices carry is never actually consulted for a biologic arm in a real
simulation, since occupancy is swept out before the next cycle's transition would use it.

**The 2-year cap** (`cap_cycle = 52` at the native 2-week cycle, `apply_cap = TRUE` default):
on the specific cycle `t == cap_cycle`, in addition to the ordinary per-cycle CT-switch above,
*all* remaining biologic-arm occupancy (not just the Moderate-Severe share) is moved to CT in one
step — modelling a maintenance drug's assumed maximum duration of use. `apply_cap = FALSE`
reproduces the no-cap structural scenario (S1, §9) for free, with no separate code path.

**Discounting:** `discount_factor(cycle, cycle_weeks, annual_rate)` —
`1 / (1 + annual_rate)^(cycle * cycle_weeks / 52)`, i.e. ordinary compound annual discounting
converted to whatever fraction of a year the cycle length represents. `annual_rate = 0.03` is the
base case; scenario S10 (§9) sweeps 0%/1.5%/5% through the same parameter, no separate code path.

**Half-cycle correction** (`half_cycle_weights(n_cycles)`, standard trapezoidal method): both
trace endpoints (cycle 0 and the final cycle) get weight 0.5, every interior cycle gets weight 1
— applied to continuously-accruing, state-occupancy-driven cost/QALY calculations only, never to
one-time/discrete events (a maintenance dose actually given on a specific cycle, Treg's own
acquisition cost). `half_cycle_correction = TRUE` is the default from `R/04_costs_utilities.R`
through every downstream caller; every PSA/EVPI/EVPPI/EJP result inherits it automatically.

**Lifetime-horizon mortality** (`run_maintenance_arm_with_mortality()`,
`R/utils/life_table.R`/`age_adjust_matrix()` in `R/utils/transition_matrix.R`): Aliyev's own
matrices carry a flat, non-age-varying ~0.00006-per-cycle all-cause death probability, realistic
for his actual (young, short-follow-up) trial cohort but wrong once the horizon runs into old
age. `age_adjust_matrix()` **replaces** (not adds to) each row's Death-column entry with the
US-life-table figure (NCHS 2021, age- and sex-specific) for the cohort's current attained age,
rescaling that row's other entries proportionally so it still sums to 1 — the "no CD excess
mortality" assumption `docs/analysis_plan.md` §7.1 item 7 specifies (total mortality =
general-population background only). Male and female sub-cohorts (50/50) run separately against
their own life-table curve and are summed — exact by linearity, not an approximation, since both
sub-cohorts run through structurally identical matrices differing only in the Death column.

**Horizons actually run** (`R/05_deterministic_results.R`): `HORIZON_CYCLES_6YR` = 160 cycles
(the current-draft/primary reporting horizon, Aliyev's own embedded mortality, no age-tracking);
`HORIZON_CYCLES_10YR` = 260 cycles (same mortality treatment, a comparability scenario, S5);
`HORIZON_CYCLES_LIFETIME` = 1691 cycles (age-adjusted mortality from `ASSUMED_PATIENT_AGE_YEARS`
= 35 to the life table's terminal age 100, S5's lifetime arm). PSA/EVPI/EVPPI/probabilistic EJP
still run at the 6.15-year horizon only — extending age-adjusted mortality to a ~40,000-call PSA
is a real, deliberately deferred performance question, not attempted.

---

## 5. Structural population variants: refractory multipliers (Decision 3, scenario S3)

**Implementation:** `R/utils/refractory_multipliers.R`, applied inside
`run_comparator_arm_lifetime()` (`R/05_deterministic_results.R`) when `refractory = TRUE`.

Rather than re-deriving a refractory population's transition probabilities from first principles,
this project applies explicit, sourced multipliers to the biologic-naive matrices above:

- **Induction**: `apply_refractory_multiplier_induction()` shrinks total response and remission
  (sourced from UNITI-1 vs. UNITI-2, Feagan et al. 2016) separately, then rebuilds the
  response-but-not-remission share as the difference, floor-clamped so remission can never exceed
  total response.
- **Maintenance**: `apply_refractory_multiplier_maintenance()` shrinks each cycle's flow into
  Remission (converted from a 36-week cumulative IM-UNITI figure to a per-cycle multiplier) and
  redistributes the removed mass into Mild; an optional second step elevates each row's Surgery
  probability (a deliberately upper-bound proxy, sourced from a population-severity-mismatched
  pair — see that file's own header) via the same proportional-rescale mechanic `age_adjust_matrix()`
  uses for installing a new Death probability.

Applied to UST/IFX/ADA's own matrices only, on a **local copy** — never to the shared matrices
object every arm reads from, so the real (non-refractory) comparator runs in the same call are
unaffected. CT's matrix and Treg's own parameterisation are untouched (Treg's reference clinical
programme is already a refractory population by design, per `docs/analysis_plan.md` §5.1).

---

## 6. The mixture-cure extension (Treg arm)

**Implementation:** `R/03_cure_fraction_module.R`'s `run_treg_arm()` /
`run_treg_arm_with_mortality()`.

Treg's non-cured population is **efficacy-equivalent to ustekinumab** — it runs on UST's own
induction split and maintenance matrix, unmodified, with no separate Treg transition matrix
anywhere in this codebase (the "20% adverse-transition reduction" the original workbook applied
was dropped entirely; all of Treg's incremental benefit flows through π and h below, not through
a separately-adjusted efficacy figure).

**The landmark split.** At cycle 28 (week 56 — chosen as the nearest cycle boundary to the plan's
original 52-week target under its now-retired 8-week-cycle design, where it was cycle 7; simply
recomputed as cycle 28 = 56/2 for this project's current native 2-week cycle, not re-derived from
scratch), patients who are in Remission under the UST-equivalent trace are split: a fraction π
(`pi_sdr`) enters Sustained Deep Remission; the remainder continues on the ordinary UST-equivalent
Markov trace exactly as before. π has no
sourced value — Treg has no efficacy data — and is swept as a free parameter across its full
0–100% range, with named anchor points at 0/50/75/90% (scenario S4, `PI_ANCHOR_GRID`,
`R/09_scenarios.R`).

**Sustained Deep Remission (SDR).** A scalar occupancy per cycle, not a column of
`MAINTENANCE_STATES` — a patient in SDR has no internal sub-structure to be in. Each cycle, SDR
occupancy is split three ways: continues in SDR, relapses (competing against death first), or
dies. The relapse hazard h (`relapse_hazard_annual`, converted to a per-cycle probability the
same way `discount_factor()` converts an annual rate) has no sourced value either — parameterised
by its implied median duration T (`duration_to_hazard()`/`hazard_to_duration()`,
`RELAPSE_DURATION_GRID_YEARS` = {2, 5, 10, 20, ∞}, base case T = 10 years) and swept jointly with
π and price for the (π, T, price) headroom surface, since h interacts materially with a lifetime
horizon (negligibly at 6.15 years).

**Relapse re-entry.** A relapsed SDR patient re-enters the ordinary Markov trace at Mild (base
case) or Moderate-Severe Responder (scenario S11, `relapse_destination` parameter) — real but
small in effect (a few thousandths of a QALY at π=0.9; exactly zero at π=0, since no patient is
ever in SDR to relapse there).

**Cap interaction.** The 2-year cap (§4) applies identically to non-cured Treg responders as to
UST/IFX/ADA, specifically so Treg gets no unearned advantage over the comparators it's compared
against. A relapse that happens *after* the cap has already fired routes directly to CT rather
than back to the (no-longer-available) biologic track — consistent with the cap's own
"exhaustion of a maintenance drug" rationale.

**Non-cured HR-advantage (scenario S12).** `apply_non_cured_hazard_ratio()`
(`R/03_cure_fraction_module.R`) is the mechanism `docs/analysis_plan.md` §6.2 itself names for
retaining some non-cure advantage instead of full UST-equivalence: a hazard ratio on the
Moderate-Severe and Surgery transitions, via the discrete-time proportional-hazards transform
`p_new = 1 - (1 - p_old)^HR` (the correct translation of a hazard ratio into a per-cycle
probability — not a naive `p_old * HR`, which has no such guarantee of staying in [0,1]). Applied
to a **local copy** of UST's matrix, used only for Treg's own non-cured track; the real UST
comparator arm's matrix is never mutated. `HR = 1` is the exact identity transform (the §6.2 base
case). **A genuinely counterintuitive, verified-correct finding**: strengthening the advantage
(HR < 1) very slightly *lowers* aggregate QALYs, because Aliyev's own Surgery row has an 86.7%
one-step chance of landing back in Remission the very next cycle — faster than continuing
through Moderate-Severe Responder or Mild offers — so reducing Surgery entry isn't unambiguously
a QALY win in this specific matrix, even though it is unambiguously a cost win. NMB still rises
monotonically as the advantage strengthens at every HR and WTP threshold tested, since the cost
saving dominates the tiny QALY effect.

---

## 7. Structural variants of utility/cost/perspective

Two further scenarios modify how costs/utilities are attached, not the transition structure
itself:

**SDR utility source (scenario S7).** `attach_sdr_costs_utilities()`
(`R/04_costs_utilities.R`) uses Remission's own utility for SDR by default (§6.2's base case: a
durable cure is "as good as" ordinary remission, no more, no less). The general-population
alternative (`sdr_utility_source = "general_population"`, `R/05`'s `run_treg_arm_lifetime()`)
substitutes a sex-averaged US EQ-5D-5L population norm (Jiang et al. 2021,
`R/utils/population_utility.R`) — a single reference-age value for the 6.15-year/10-year
horizons, or a genuinely age-varying per-cycle schedule for the lifetime horizon. Population
utility at the cohort's baseline age (0.843) modestly exceeds Remission's own (0.82), so this
scenario raises Treg's QALYs slightly wherever any patient is actually cured.

**Perspective (scenario S9).** `health_state_monitoring_costs()` (healthcare-sector, the base
case) vs. `societal_monitoring_costs()` — the latter adds a flat per-cycle productivity/indirect
cost (Manceur et al. 2020's $2,168/year incremental work-loss cost per CD patient) to every
living state uniformly, identically across every arm. `perspective = "healthcare_sector"`
(default) vs. `"societal"` selects which `monitoring_costs` vector every cost-attachment call
uses — no other code path changes.

---

## 8. Cost and utility attachment

**Implementation:** `R/04_costs_utilities.R`'s `attach_maintenance_costs_utilities()` (UST/IFX/
ADA/CT), `attach_sdr_costs_utilities()` + `attach_treg_costs_utilities()` (Treg, combining the
Markov-portion logic above with the SDR-specific rule).

**Monitoring cost**: a flat per-cycle figure by health state (`ALIYEV_MONITORING_COST_2017_USD`,
inflated to 2025 USD), identical across every arm — Moderate-Severe and Moderate-Severe Responder
share a figure (Aliyev does not report M-SR separately), Death is always 0. Surgery's own
monitoring cost is markedly higher than every other live state's, which combines with its own
row's fast return-to-Remission property (§6) to produce the cost-saving/QALY-neutral pattern seen
in several scenarios above.

**Drug cost**: UST/IFX/ADA each have an induction (one-time, decision-tree-level) and maintenance
(recurring, dosed at each therapy's own cadence — q8w for UST/IFX, every-other-week for ADA, read
from `data/raw/biologic_dosing_schedule.csv` rather than hardcoded) drug cost, priced from CMS
Part B ASP or NADAC depending on which benefit pathway actually applies to that product's
administration route. CT has its own flat per-cycle drug cost, added on top of its own monitoring
cost for every living CT-track state. Treg's own one-time acquisition + administration +
(optional) cyclophosphamide + (optional) observation-stay cost is `treg_dose_cost()`
(base case) / `treg_price_dependent_dose_cost()` (price swept as a variable, headroom/EJP use).
Non-cured Treg patients are never charged UST's own drug cost, despite running on its transition
matrix — they were never actually given ustekinumab.

**SDR cost**: Remission's own monitoring cost, halved after the 2-year cap boundary (an explicit,
flagged assumption, not an independently sourced figure) — SDR patients accrue no drug cost at
all, by definition of being off therapy.

---

## 9. Cross-reference: structural scenarios (S1–S12) to the mechanism each modifies

| # | Scenario | Modifies | Implemented |
|---|---|---|---|
| S1 | 2-year cap on/off | §4 (`apply_cap`) | Yes |
| S2 | ADA in/out of comparator set | §2 (which arms run) | Yes |
| S3 | Biologic-naive vs. refractory | §3, §4 (induction/maintenance matrices) | Yes |
| S4 | Cure-fraction anchor (π × T grid) | §6 (π, h) | Yes |
| S5 | Time horizon | §4 (`n_cycles`, mortality treatment) | Yes |
| S6 | Treg 1 dose vs. 2 doses | §8 (drug cost timing) | **Dropped** — no dosing-interval anchor exists; not worth inventing one |
| S7 | SDR utility source | §7 | Yes |
| S8 | Comparator pricing basis | §8 (`load_drug_prices(pricing_basis = ...)`) | Yes |
| S9 | Perspective | §7 | Yes |
| S10 | Discount rate | §4 (`annual_rate`) | Yes |
| S11 | SDR relapse destination | §6 (`relapse_destination`) | Yes |
| S12 | Non-cured Treg HR-advantage | §6 (`apply_non_cured_hazard_ratio()`) | Yes |

---

## 10. Module index

| File | Responsibility |
|---|---|
| `R/00_derive_transition_probs.R` | Load/validate Aliyev's published induction and maintenance transition probabilities |
| `R/01_decision_tree.R` | Induction-phase decision tree (§3) |
| `R/02_markov_engine.R` | Cohort Markov maintenance engine, CT-switch, 2-year cap, discounting (§4) |
| `R/03_cure_fraction_module.R` | Mixture-cure extension: landmark split, SDR, relapse, non-cured HR-advantage (§6) |
| `R/04_costs_utilities.R` | Cost/utility attachment, half-cycle correction, SDR cost/utility rule, societal perspective (§7–8) |
| `R/05_deterministic_results.R` | Base case, headroom frontier (Aim 4), refractory scenario orchestration, horizon/pi/duration constants |
| `R/06_psa.R` | 10,000-draw probabilistic sensitivity analysis |
| `R/07_evpi_evppi.R` | EVPI/EVPPI by parameter subset (Aims 2–3) |
| `R/08_ejp.R` | Economically justifiable price (deterministic + probabilistic), gross margin over COGS |
| `R/09_scenarios.R` | Structural scenario orchestration (S1, S2, S4, S5, S7, S9, S10, S11, S12) |
| `R/utils/transition_matrix.R` | Shared matrix helpers: validation, construction, `age_adjust_matrix()` |
| `R/utils/life_table.R` | US life-table background mortality for the lifetime horizon |
| `R/utils/population_utility.R` | US general-population utility norms (S7) |
| `R/utils/refractory_multipliers.R` | UNITI-1-vs-UNITI-2-sourced refractory-population multipliers (S3) |

For parameter provenance (which CSV backs which figure, and why it's trusted as a parameter
rather than workbook snapshot state), see `data/data_dictionary.md` and
`tests/testthat/test-parameter-provenance.R`. For the design rationale behind each choice above —
the alternatives considered, the co-author decisions, and the open items — see
`docs/analysis_plan.md`.
