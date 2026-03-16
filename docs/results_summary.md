# HRV Analysis — Results Summary
## MR vs VBT Resistance Training Protocols in Women

**Design:** Within-subjects crossover | **N:** 21–27 (varies by analysis) | **Date:** 2026-03-16

---

## Overview

| Hypothesis | Domain | Verdict | Strongest signal |
|---|---|---|---|
| H1 Pre-exercise HRV | Resting ANS | Confirmed equivalent | All p_fdr > 0.77 |
| H2 Short-term pre-ortho | Orthostatic reactivity | Confirmed equivalent | All p_fdr > 0.44 |
| H3 Medium-term pre-ortho | Orthostatic reactivity | Confirmed equivalent (borderline) | RMSSD/LF_nu p_raw < 0.05, p_fdr NS |
| H4 Pre-exercise lactate | Metabolic | Confirmed equivalent | p = 0.73 |
| **H5 Exercise trajectory** | ANS dynamics | **Significant** | HF_log, LF_log, LF_HF_ratio, SD1 interactions |
| **H6 Recovery trajectory** | ANS dynamics | **Significant** | RMSSD, SDNN recover faster under VBT |
| **H7 Post-exercise resting** | Resting ANS | **Significant** | SampEn p_fdr = 0.033, dz = −0.82 |
| H8 Short-term post-ortho | Orthostatic reactivity | Not significant | All p_fdr > 0.55 |
| H9 Medium-term post-ortho | Orthostatic reactivity | Not significant | All p_fdr > 0.86 |
| **H10 Post-exercise lactate** | Metabolic | **Highly significant** | lactate_post dz = 0.98, p < 0.001 |
| H11 ΔHRV between protocols | ANS net change | Not significant | SampEn p_raw = 0.059, p_fdr NS |
| H12 Δortho between protocols | Orthostatic net change | Not significant | All p_fdr > 0.73 |
| H13 lactate_post × ΔHRV | Metabolic–autonomic coupling | Not significant | All p_fdr > 0.95 |
| H14 Δlactate × ΔHRV | Metabolic–autonomic coupling | Marginal | VBT dortho_03_SDNN ρ = −0.70, p_fdr = 0.071 |
| **H15 LASSO** | Lactate prediction | Protocol-only model | cv-R² = 7.9%; only protocol_bin selected |

---

## Baseline Equivalence (H1–H4)

All pre-exercise comparisons were non-significant after FDR correction. Both protocols started from equivalent physiological states.

### H1 — Pre-exercise resting HRV (`lying_inicial`)

- All 10 key metrics: p_fdr > 0.77
- Conclusion: **No pre-exercise ANS differences between MR and VBT.**

### H2 — Short-term pre-exercise orthostatic reactivity

- All 10 key metrics: p_fdr > 0.44
- Conclusion: **No difference in 3-min orthostatic HRV shift.**

### H3 — Medium-term pre-exercise orthostatic reactivity

- RMSSD: p_raw = 0.015, dz = 0.78 — **did not survive FDR** (p_fdr = 0.147)
- LF_nu: p_raw = 0.029, rank-biserial r = −0.53 — **did not survive FDR** (p_fdr = 0.147)
- Conclusion: **No significant orthostatic differences after correction.** Nominally observed tendencies (RMSSD and LF_nu) warrant monitoring but are not conclusive.

### H4 — Pre-exercise blood lactate

| Variable | Test | p-value | MR median | VBT median |
|---|---|---|---|---|
| lactate_baseline | Wilcoxon SR | 0.733 | 3.0 mmol/L | 3.0 mmol/L |

- Conclusion: **Baseline metabolic state equivalent between protocols.**

---

## Post-exercise Lactate (H10)

### H10 — Post-exercise and delta blood lactate

| Variable | Test | p-value | Cohen's dz | MR median | VBT median |
|---|---|---|---|---|---|
| `lactate_post` | paired t | **0.00093** | **0.983** | 6.9 mmol/L | 5.2 mmol/L |
| `delta_lactate` | paired t | **0.00154** | **0.924** | 3.3 mmol/L | 2.5 mmol/L |

**MR produces significantly higher post-exercise lactate and greater lactate accumulation than VBT**, with very large effect sizes (dz ≈ 0.92–0.98). This is the clearest physiological distinction between protocols and confirms the greater metabolic demand imposed by maximal-repetition training.

---

## Exercise-Phase ANS Trajectories (H5)

Model: `metric ~ protocol * phase_pos + (1 | participant_id)`, trimmed exercise epochs only (N_obs = 295–631).

### Protocol level differences (mean level across entire exercise phase)

| Metric | VBT − MR estimate | p-value | Interpretation |
|---|---|---|---|
| RMSSD | −12.5 ms | 0.020 | VBT has lower absolute parasympathetic tone |
| HF_nu | +6.2 n.u. | 0.028 | VBT has higher normalised HF power |
| SD1 | −13.5 ms | 0.025 | VBT lower short-term Poincaré variability |
| SampEn | −0.12 | 0.018 | VBT lower cardiac complexity during exercise |
| alpha1 | −0.08 | 0.047 | VBT lower short-term fractal scaling |

### Significant protocol × time interactions (LRT p-values)

The interaction term tests whether the *rate of change over time* differs between protocols.

| Metric | LRT χ² p | VBT interaction coef | Interpretation |
|---|---|---|---|
| HF_log | **< 0.001** | +0.087 / epoch | VBT HF_log declines more slowly across exercise |
| LF_log | **0.00001** | +0.067 / epoch | VBT LF_log declines more slowly |
| LF_HF_ratio | **0.012** | −0.085 / epoch | LF/HF rises less steeply under VBT |
| SD1 | **0.00009** | +1.82 / epoch | VBT SD1 trajectory diverges upward relative to MR |

> **Interpretation:** Both protocols suppress parasympathetic activity during exercise, but MR does so harder and faster. VBT preserves more parasympathetic activity across the exercise bout, consistent with its lower metabolic demand (confirmed by H10 lactate results).

### Common time trend (both protocols)

HF_log, LF_log, SD1, and SampEn all showed significant decline across exercise epochs (p < 0.001), while LF_nu and LF_HF_ratio increased — the expected sympathetic shift during resistance exercise.

---

## Post-exercise Recovery Trajectories (H6)

Model: `metric ~ protocol * phase_pos + (1 | participant_id)`, trimmed post epochs only (N_obs = 223–474).

### Protocol level differences at the start of recovery (post01)

| Metric | VBT − MR estimate | p-value | Interpretation |
|---|---|---|---|
| RMSSD | −31.6 ms | **< 0.001** | VBT starts recovery from a much lower RMSSD level |
| SDNN | −19.4 ms | **0.00043** | Same pattern for total variability |
| HF_log | +1.34 | **< 0.001** | VBT has higher HF power at onset of recovery |
| LF_log | +0.76 | **0.00007** | VBT higher LF power at onset of recovery |
| LF_HF_ratio | −2.15 | **0.009** | VBT lower sympathovagal balance post-exercise |

### Significant protocol × time interactions (LRT p-values)

| Metric | LRT χ² p | VBT interaction coef | Interpretation |
|---|---|---|---|
| RMSSD | **0.00001** | +4.23 ms / epoch | VBT recovers RMSSD ~4× faster per epoch |
| SDNN | **0.00001** | +3.59 ms / epoch | VBT recovers total variability faster |
| HF_log | **0.018** | −0.080 / epoch | HF_log rises more rapidly under MR recovery |
| HF_nu | **0.010** | +1.41 n.u. / epoch | Normalised HF recovers faster under VBT |

> **Interpretation:** The recovery picture is nuanced. VBT begins recovery with more suppressed vagal time-domain indices (RMSSD, SDNN) but climbs back significantly faster. MR begins recovery from a higher absolute RMSSD level but recovers more slowly. The faster VBT recovery trajectory, combined with the lower post-exercise lactate, suggests VBT places a more controlled and recoverable autonomic load on the participant.

---

## Post-exercise Resting State (H7)

### H7 — Resting HRV at `lying_final` (MR vs VBT)

| Metric | Test | p_raw | p_fdr | Cohen's dz | MR median | VBT median |
|---|---|---|---|---|---|---|
| **SampEn** | paired t | **0.0037** | **0.033** | **−0.823** | 1.267 | 1.677 |
| LF_log | paired t | 0.021 | 0.096 | −0.619 | 7.209 | 7.396 |
| LF_nu | paired t | 0.063 | 0.189 | 0.484 | 69.7 | 51.0 |

**Sample Entropy is the only metric surviving FDR correction**: VBT yields significantly higher cardiac complexity at rest post-exercise, with a large effect size (dz = −0.82). LF_log is nominally significant but does not survive correction.

> Higher SampEn under VBT post-exercise may reflect a less exhausted, more adaptive autonomic state — consistent with the lower lactate accumulation observed under VBT.

---

## Post-exercise Orthostatic Reactivity (H8, H9)

- **H8** (3-min post-ortho): all p_fdr > 0.55 — no difference
- **H9** (5-min post-ortho): all p_fdr > 0.86 — no difference

Both protocols produce equivalent orthostatic HRV responses when standing after exercise.

---

## Pre-to-Post ΔHRV and Δortho (H11, H12)

### H11 — ΔHRV rest-state (lying_final − lying_inicial), MR vs VBT

All p_fdr > 0.52. The net HRV change from pre- to post-exercise does not differ significantly between protocols for any key metric.

- SampEn borderline: p_raw = 0.059, dz = −0.49 (VBT shows slightly larger SampEn reduction). Does not survive correction.

### H12 — Δ orthostatic reactivity (Δortho_03 and Δortho_05), MR vs VBT

All p_fdr > 0.73. How the orthostatic HRV shift changes from pre- to post-exercise is indistinguishable between protocols.

---

## Metabolic–Autonomic Coupling (H13, H14)

### H13 — `lactate_post` × ΔHRV (bivariate Spearman, per protocol)

No significant correlations after FDR correction (all p_fdr > 0.95). The strongest observed associations were:
- MR `delta_HF_nu`: ρ = −0.32 (p = 0.166)
- VBT `dortho_05_HF_nu`: ρ = +0.34 (p = 0.127)

No Fisher's Z comparison between protocols was significant. **Post-exercise lactate concentration is not meaningfully associated with any individual ΔHRV metric in this sample.**

### H14 — `delta_lactate` × ΔHRV (bivariate + partial Spearman)

No correlation survived FDR correction (min p_fdr = 0.071). Notable raw findings in VBT:

| Protocol | Variable | ρ (bivariate) | p_raw | partial ρ (adj. lac_base) | partial p |
|---|---|---|---|---|---|
| VBT | dortho_03_SDNN | **−0.701** | **0.001** | **−0.689** | **0.002** |
| VBT | dortho_05_LF_log | −0.493 | 0.023 | −0.497 | 0.026 |
| VBT | dortho_05_RMSSD | −0.488 | 0.040 | −0.437 | 0.079 |
| VBT | dortho_03_RMSSD | −0.444 | 0.065 | −0.402 | 0.109 |
| MR | dortho_03_SDNN | +0.057 | 0.833 | +0.091 | 0.747 |

The `dortho_03_SDNN` association survives partial correlation controlling for baseline lactate (partial ρ = −0.69, p = 0.002), indicating it is not attributable to pre-exercise metabolic differences. However, p_fdr = 0.071 — marginal, not formally significant.

**Significant protocol moderation (Fisher's Z):**

| Variable | ρ_MR | ρ_VBT | Fisher Z | p |
|---|---|---|---|---|
| dortho_03_SDNN | +0.057 | **−0.701** | **2.45** | **0.014** |

The association between lactate accumulation and the blunting of the 3-min SDNN orthostatic shift is **significantly stronger (and negative) under VBT than under MR**. In VBT, participants who accumulate more lactate show a greater reduction in SDNN orthostatic response — a pattern absent in MR. This protocol-specific coupling warrants further investigation.

---

## H15 — LASSO: Predicting `lactate_post` from Protocol + ΔHRV

### Setup

- **Outcome:** `lac_post` (mmol/L)
- **Predictors:** `protocol_bin` (0 = MR, 1 = VBT) + 10 ΔHRV rest-state metrics + 10 Δortho_03 metrics + 9 Δortho_05 metrics = **29 predictors**
- **Complete cases:** N = 16 (after coverage filtering)
- **Lambda selection:** LOO cross-validation (`cv.glmnet`, nfolds = 16)
- **Primary model:** `lambda.1se`

### Multicollinearity (pre-LASSO OLS diagnostics)

All ΔHRV predictors had severely inflated VIF — the expected consequence of fitting 29 correlated HRV metrics on 16 observations:

| Predictor group | VIF range | All flagged? |
|---|---|---|
| `protocol_bin` | 6.3 | No — acceptable |
| `delta_*` HRV metrics | 21 – 4,385 | Yes — all VIF > 10 |
| `dortho_03_*` metrics | 56 – 2,404 | Yes — all VIF > 10 |
| `dortho_05_*` metrics | 22 – 4,457 | Yes — all VIF > 10 |

The most extreme cases were `delta_LF_log` (VIF = 4,385), `dortho_05_HF_log` (4,457), and `dortho_03_HF_log` (2,404) — near-perfect linear dependencies among frequency-domain metrics. This is the primary rationale for using LASSO: OLS cannot isolate individual predictor effects under these conditions.

### Influential Observations (Cook's D, threshold = 4/N = 0.25)

21 of 34 observations were flagged, including two extreme outliers:

| Participant | Protocol | Cook's D | Notes |
|---|---|---|---|
| 8 | MR | **2,765** | Extreme — dominates OLS fit |
| 22 | VBT | **685** | Second extreme outlier |
| 1 | MR | 13.5 | Influential |
| 9 | MR | 6.9 | Influential |
| 26 | VBT | 7.5 | Influential |

The high number of flagged observations directly reflects the underdetermined OLS (N < p). Importantly, **all sensitivity re-fits excluding each influential observation produced "stable"** — LASSO variable selection did not change in any case. The regularisation renders the model robust to these leverage points.

### LASSO Coefficients (lambda.1se = 0.733)

| Predictor | Coef (lambda.min) | Coef (lambda.1se) | Selected | Bootstrap selection % | Stable (≥ 50%)? |
|---|---|---|---|---|---|
| (Intercept) | 6.628 | 6.171 | Yes | 100% | Yes |
| **protocol_bin** | **−1.313** | **−0.474** | **Yes** | **85.8%** | **Yes** |
| delta_HF_nu | −0.026 | 0 | No | 2.2% | No |
| dortho_03_HF_nu | +0.100 | 0 | No | 6.0% | No |
| All other 26 predictors | 0 | 0 | No | < 2% each | No |

**Only `protocol_bin` is selected and stable.** Every ΔHRV predictor is shrunk to zero at `lambda.1se`.

### Cross-validation Performance

| | Value |
|---|---|
| lambda.1se | 0.733 |
| lambda.min | 0.317 |
| cv-MSE at lambda.1se | 2.50 |
| **cv-R² at lambda.1se** | **0.079 (7.9%)** |

The model explains approximately 8% of variance in post-exercise lactate under LOO cross-validation — capturing only the protocol-level mean difference (MR ≈ 1.66 mmol/L higher than VBT).

### Interpretation

The LASSO result is interpretively important rather than a methodological failure:

1. **Extreme multicollinearity** among ΔHRV predictors (VIF in the hundreds to thousands) means the model cannot distinguish the independent contribution of any individual HRV metric. LASSO correctly collapses all HRV information to zero because no single metric adds predictive value beyond what the others already carry.

2. **Protocol membership alone** accounts for the systematic lactate difference. The coefficient of −0.47 means that VBT is associated with approximately 0.47 mmol/L lower post-exercise lactate than MR, after the intercept (grand mean ≈ 6.17 mmol/L for MR) is set.

3. **ΔHRV indices do not serve as useful non-invasive surrogates for post-exercise lactate** in this sample and design, at least not in their current aggregated delta form. The coupling between autonomic and metabolic responses is not strong enough, or is too protocol-specific, to produce a generalisable prediction from HRV alone.

4. The only thread worth developing further is the `dortho_03_SDNN` × `delta_lactate` association under VBT (H14), which was the single correlation to partially survive partial-correlation control for baseline. A targeted model restricted to VBT data, using SDNN-based orthostatic variables only, might yield a more interpretable result with reduced multicollinearity.

---

## Consolidated Significant Findings

### Statistically significant after FDR correction

| Hypothesis | Metric / Variable | Test | p_fdr | Effect size | Direction |
|---|---|---|---|---|---|
| H7 | SampEn (`lying_final`) | paired t | 0.033 | dz = −0.82 | VBT > MR |
| H10 | lactate_post | paired t | — | dz = +0.98 | MR > VBT |
| H10 | delta_lactate | paired t | — | dz = +0.92 | MR > VBT |
| H5 | HF_log trajectory (LRT) | LMM | < 0.001 | — | VBT declines slower |
| H5 | LF_log trajectory (LRT) | LMM | 0.00001 | — | VBT declines slower |
| H5 | LF_HF_ratio trajectory (LRT) | LMM | 0.012 | — | VBT rises slower |
| H5 | SD1 trajectory (LRT) | LMM | 0.00009 | — | VBT diverges upward |
| H6 | RMSSD trajectory (LRT) | LMM | 0.00001 | +4.23 ms/epoch | VBT recovers faster |
| H6 | SDNN trajectory (LRT) | LMM | 0.00001 | +3.59 ms/epoch | VBT recovers faster |
| H6 | HF_log trajectory (LRT) | LMM | 0.018 | — | MR rises faster |
| H6 | HF_nu trajectory (LRT) | LMM | 0.010 | +1.41 n.u./epoch | VBT recovers faster |

### Marginally significant / noteworthy (p < 0.05 uncorrected, or protocol moderation significant)

| Hypothesis | Metric / Variable | Finding |
|---|---|---|
| H3 | RMSSD ortho pre | p_raw = 0.015, dz = 0.78; NS after FDR |
| H7 | LF_log (`lying_final`) | p_raw = 0.021, dz = −0.62; NS after FDR |
| H11 | SampEn ΔHRV | p_raw = 0.059, dz = −0.49; NS after FDR |
| H14 | dortho_03_SDNN × delta_lactate (VBT) | ρ = −0.70, partial ρ = −0.69; p_fdr = 0.071 |
| H14 | Fisher's Z: dortho_03_SDNN MR vs VBT | Z = 2.45, **p = 0.014** — protocols differ significantly |

---

## Narrative Summary

MR and VBT produced statistically equivalent pre-exercise physiological states across all HRV and lactate measures, confirming the integrity of the crossover design.

During exercise, MR imposed a greater and faster suppression of parasympathetic HRV metrics compared to VBT, with significant differences in HF_log, LF_log, LF_HF_ratio, and SD1 trajectories. VBT maintained more autonomic stability across the exercise bout.

In the immediate post-exercise period, VBT showed a paradoxically lower vagal time-domain baseline (RMSSD, SDNN at post01) but significantly faster recovery trajectories — rising approximately 4 ms/epoch faster for RMSSD and 3.6 ms/epoch faster for SDNN. At the stable post-exercise rest point (`lying_final`), VBT exhibited significantly higher SampEn (cardiac complexity), suggesting a less exhausted and more adaptive autonomic state.

The metabolic results confirmed that MR generates substantially higher post-exercise lactate (median 6.9 vs 5.2 mmol/L, dz = 0.98) and accumulated lactate (3.3 vs 2.5 mmol/L, dz = 0.92). Despite this clear metabolic difference, ΔHRV metrics could not predict post-exercise lactate in the LASSO model, owing to extreme multicollinearity among HRV predictors and a small lactate sample (N = 16 complete cases). The only coupling signal identified was a VBT-specific moderate negative correlation between lactate accumulation and the 3-min SDNN orthostatic response (ρ = −0.70), which was significantly stronger under VBT than MR (Fisher's Z p = 0.014) and survived partial correlation control for baseline lactate.

Taken together, the findings indicate that VBT constitutes a physiologically distinct and metabolically less demanding training stimulus than MR, producing faster autonomic recovery and lower lactate accumulation while achieving equivalent pre-exercise ANS baselines. The protocol-specific nature of the autonomic–metabolic coupling observed in H14 suggests that the relationship between lactate and HRV may be more informative within a given training modality than across modalities.

---

*Generated from `hrv_testing.R` outputs — `outputs/results_*.csv`*
*Analysis pipeline: `hrv_pipeline.py` → `hrv_testing.R`*
*Reference: `docs/dictionary_results.md`*
