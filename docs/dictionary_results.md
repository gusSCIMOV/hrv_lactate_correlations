### Output Files and Column Dictionaries


[Back to README](../README.md)
---

#### `results_00_shapiro_normality.csv`
Normality of the within-subject paired difference (MR − VBT) per metric and epoch slot.
Used internally by the script to select the appropriate test for each metric.

| Column | Description |
|---|---|
| `hypothesis` | Hypothesis label the check applies to (e.g. `H1_lying_inicial`) |
| `metric` | HRV metric name (canonical) |
| `n_pairs` | Number of complete MR–VBT pairs available |
| `mean_diff` | Mean of (MR − VBT) differences |
| `SD_diff` | Standard deviation of (MR − VBT) differences |
| `SW_W` | Shapiro-Wilk W statistic |
| `SW_p` | Shapiro-Wilk p-value |
| `normal` | `TRUE` if SW_p > 0.05 (difference distribution consistent with normality) |

---

#### `results_00_levene_variance.csv`
Variance homogeneity check between MR and VBT groups at fixed supine epochs (H1 and H7).

| Column | Description |
|---|---|
| `hypothesis` | Hypothesis label (`H1_lying_inicial` or `H7_lying_final`) |
| `metric` | HRV metric name |
| `Levene_F` | Levene F-statistic |
| `Levene_p` | p-value of Levene's test |
| `equal_var` | `TRUE` if Levene_p > 0.05 (variances not significantly different) |

---

#### `results_H1_H7_fixed_epochs_hrv.csv`
Paired tests comparing MR vs VBT at fixed supine rest epochs.
H1 = pre-exercise resting ANS (`lying_inicial`); H7 = post-exercise resting ANS (`lying_final`).

| Column | Description |
|---|---|
| `hypothesis` | `H1` or `H7` |
| `metric` | HRV metric name |
| `n` | Number of complete MR–VBT pairs |
| `test` | `paired_t` or `wilcoxon_SR` (chosen per SW result) |
| `statistic` | t-statistic (paired t) or W-statistic (Wilcoxon) |
| `df_approx` | Degrees of freedom (paired t only; NA for Wilcoxon) |
| `p_value` | Raw p-value |
| `effect_size` | Cohen's dz (paired t) or rank-biserial r (Wilcoxon) |
| `effect_label` | `Cohen_dz` or `rank_biserial_r` |
| `median_MR` | Median of MR values at this epoch |
| `median_VBT` | Median of VBT values at this epoch |
| `mean_MR` | Mean of MR values at this epoch |
| `mean_VBT` | Mean of VBT values at this epoch |
| `p_fdr` | BH-adjusted p-value (FDR correction across all metrics within the hypothesis) |
| `sig_raw` | `TRUE` if p_value < 0.05 |
| `sig_fdr` | `TRUE` if p_fdr < 0.05 |

---

#### `results_H2_H3_H8_H9_ortho_pct.csv`
Paired tests comparing MR vs VBT orthostatic reactivity, expressed as **% change** relative to the lying baseline.

- H2 = short-term pre-exercise ortho (`stand_up_inicial_03` vs `lying_inicial`)
- H3 = medium-term pre-exercise ortho (`stand_up_inicial_05` vs `lying_inicial`)
- H8 = short-term post-exercise ortho (`stand_up_final_03` vs `lying_final`)
- H9 = medium-term post-exercise ortho (`stand_up_final_05` vs `lying_final`)

Columns are identical to `results_H1_H7_fixed_epochs_hrv.csv` with `hypothesis` values `H2`, `H3`, `H8`, `H9`.

> Note: values here are **% changes**, not raw HRV units. A positive value means the metric increased on standing relative to lying; a negative value means it decreased. The test compares whether that % shift differs between MR and VBT.

---

#### `results_H4_H10_lactate.csv`
Paired tests on blood lactate variables. H4 = pre-exercise baseline; H10 = post-exercise and delta.

| Column | Description |
|---|---|
| `hypothesis` | `H4` or `H10` |
| `variable` | `lactate_baseline`, `lactate_post`, or `delta_lactate` |
| `n` | Number of complete MR–VBT pairs |
| `SW_p_diff` | Shapiro-Wilk p-value on the paired difference (reported for transparency) |
| `test` | `paired_t` or `wilcoxon_SR` |
| `statistic` | t or W statistic |
| `df_approx` | Degrees of freedom (paired t only; NA for Wilcoxon) |
| `p_value` | Raw p-value (no FDR applied — single variable per hypothesis) |
| `effect_size` | Cohen's dz or rank-biserial r |
| `effect_label` | `Cohen_dz` or `rank_biserial_r` |
| `median_MR` | Median lactate value (MR) |
| `median_VBT` | Median lactate value (VBT) |
| `mean_MR` | Mean lactate value (MR) |
| `mean_VBT` | Mean lactate value (VBT) |

> `delta_lactate` = `lactate_post − lactate_baseline` within each protocol.

---

#### `results_H5_exercise_LMM.csv`
Linear Mixed Model results for the **exercise-phase** ANS trajectory (H5).
Model: `metric ~ protocol * phase_pos + (1 | participant_id)`, fitted on `trimmed == TRUE` exercise epochs only.
One row per model coefficient per metric.

| Column | Description |
|---|---|
| `hypothesis` | `H5_exercise_trajectory` |
| `metric` | HRV metric name |
| `term` | Model term: `(Intercept)`, `protocolVBT`, `phase_pos`, `protocolVBT:phase_pos` |
| `estimate` | Coefficient estimate (unstandardised) |
| `SE` | Standard error of the estimate |
| `df_satterthwaite` | Satterthwaite degrees of freedom (lmerTest) |
| `t_value` | t-statistic |
| `p_value` | p-value based on Satterthwaite approximation |
| `LRT_interaction_chisq` | Chi-square statistic from Likelihood Ratio Test comparing full vs no-interaction model |
| `LRT_interaction_p` | LRT p-value for the `protocol × phase_pos` interaction (same value repeated across all rows for a given metric) |
| `singular` | `TRUE` if the random-effects structure is singular (interpret with caution) |
| `n_obs` | Total number of observations used |
| `n_subjects` | Number of unique participants |

**Key terms to focus on:**

| Term | Interpretation |
|---|---|
| `protocolVBT` | Overall level difference: is VBT higher/lower than MR across the exercise phase? |
| `phase_pos` | Common time trend: does the metric change over exercise epochs regardless of protocol? |
| `protocolVBT:phase_pos` | **Trajectory difference**: does the rate of change over time differ between protocols? |

---

#### `results_H6_post_LMM.csv`
Same structure as `results_H5_exercise_LMM.csv` but for the **post-exercise recovery trajectory** (H6).
Fitted on `trimmed == TRUE` rows where `epoch` matches `^post` (i.e., `post01`–`postN`).

Column definitions are identical to H5 above; `hypothesis` value is `H6_post_trajectory`.

> `lying_final` and `stand_up_final_*` are **excluded** from this model — they are tested separately as fixed-point comparisons in H7–H9.

---

#### `results_H11_delta_hrv.csv`
Paired tests comparing **ΔHRV** (lying_final − lying_inicial) between MR and VBT protocols (H11).

| Column | Description |
|---|---|
| `hypothesis` | `H11` |
| `metric` | HRV metric name |
| `n` | Number of complete MR–VBT pairs |
| `test` | `paired_t` or `wilcoxon_SR` |
| `statistic` | t or W statistic |
| `df_approx` | Degrees of freedom (paired t only; NA for Wilcoxon) |
| `p_value` | Raw p-value |
| `effect_size` | Cohen's dz or rank-biserial r |
| `effect_label` | `Cohen_dz` or `rank_biserial_r` |
| `median_delta_MR` | Median ΔHRV under MR |
| `median_delta_VBT` | Median ΔHRV under VBT |
| `p_fdr` | BH-adjusted p-value across all metrics |
| `sig_fdr` | `TRUE` if p_fdr < 0.05 |

---

#### `results_H12_delta_ortho.csv`
Paired tests comparing **delta orthostatic reactivity** (ortho_post − ortho_pre, for _03 and _05 windows) between MR and VBT protocols (H12).

Columns identical to `results_H11_delta_hrv.csv`; `hypothesis` values are `H12_ortho_03` and `H12_ortho_05`.

> Values here are differences of % changes (i.e., second-order deltas): how much the orthostatic HRV shift changes from pre- to post-exercise, and whether that change differs by protocol.

---

#### `results_H13_lactate_post_spearman.csv`
Bivariate Spearman correlations between `lactate_post` and each ΔHRV index, run per protocol and tested for protocol difference via Fisher's Z (H13).

| Column | Description |
|---|---|
| `hypothesis` | `H13` |
| `protocol` | `MR`, `VBT`, or `pooled` |
| `lactate_var` | `lactate_post` |
| `hrv_delta_var` | ΔHRV variable name (e.g. `delta_RMSSD`, `Δortho_03_RMSSD`) |
| `n` | Number of pairs |
| `rho` | Spearman ρ |
| `rho_CI_low` | Bootstrap 2.5th percentile (2 000 iterations) |
| `rho_CI_high` | Bootstrap 97.5th percentile |
| `p_value` | Raw p-value |
| `p_fdr` | BH-adjusted p-value across all metrics within H13 |
| `sig_fdr` | `TRUE` if p_fdr < 0.05 |
| `fisher_Z_MR_vs_VBT` | Fisher's Z statistic comparing ρ_MR vs ρ_VBT |
| `fisher_p_MR_vs_VBT` | p-value for protocol difference in correlation |

---

#### `results_H14_delta_lactate_spearman.csv`
Bivariate and partial Spearman correlations between `delta_lactate` and each ΔHRV index, with and without controlling for `lactate_baseline` (H14).

| Column | Description |
|---|---|
| `hypothesis` | `H14` |
| `protocol` | `MR`, `VBT`, or `pooled` |
| `lactate_var` | `delta_lactate` |
| `hrv_delta_var` | ΔHRV variable name |
| `n` | Number of pairs |
| `rho` | Spearman ρ (bivariate) |
| `rho_CI_low` | Bootstrap 2.5th percentile |
| `rho_CI_high` | Bootstrap 97.5th percentile |
| `p_value` | Raw bivariate p-value |
| `partial_rho` | Partial Spearman ρ controlling for `lactate_baseline` (via `ppcor::pcor`) |
| `partial_p` | p-value for partial correlation |
| `p_fdr` | BH-adjusted p-value (bivariate) across all metrics within H14 |
| `sig_fdr` | `TRUE` if p_fdr < 0.05 |
| `fisher_Z_MR_vs_VBT` | Fisher's Z statistic |
| `fisher_p_MR_vs_VBT` | p-value for protocol difference |

---

#### `results_H15_lasso_lactate_post.csv`
LASSO regression results predicting `lactate_post` from `protocol` and all ΔHRV predictors (H15).

**Sub-file A — `results_H15_lasso_coefficients.csv`:** Final model coefficients.

| Column | Description |
|---|---|
| `predictor` | Predictor variable name (z-scored continuous; `protocol` as binary) |
| `lambda_min` | Coefficient at `lambda.min` |
| `lambda_1se` | Coefficient at `lambda.1se` (primary reported model) |
| `selected` | `TRUE` if coefficient is non-zero at `lambda.1se` |
| `boot_selection_pct` | % of 500 bootstrap samples in which this predictor was selected (non-zero) |
| `stable` | `TRUE` if `boot_selection_pct` ≥ 50 % |

**Sub-file B — `results_H15_lasso_diagnostics.csv`:** Pre-LASSO OLS diagnostic checks.

| Column | Description |
|---|---|
| `predictor` | Predictor name |
| `VIF` | Variance Inflation Factor from OLS fit |
| `VIF_flag` | `TRUE` if VIF > 10 |

**Sub-file C — `results_H15_lasso_cv.csv`:** Cross-validation performance.

| Column | Description |
|---|---|
| `lambda` | Lambda value tested |
| `cvm` | Mean cross-validated MSE |
| `cvsd` | SD of cross-validated MSE |
| `lambda_min` | `TRUE` for the lambda minimising MSE |
| `lambda_1se` | `TRUE` for the 1-SE rule lambda |
| `cv_R2_lambda_1se` | Cross-validated R² at `lambda.1se`: `1 − (MSE / Var(lactate_post))` |

**Sub-file D — `results_H15_influential_obs.csv`:** Cook's D flagged observations.

| Column | Description |
|---|---|
| `participant_id` | Participant ID |
| `protocol` | `MR` or `VBT` |
| `cooks_d` | Cook's distance from OLS fit |
| `influential` | `TRUE` if Cook's D > 4/N |
| `sensitivity_excluded` | Reports whether LASSO result changes substantially when this observation is left out |
