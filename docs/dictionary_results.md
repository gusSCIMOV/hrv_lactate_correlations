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
Paired tests comparing **ΔHRV rest-state** (`lying_final − lying_inicial`) between MR and VBT protocols (H11).
`median_MR` / `median_VBT` here represent the median of the *delta* values under each protocol (the input to the test is already the difference score, so the columns retain the generic names produced by `paired_test_table()`).

| Column | Description |
|---|---|
| `hypothesis` | `H11` |
| `metric` | HRV metric name (canonical) |
| `n` | Number of complete MR–VBT pairs |
| `test` | `paired_t` or `wilcoxon_SR` (chosen per Shapiro-Wilk on the paired difference) |
| `statistic` | t-statistic (paired t) or W-statistic (Wilcoxon) |
| `df_approx` | Degrees of freedom (paired t only; `NA` for Wilcoxon) |
| `p_value` | Raw p-value |
| `effect_size` | Cohen's dz (paired t) or rank-biserial r (Wilcoxon) |
| `effect_label` | `Cohen_dz` or `rank_biserial_r` |
| `median_MR` | Median of ΔHRV values (lying_final − lying_inicial) under MR |
| `median_VBT` | Median of ΔHRV values (lying_final − lying_inicial) under VBT |
| `mean_MR` | Mean of ΔHRV values under MR |
| `mean_VBT` | Mean of ΔHRV values under VBT |
| `p_fdr` | BH-adjusted p-value (FDR correction across all key metrics within H11) |
| `sig_raw` | `TRUE` if p_value < 0.05 |
| `sig_fdr` | `TRUE` if p_fdr < 0.05 |

---

#### `results_H12_delta_ortho.csv`
Paired tests comparing **Δ orthostatic reactivity** (`ortho_post − ortho_pre`, both the 3-min and 5-min windows) between MR and VBT protocols (H12).
Two hypothesis blocks are stacked in one file: `H12_ortho_03` (3-min window) and `H12_ortho_05` (5-min window).

> Values are second-order deltas: the difference of % orthostatic shifts from pre- to post-exercise. A positive value means standing caused a larger HRV change post-exercise than pre-exercise.

| Column | Description |
|---|---|
| `hypothesis` | `H12_ortho_03` or `H12_ortho_05` |
| `metric` | HRV metric name (canonical) |
| `n` | Number of complete MR–VBT pairs |
| `test` | `paired_t` or `wilcoxon_SR` |
| `statistic` | t-statistic or W-statistic |
| `df_approx` | Degrees of freedom (paired t only; `NA` for Wilcoxon) |
| `p_value` | Raw p-value |
| `effect_size` | Cohen's dz or rank-biserial r |
| `effect_label` | `Cohen_dz` or `rank_biserial_r` |
| `median_MR` | Median of Δortho values under MR |
| `median_VBT` | Median of Δortho values under VBT |
| `mean_MR` | Mean of Δortho values under MR |
| `mean_VBT` | Mean of Δortho values under VBT |
| `p_fdr` | BH-adjusted p-value across all key metrics within each hypothesis block |
| `sig_raw` | `TRUE` if p_value < 0.05 |
| `sig_fdr` | `TRUE` if p_fdr < 0.05 |

---

#### `results_H13_lactate_post_spearman.csv`
Bivariate Spearman correlations between `lac_post` and each ΔHRV / Δortho index, run separately for MR and VBT, with protocol comparison via Fisher's Z (H13).

> `protocol` takes values `"MR"` or `"VBT"` only — no pooled row is produced. Inter-protocol differences are captured by `fisher_Z_MR_vs_VBT` / `fisher_p_MR_vs_VBT`, which repeat the same value on both protocol rows for a given `hrv_delta_var`.

| Column | Description |
|---|---|
| `hypothesis` | `H13` |
| `protocol` | `MR` or `VBT` |
| `lactate_var` | `lac_post` (post-exercise blood lactate as stored in the coupling table) |
| `hrv_delta_var` | ΔHRV predictor name (e.g. `delta_RMSSD`, `dortho_03_RMSSD`, `dortho_05_HF_log`) |
| `n` | Number of complete pairs for this protocol × variable combination |
| `rho` | Spearman ρ |
| `rho_CI_low` | Percentile bootstrap 2.5th percentile (2 000 iterations) |
| `rho_CI_high` | Percentile bootstrap 97.5th percentile |
| `p_value` | Raw p-value (normal approximation on the Spearman t-statistic) |
| `p_fdr` | BH-adjusted p-value across all protocol × variable combinations within H13 |
| `sig_fdr` | `TRUE` if p_fdr < 0.05 |
| `fisher_Z_MR_vs_VBT` | Fisher's Z statistic testing ρ_MR ≠ ρ_VBT for this variable |
| `fisher_p_MR_vs_VBT` | Two-tailed p-value for the Fisher's Z test |

---

#### `results_H14_delta_lactate_spearman.csv`
Bivariate and partial Spearman correlations between `lac_delta` and each ΔHRV / Δortho index, with and without controlling for `lac_base` (H14).
Output is a merged table: bivariate results from `run_spearman()` joined with partial-correlation results from `ppcor::pcor()`.

> `protocol` takes values `"MR"` or `"VBT"` only. `lactate_var` in the output is `"lac_delta"` (the internal coupling-table column name for `lactate_post − lactate_baseline`).

| Column | Description |
|---|---|
| `protocol` | `MR` or `VBT` |
| `hrv_delta_var` | ΔHRV predictor name (e.g. `delta_RMSSD`, `dortho_03_RMSSD`) |
| `hypothesis` | `H14` |
| `lactate_var` | `lac_delta` (= `lac_post − lac_base` within each protocol) |
| `n` | Number of complete pairs for the bivariate correlation |
| `rho` | Bivariate Spearman ρ |
| `rho_CI_low` | Percentile bootstrap 2.5th percentile (2 000 iterations) |
| `rho_CI_high` | Percentile bootstrap 97.5th percentile |
| `p_value` | Raw bivariate p-value |
| `p_fdr` | BH-adjusted p-value (bivariate) across all combinations within H14 |
| `sig_fdr` | `TRUE` if p_fdr < 0.05 |
| `fisher_Z_MR_vs_VBT` | Fisher's Z statistic comparing ρ_MR vs ρ_VBT |
| `fisher_p_MR_vs_VBT` | Two-tailed p-value for the Fisher's Z test |
| `partial_rho` | Partial Spearman ρ of `lac_delta` ~ `hrv_delta_var` controlling for `lac_base` (via `ppcor::pcor`); `NA` if fewer than 5 complete triples |
| `partial_p` | p-value for the partial correlation |
| `n_partial` | Number of complete triples (`lac_delta`, `hrv_delta_var`, `lac_base`) used for the partial correlation |

---

#### H15 output files (4 separate CSVs)

There is no single `results_H15_lasso_lactate_post.csv`. H15 produces four independent files described below.

---

#### `results_H15_lasso_coefficients.csv`
Final LASSO model coefficients at both `lambda.min` and `lambda.1se`, plus bootstrap selection stability.
The first row is always `(Intercept)`; remaining rows are the predictors (z-scored continuous variables and `protocol_bin`).

| Column | Description |
|---|---|
| `predictor` | Predictor name; first row is `(Intercept)` |
| `lambda_min` | Coefficient at `lambda.min` (lowest CV error) |
| `lambda_1se` | Coefficient at `lambda.1se` (primary reported model; more parsimonious) |
| `selected` | `TRUE` if coefficient is non-zero at `lambda.1se` |
| `boot_selection_pct` | % of 500 bootstrap resamples in which this predictor had a non-zero coefficient at the same `lambda.1se` |
| `stable` | `TRUE` if `boot_selection_pct` ≥ 50 % |

---

#### `results_H15_lasso_diagnostics.csv`
Pre-LASSO OLS Variance Inflation Factors for multicollinearity screening.
**Produced only when N > p + 2** (OLS is estimable); absent otherwise.

| Column | Description |
|---|---|
| `predictor` | Predictor name (matches predictors in `results_H15_lasso_coefficients.csv`, excluding intercept) |
| `VIF` | Variance Inflation Factor from the saturated OLS fit |
| `VIF_flag` | `TRUE` if VIF > 10 (severe multicollinearity; consider removing or combining) |

---

#### `results_H15_lasso_cv.csv`
Full regularisation path cross-validation performance from `glmnet::cv.glmnet()` (LOO-CV, `nfolds = N`).

| Column | Description |
|---|---|
| `lambda` | Lambda value on the regularisation path |
| `cvm` | Mean cross-validated MSE at this lambda |
| `cvsd` | Standard deviation of cross-validated MSE |
| `lambda_is_min` | `TRUE` for the lambda that minimises mean CV error (`lambda.min`) |
| `lambda_is_1se` | `TRUE` for the 1-SE rule lambda (`lambda.1se`; primary model) |
| `cv_R2_lambda_1se` | Cross-validated R² at `lambda.1se`: `1 − (MSE / Var(lac_post))`; populated only on the `lambda_is_1se = TRUE` row, `NA` elsewhere |

---

#### `results_H15_influential_obs.csv`
Cook's distance diagnostic from the pre-LASSO OLS fit, flagging potentially influential observations.
**Produced only when OLS is estimable (N > p + 2) and at least one observation exceeds the Cook's D threshold.**

| Column | Description |
|---|---|
| `participant_id` | Participant ID |
| `protocol` | `MR` or `VBT` |
| `cooks_d` | Cook's distance from the OLS fit |
| `influential` | `TRUE` if Cook's D > 4/N |
| `sensitivity_excluded` | `"stable"` if LASSO variable selection is unchanged when this observation is excluded; `"selection_changed"` if the selected predictor set differs; `NA` if the sensitivity refit failed |
