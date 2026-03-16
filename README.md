# HRV Analysis Project — MR (maximal resistance) vs VBT (velocity-based Training) Strength Protocols in Women
## Project Documentation & Analysis Plan

---

## 1. HRV ANALYSIS

---

## 2. Project Goals

Compare the autonomic nervous system (ANS) response, excercise, recovery phases, and lactate clearence between two resistance training protocols in women:

- **MR** — Maximal Repetitions (traditional failure-based loading)
- **VBT** — Velocity-Based Training (load controlled by movement velocity)

Specifically:
1. Characterize **ANS temporal dynamics** (pre → exercise → post) for each protocol via HRV metrics
2. Identify **between-protocol differences** in HRV and lactate at matched time points
3. Assess **autonomic recovery** completeness and speed post-exercise
4. Quantify **metabolic-autonomic coupling** via lactate ↔ HRV correlations

---

## 3. Study Design

| Feature | Detail |
|---|---|
| Design | **Within-subjects crossover** — every participant completed both MR and VBT |
| Total participants | 27 (lactate data); 21 with MR HRV files; 24 with VBT HRV files |
| Sex | Women |
| Unit of analysis | Participant × Protocol (paired observations) |

> MR and VBT observations per participant are not independent. All inferential tests are paired/within-subject approaches.

---

## 4. Data

### 4.1 HRV Files
- **Location:** `RM_xlsx/` → `{ID}_{Name}_{Surname}_RM.xlsx`; `VPM_xlsx/` → `{ID}_{Name}_{Surname}_VPM.xlsx`
- **Structure:** One file per participant × protocol. Each row = one 5-minute epoch window. Columns = HRV metrics exported from HRV analysis software.
- **Epochs:** Up to 37 named windows spanning baseline → exercise → recovery → final rest + orthostatic tests

### 4.2 Epoch Structure

| Phase | Epochs | Notes |
|---|---|---|
| Baseline | `lying_inicial`, `stand_up_inicial_03`, `stand_up_inicial_05` | Fixed; pre-exercise resting + orthostatic |
| Exercise | `exer01` – `exerN` | Variable length per subject (protocol-dependent) |
| Recovery | `post01` – `postN` | Variable length per subject |
| Final rest | `lying_final`, `stand_up_final_03`, `stand_up_final_05` | Fixed; post-exercise resting + orthostatic |

### 4.3 HRV Metrics Extracted

- **As specified within ./config as KEY_METRICS or key HRV metrics (see hypothesis box below)**

### 4.4. Orthostatic reactivity (Computed within R)

- pre-exercise orthostatic reactivity (stand_up_inicial_03 and _05 vs lying_inicial)
- post-exercise orthostatic reactivity (stand_up_final_03 and _05 vs lying_final)

- **Orthostatic reactivity**: computed as % change = (stand_up − lying) / |lying| × 100 per participant × protocol, then the MR vs VBT difference is tested — not the raw stand_up values.

### 4.5 Blood Lactate
- **File:** `Lactate.xlsx` — 27 participants, 4 columns: `RM_baseline_Lactate`, `RM_post_Lactate`, `VPM_baseline_Lactate`, `VPM_post_Lactate`
- **Mapped to:** `lying_inicial` (baseline) and `lying_final` (post-exercise) epochs only
- **Joined on:** numeric participant ID 

---

## 5. Hypotheses to Test

### 5.1 Phases Differences (MR vs VBT)

For each HRV metric and lactate variable, at each epoch slot:

| # | Hypothesis | Variables | Epoch(s) |
|---|---|---|---|
| H1 | Pre-exercise ANS state does not differ between protocols pre-exercise | key HRV metrics (apply FDR) | `lying_inicial` |
| H2 | Short-term Pre-exercise Orthostatic reactivity does not differ between protocols | key HRV metrics (appy FDR) |  `stand_up_inicial_03`, `lying_inicial` |
| H3 | Medium-term Pre-exercise Orthostatic reactivity does not differ between protocols | key HRV metrics (appy FDR) |  `stand_up_inicial_05`, `lying_inicial` |
| H4 | Pre-exercise blood lactate does not differ between protocols | `lactate_baseline` | `lying_inicial` |
| H6 | Exercise-phase ANS trajectory differs between protocols | key HRV metrics (Linear Mixed Model) | `exer n `–`exerN` (trimmed) |
| H6 | Post-exercise ANS trajectory differs between protocols | key HRV metrics (Linear Mixed Model) | `post n `–`postN` (trimmed) |
| H7 | Post-exercise resting ANS state differs between protocols | key HRV metrics (apply FDR) | `lying_final` |
| H8 | Short-term Post-exercise Orthostatic reactivity does not differ between protocols | key HRV metrics (appy FDR) |  `stand_up_final_03`, `lying_final` |
| H9 | Medium-term Post-exercise Orthostatic reactivity does not differ between protocols | key HRV metrics (appy FDR) |  `stand_up_final_05`, `lying_final` |
| H10 | Post-exercise  and delta blood lactate (Δlactate) does not differ between protocols | `lactate_post`, `delta_lactate` | `lying_final` |

### 5.2 Pre–Post Changes between Protocol (delta measurements)

| # | Hypothesis | Variables | Epoch comparison |
|---|---|---|---|
| H11 | HRV changes from pre to post exercise (delta HRV measurements or ΔHRV) doesnt´t difer between protocol | All HRV metrics | `lying_final`  to `lying_inicial`|
| H12 | Delta Orthostatic reactivity (from pre to post exercise) doesnt´t difer between protocol | All HRV metrics | (`ortho_pre_03, ortho_post_03`) and (`ortho_pre_05, ortho_post_05`)|

### 5.4 Metabolic–Autonomic Coupling (including delta (Δ) measurements)

| # | Hypothesis | Variables |
|---|---|---|
| H13 | Post-exercise lactate correlates with ΔHRV indices  between rest states (`lying_final − lying_inicial`) and Orthostatic reactivity short (`ortho_pre_03, ortho_post_03`) and medium (`ortho_pre_05, ortho_post_05`) deltas | `lactate_post` × ΔHRV |
| H14 | Δlactate correlates with ΔHRV indices between rest states (`lying_final − lying_inicial`) and Orthostatic reactivity short (`ortho_pre_03, ortho_post_03`) and medium (`ortho_pre_05, ortho_post_05`) deltas| `delta_lactate` × all ΔHRV |
| H15 | Association between lactate and HRV and : Multivariate fatigue prediction (lactate) | LASSO regression |

---

## 6. Python Pipeline — Current Outline (`hrv_pipeline.py`)

### Overview

```
project_root/
├── config.py          # ← edit PROJECT_ROOT here before running
├── hrv_utils.py       # helper functions (imported by pipeline)
├── hrv_pipeline.py    # main script — run this
├── RM_xlsx/           # *_RM.xlsx files
├── VPM_xlsx/          # *_VPM.xlsx files
├── Lactate.xlsx
└── outputs/           # all CSVs and figures written here
```

| Function | Purpose |
|---|---|
| `hrv_pipeline.py` | main script: performs EDA, plots, exports `hrv_long_master.csv` — full long-format DataFrame, `hrv_wide_delta.csv` — one row per participant × protocol with phase-mean HRV + lactate columns |
| `config.py` | plain Python module-level variables used  in the hrv_pipeline icluding input paths/files. Set `KEY_METRICS` to include specific HRV metrics  |
| `hrv_utils.py` | Exploratory data analysis  and plot helpers |


###  EDA & Visualization (hrv_utils.py)

| Output file | Content |
|---|---|
| `01_missingness_MR.png` | Heatmap of missing data proportion by epoch × metric |
| `02_dist_qq_{slot}.png` | Histograms + Q-Q plots for 6 fixed epochs + pooled exercise/post (trimmed and untrimmed) |
| `03_violin_{slot}.png` | Ggplot-style violins + strip for all HRV metrics across 10 slots (6 fixed + 4 pooled) |
| `04_temporal_trajectories_untrimmed.png` | Median ± IQR trajectory per metric, all epochs, unbalanced N |
| `04_temporal_trajectories_trimmed.png` | Same with trimmed epochs (balanced N per phase position) |
| `04b_boxplot_traj_untrimmed.png` | Narrow side-by-side boxplots per epoch + dotted median trend line, untrimmed |
| `04b_boxplot_traj_trimmed.png` | Same, trimmed |
| `05_pairplot_{slot}.png` | Pairplots for 6 fixed epochs + pooled exercise/post |
| `06_corr_{epoch}_{protocol}.png` | Spearman correlation matrix, upper triangle + dendrogram, per epoch × protocol |
| `06_corr_all_epochs_{protocol}.png` | Same, participant-median aggregated across all epochs |
| `07_lactate_{pair}_{x_col}.png` | Lactate vs ΔHRV scatterplots with regression lines + confidence ellipses per protocol |
| `09_violin_lactate.png` | Lactate pre / post / delta violin comparison MR vs VBT |

---

## Usage :

- **Python 3.10+** recommended (tested on 3.10)
- Install dependencies:
```bash
pip install pandas numpy matplotlib seaborn scipy openpyxl
```

Open `config.py` and set your local project root:
```python
# config.py  — line 12
PROJECT_ROOT = Path("C:/your/path/to/project_root")
```

### Run
```bash
cd /path/to/project_root
python hrv_pipeline.py
```
---

## 7. R - analysis (version-arch-installed-packages)

```bash
R.version
               _                                
platform       x86_64-w64-mingw32               
arch           x86_64                           
os             mingw32                          
crt            ucrt                             
system         x86_64, mingw32                  
status                                          
major          4                                
minor          5.1                              
year           2025                             
month          06                               
day            13                               
svn rev        88306                            
language       R                                
version.string R version 4.5.1 (2025-06-13 ucrt)
nickname       Great Square Root  

```
### 7.1 Input Files

**`hrv_long_master.csv`** — primary file for mixed models and time-series tests. Take `epoch` as unique epoch label

| Column | Description |
|---|---|
| `participant_id` | Numeric ID (1–27) |
| `participant_name` | Full name |
| `protocol` | `MR` or `VBT` |
| `epoch` | Canonical UNIQUE epoch label (`lying_inicial`, `exer01`, etc.) |
| `epoch_index` | Global integer position (1–37) |
| `phase` | descriptive labels for  `baseline`, `exercise`, `recovery` |
| `RMSSD`, `SDNN`, ... | All HRV metrics (canonical names) |
| `lactate` | Blood lactate (mmol/L); populated only at `lying_inicial` and `lying_final` |


### 7.2 Assumption Checks (per metric × epoch/phase)

```r
# Normality
shapiro.test(x)                          # per group, n < 50
ggplot + stat_qq()                       # visual

# Variance homogeneity
car::leveneTest(metric ~ protocol)

# Sphericity (repeated measures ANOVA only)
rstatix::anova_test()                    # includes Mauchly + GG correction
```

### 7.3 Hypothesis Testing Decision Tree

```
Single epoch and pre-post (delta) comparison between MR vs VBT:
  Normal + homogeneous variance  → paired t-test
  Non-normal or small N          → Wilcoxon signed-rank (paired)
  Apply FDR correction           → (`p.adjust(method="BH")`) across all the key HRV metrics for each epoch
  Effect size                    → Cohen's dz (effsize) / rank-biserial r (rstatix)

Multiple epochs within/between protocols:
  Primary recommendation         → Linear Mixed Model
    lme4: metric ~ protocol * epoch_index + (protocol | participant_id)
  Non-parametric alternative     → Friedman + Scheirer-Ray-Hare
  Effect size                    → partial η² (effectsize) / Kendall's W

```

### 7.4 Correlation & Regression

Analyses in this section correspond directly to **H13–H15** (§5.4 Metabolic–Autonomic Coupling). All correlations are run **per protocol** (MR and VBT separately), with a Fisher's Z test to evaluate whether the association differs between protocols. 

---

#### H13 — Post-exercise lactate × ΔHRV bivariate associations

`lactate_post` is the **dependent variable of interest** (metabolic stress outcome). ΔHRV indices and orthostatic reactivity deltas are the non-invasive surrogates tested as correlates.

| Goal | Method | Assumption checks |
|---|---|---|
| `lactate_post` × ΔHRV rest-state (`lying_final − lying_inicial`) per key metric | Spearman ρ + 95% CI via bootstrap (`boot`, 2 000 iterations) | Scatter plot linearity, outlier inspection, Cook's D |
| `lactate_post` × short-term orthostatic reactivity delta (`Δortho_03 = ortho_post_03 − ortho_pre_03`) | Spearman ρ + 95% CI via bootstrap | As above |
| `lactate_post` × medium-term orthostatic reactivity delta (`Δortho_05 = ortho_post_05 − ortho_pre_05`) | Spearman ρ + 95% CI via bootstrap | As above |
| Protocol comparison (ρ_MR vs ρ_VBT) | Fisher's Z transformation | Adequate N per group (≥ 10) |
| FDR correction | `p.adjust(method = "BH")` across all metrics within H13 | — |

---

#### H14 — Δlactate × ΔHRV bivariate and partial associations

`delta_lactate` (`lactate_post − lactate_baseline`) captures within-session metabolic load change. Correlations mirror H13 but use the change score as the lactate index.

| Goal | Method | Assumption checks |
|---|---|---|
| `delta_lactate` × ΔHRV rest-state (`lying_final − lying_inicial`) per key metric | Spearman ρ + 95% CI via bootstrap | Scatter linearity, outliers, Cook's D |
| `delta_lactate` × `Δortho_03` and `Δortho_05` | Spearman ρ + 95% CI via bootstrap | As above |
| Controlling for `lactate_baseline` or baseline HRV | Partial Spearman correlation (`ppcor::pcor`) | Residual normality (Shapiro-Wilk on partial residuals) |
| Protocol comparison (ρ_MR vs ρ_VBT) | Fisher's Z transformation | Adequate N per group |
| FDR correction | `p.adjust(method = "BH")` across all metrics within H14 | — |

---

#### H15 — Multivariate prediction of post-exercise lactate (LASSO)

**Framing:** `lactate_post` is the **outcome**. The goal is to identify which combination of non-invasive HRV-based delta measurements (and protocol type) best predicts post-exercise blood lactate concentration — establishing ΔHRV indices as metabolic surrogates.

| Element | Detail |
|---|---|
| **Outcome** | `lactate_post` (continuous, mmol/L) |
| **Predictors** | `protocol` (binary: 0 = MR, 1 = VBT); ΔHRV rest-state per key metric (`lying_final − lying_inicial`); `Δortho_03` and `Δortho_05` per key metric |
| **Reason for LASSO** | N ≈ 27 (lactate sample); predictor count easily exceeds N/10 threshold with 10 key metrics × delta types → regularisation required. LASSO performs automatic variable selection by shrinking non-informative coefficients to zero |
| **Lambda selection** | `glmnet::cv.glmnet`, leave-one-out cross-validation (`nfolds = nrow(df)`) given small N; `lambda.1se` preferred for parsimony |
| **Predictor scaling** | All continuous predictors z-scored before fitting (mean = 0, SD = 1) so LASSO penalties are comparable |
| **LASSO stability** | Bootstrap resampling (500 iterations) to report selection frequency per predictor; only predictors selected in ≥ 50 % of bootstrap samples are considered stable |

**Assumption & diagnostic checks (run on OLS fit prior to LASSO):**

| Check | Tool | Threshold / action |
|---|---|---|
| Multicollinearity | VIF via `car::vif()` on OLS | VIF > 10 → flag predictor; consider removing or combining |
| Influential observations | Cook's D via `stats::cooks.distance()` | Cook's D > 4/N → inspect and report; re-run LASSO excluding the observation as sensitivity check |
| Residual normality | Shapiro-Wilk on OLS residuals | If violated, consider rank-based (Spearman-based) variable selection or report as caveat |
| Heteroscedasticity | Breusch-Pagan test (`lmtest::bptest`) | Significant → use robust standard errors on the final OLS interpretation |

> **Note:** Because all predictors are ΔHRV measurements (non-invasive, field-accessible), a parsimonious LASSO model with good cross-validated R² would support the clinical claim that HRV monitoring can substitute for blood draws in estimating post-exercise lactate.


### 7.5 R Packages

```r
library(lme4)        # linear mixed models
library(lmerTest)    # p-values for LMM
library(rstatix)     # paired tests, effect sizes, anova_test
library(effectsize)  # Cohen's dz, partial eta²
library(car)         # Levene's test, VIF
library(ppcor)       # partial correlations
library(mgcv)        # GAMM for nonlinear trajectories
library(boot)        # bootstrap CIs for correlations
library(glmnet)      # LASSO / ridge regression (H15)
library(lmtest)      # Breusch-Pagan heteroscedasticity test (H15)
library(ggplot2)     # visualization
library(ggpubr)      # publication-ready stat plots
```

---

## 8. Caveats for Hypothesis Testing

1. **Crossover / within-subjects design** — use paired or mixed-model approaches with `participant_id` as random effect.

2. **Multiple comparisons** — testing 13+ HRV metrics × 8 epoch slots generates many simultaneous tests. Apply FDR correction (`p.adjust(method="BH")`) within each hypothesis family.

3. **Non-normality expected** — `RMSSD`, `HF_ms2`, `LF_ms2`, and `Total_power` are right-skewed by nature. Use log-transformed versions for parametric tests or default to Wilcoxon signed-rank. Confirm per Shapiro-Wilk before deciding.

4. **Unequal epoch counts** — even after trimming, subjects differ in total exercise/post epochs. The trimmed `master_trimmed` dataset standardises this for trajectory analysis, but the long-format file still reflects real observed data for mixed models (which handle unbalanced designs natively).

5. **Lactate N is smaller** — not all HRV participants have lactate data (27 lactate records vs 21/24 HRV files). Lactate correlations have lower power; interpret effect sizes alongside p-values.

6. **Exercise epoch trimming affects generalisability** — the trimmed dataset represents the shortest protocol's duration. Findings about exercise-phase HRV apply to that common window.

---

## 9. R Hypothesis Testing Script (`hrv_testing.R`)

### 9.1 How to Run

**Prerequisite:** run `hrv_pipeline.py` first so that `outputs/hrv_long_master.csv` exists.

```r
# Option A — from the R console or RStudio
setwd("C:/Users/dongu/Documents/UNAL_DMCH/HRV/UPN_RM_VPM")
source("hrv_testing.R")

```

### 9.2 What the Script Does (execution order)

| Step | Action |
|---|---|
| 1 | Loads `hrv_long_master.csv`; coerces `protocol` to factor, `phase_pos` to integer, `trimmed` to logical |
| 2 | Computes **orthostatic reactivity** as % change: `(stand_up − lying) / |lying| × 100` for all 4 orthostatic slots (pre_03, pre_05, post_03, post_05) |
| 3 | Runs **Shapiro-Wilk** on the paired difference (MR − VBT) per metric × epoch slot → normality decision for each test |
| 4 | Runs **Levene's test** on raw values at `lying_inicial` and `lying_final` → variance homogeneity check |
| 5 | Applies **decision tree** per metric: SW p > 0.05 → paired t-test + Cohen's dz; SW p ≤ 0.05 → Wilcoxon signed-rank + rank-biserial r; BH-FDR applied within each hypothesis family |
| 6 | Fits **Linear Mixed Models** (`lmerTest`) for exercise and post-exercise trajectories: `metric ~ protocol * phase_pos + (1 | participant_id)`, plus LRT for the interaction term |
| 7 | Exports all results as CSV files to `outputs/` |

[Results_description](docs/dictionary_results.md)
