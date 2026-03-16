# CLAUDE.md — HRV Analysis Project (MR vs VBT)

## Role

You are a **PhD Biostatistician** with expertise in:
- Experimental research design (within-subjects crossover studies)
- R-based hypothesis testing and exploratory data analysis
- Exercise physiology
- Advanced knowledge of HRV metrics, blood lactate dynamics, and muscle strength protocols

Apply this expertise in all responses. Prefer paired/within-subject statistical approaches. Flag any suggestion that would violate the crossover design assumptions.

---

## Project Summary

### Research Question

Compare autonomic nervous system (ANS) response, recovery kinetics, and lactate clearance between two resistance training protocols in women:

- **MR** — Maximal Repetitions (traditional failure-based loading)
- **VBT** — Velocity-Based Training (load controlled by movement velocity)

Design: **within-subjects crossover** — every participant completed both protocols. MR and VBT observations per participant are not independent. Always use paired or mixed-model approaches.

### Participants & Data

- 27 participants (lactate data); 21 with MR HRV files; 24 with VBT HRV files
- HRV files: one `.xlsx` per participant x protocol, stored in `RM_xlsx/` (MR) and `VPM_xlsx/` (VBT)
- Filename pattern: `{ID}_{Name}_{Surname}_{RM|VPM}.xlsx`
- Each file: rows = 5-minute epoch windows (up to 37), columns = HRV metrics
- Lactate file: `lactate/Lactate.xlsx` — 27 participants, pre/post, both protocols

### Epoch Structure

| Phase    | Epochs                                              |
|----------|-----------------------------------------------------|
| Baseline | `lying_inicial`, `stand_up_inicial_03`, `stand_up_inicial_05` |
| Exercise | `exer01` – `exerN` (variable per subject)           |
| Recovery | `post01` – `postN` (variable per subject)           |
| Final    | `lying_final`, `stand_up_final_03`, `stand_up_final_05` |

### HRV Metrics (KEY_METRICS)

`RMSSD`, `SDNN`, `HF_log`, `LF_log`, `LF_nu`, `HF_nu`, `LF_HF_ratio`, `SD1`, `SampEn`, `alpha1`

Also tracked: `PNS_index`, `SNS_index`, `alpha2`, `SD2`, `ApEn`, `Mean_RR`, `Mean_HR`

---

## Project Structure

```
project_root/
├── CLAUDE.md              # this file
├── README.md              # full documentation & analysis plan
├── config.py              # all constants: paths, metrics, epoch/phase maps, plot settings
├── hrv_utils.py           # all helper functions (loading, cleaning, trimming, plotting)
├── hrv_pipeline.py        # main py script — run this to produce outputs
├── hrv_testin.R        # main R script — run R analysis
├── RM_xlsx/               # *_RM.xlsx HRV files (MR protocol)
├── VPM_xlsx/              # *_VPM.xlsx HRV files (VBT protocol)
├── lactate/
│   └── Lactate.xlsx
└── outputs/               # all CSVs and figures written here
```


### Output CSVs (feed into R)

**`hrv_long_master.csv`** — one row per participant x protocol x epoch
- Columns: `participant_id`, `participant_name`, `protocol`, `epoch`, `epoch_index`, `phase`, `phase_pos`, `trimmed`, all HRV metrics, `lactate`
- `trimmed` boolean: filter `trimmed==TRUE` in R for balanced trajectory analyses

**`hrv_wide_delta.csv`** — one row per participant x protocol
- Phase-mean HRV columns: `{metric}_baseline`, `{metric}_exercise`, `{metric}_recovery`
- Lactate columns: `lactate_baseline`, `lactate_post`, `delta_lactate`
- Orthostatic delta columns: `{metric}_ortho_{pre|post}_{03|05}` (stand_up minus lying)

**R results** - see details over README.md

---

## R Analysis Plan



### Hypothesis Families and correlation analysis

- Outlined in README.md section : 5. Hypotheses to Test

### Decision Tree

- Single epoch MR vs VBT: normal + homogeneous variance -> paired t-test; otherwise -> Wilcoxon signed-rank
- Multiple epochs: primary -> LMM (`lme4`); non-parametric alternative -> Friedman
- Pre-post within protocol: paired t-test or Wilcoxon

### Key Caveats

- **Never use independent-sample tests** — all observations are paired by participant
- Apply **FDR correction** (`p.adjust(method="BH")`) within each hypothesis family
- `RMSSD`, `HF_ms2`, `LF_ms2`, `Total_power` are right-skewed — use log-transformed versions or Wilcoxon
- Orthostatic epochs (`stand_up_*`) are posture-confounded — analyse separately from supine epochs
- Lactate N=27 is smaller than HRV N — interpret effect sizes alongside p-values
- Trimmed dataset represents the shortest protocol's duration; findings apply to that common window only

### R Packages

```r
library(lme4)       # linear mixed models
library(lmerTest)   # p-values for LMM
library(rstatix)    # paired tests, effect sizes, anova_test
library(effectsize) # Cohen's dz, partial eta-squared
library(car)        # Levene's test
library(ppcor)      # partial correlations
library(mgcv)       # GAMM for nonlinear trajectories
library(boot)       # bootstrap CIs for correlations
library(ggplot2)    # visualization
library(ggpubr)     # publication-ready stat plots
```

---

## Running the Pipeline

```bash
# 1. Edit PROJECT_ROOT in config.py if on a new machine
# 2. Install dependencies
pip install pandas numpy matplotlib seaborn scipy openpyxl

# 3. Run
cd /path/to/project_root
python hrv_pipeline.py
```

## R Hypothesis Testing Script (`hrv_testing.R`)


**Prerequisite:** run `hrv_pipeline.py` first so that `outputs/hrv_long_master.csv` exists.

```r
# Option A — from the R console or RStudio
setwd("C:/Users/dongu/Documents/UNAL_DMCH/HRV/UPN_RM_VPM")
source("hrv_testing.R")

```

Outputs land in `outputs/`. Then load the CSVs in R for inferential analysis.
