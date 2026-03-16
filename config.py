"""
config.py — Project-wide configuration for HRV Analysis Pipeline
Edit only this file when moving between machines or changing project structure.
Import in any module with: from config import *
"""

from pathlib import Path


# ─────────────────────────────────────────────────────────────────────────────
# KEY HRV METRICS
# ─────────────────────────────────────────────────────────────────────────────

KEY_METRICS = ["RMSSD","SDNN","HF_log","LF_log","LF_nu","HF_nu","LF_HF_ratio","SD1","SampEn","alpha1"]

# ─────────────────────────────────────────────────────────────────────────────
# PATHS  ← only section that needs editing per machine
# ─────────────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path("C:/Users/dongu/Documents/UNAL_DMCH/HRV/UPN_RM_VPM")
RM_DIR       = PROJECT_ROOT / "RM_xlsx"
VPM_DIR      = PROJECT_ROOT / "VPM_xlsx"
LACTATE_FILE = PROJECT_ROOT / "lactate/Lactate.xlsx"
OUT_DIR      = PROJECT_ROOT / "outputs"
OUT_DIR.mkdir(exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# CANONICAL EPOCH ORDER & PHASE MAP
# ─────────────────────────────────────────────────────────────────────────────
EPOCH_ORDER = (
    ["lying_inicial"] +
    [f"stand_up_inicial_{m:02d}" for m in [3, 5]] +
    [f"exer{i:02d}" for i in range(1, 25)] +
    [f"post{i:02d}" for i in range(1, 12)] +
    ["lying_final"] +
    [f"stand_up_final_{m:02d}" for m in [3, 5]]
)

PHASE_MAP = {e: "baseline" for e in
             ["lying_inicial", "stand_up_inicial_03", "stand_up_inicial_05"]}
PHASE_MAP.update({f"exer{i:02d}": "exercise" for i in range(1, 25)})
PHASE_MAP.update({f"post{i:02d}": "recovery" for i in range(1, 12)})
PHASE_MAP.update({e: "recovery" for e in
                  ["lying_final", "stand_up_final_03", "stand_up_final_05"]})

# ─────────────────────────────────────────────────────────────────────────────
# NAMED EPOCH SLOTS
# ─────────────────────────────────────────────────────────────────────────────
FIXED_EPOCHS = [
    "lying_inicial", "stand_up_inicial_03", "stand_up_inicial_05",
    "lying_final",   "stand_up_final_03",   "stand_up_final_05",
]

MATCHED_PAIRS = [
    ("lying_inicial",       "lying_final"),
    ("stand_up_inicial_03", "stand_up_final_03"),
    ("stand_up_inicial_05", "stand_up_final_05"),
]

PHASE_ORDER = ["baseline", "exercise", "recovery"]

# ─────────────────────────────────────────────────────────────────────────────
# HRV METRIC ALIASES  (maps canonical name → possible raw column names)
# ─────────────────────────────────────────────────────────────────────────────
# column name audit — select core HRV metrics robustly
COL_ALIASES = {
    "RMSSD":      ["RMSSD_ms",   "RMSSDms",   "RMSSD"],
    "SDNN":       ["SDNN_ms",    "SDNNms",     "SDNN"],
    "Mean_RR":    ["Mean_RR_ms", "MeanRR_ms",  "Mean_RR_ms_"],
    "Mean_HR":    ["Mean_HR_beatsmin", "Mean_HR"],
    "LF_HF_ratio":      ["LFHF_ratio", "LF_HF_ratio", "LF_HF"], 
    "HF_log":     ["HF_log",     "HFlog",     "HFLog"],
    "LF_log":     ["LF_log",     "LFlog",     "LFLog"],
    "LF_nu":      ["LF_nu_",      "LFnu",      "LF_nu"],
    "HF_nu":      ["HF_nu_",      "HFnu",      "HF_nu"],
    "SD1":        ["SD1_ms_",     "SD1ms",     "SD1"],
    "SD2":        ["SD2_ms_",     "SD2ms",     "SD2"],
    "SampEn":     ["Sample_entropy_SampEn", "SampEn"],
    "ApEn":       ["Approximate_entropy_ApEn", "ApEn"],
    "alpha1":     ["alpha_1","alpha_1_",    "alpha1"],
    "PNS_index":  ["PNS_index_",  "PNS_index"],
    "SNS_index":  ["SNS_index_",  "SNS_index"],
}
# ─────────────────────────────────────────────────────────────────────────────
# TRIMMING
# ─────────────────────────────────────────────────────────────────────────────
ANCHOR_N = 6          # epochs to preserve at each end of exercise/post phases

# ─────────────────────────────────────────────────────────────────────────────
# VISUALIZATION
# ─────────────────────────────────────────────────────────────────────────────
PROTO_COLORS = {"MR": "#E64B35", "VBT": "#4DBBD5"}

BOX_OFFSET = 0.3      # left/right nudge per protocol in boxplot trajectories
BOX_WIDTH  = 0.25     # width of each narrow box
MIN_N_BOX  = 3        # minimum observations to draw a box at an epoch
