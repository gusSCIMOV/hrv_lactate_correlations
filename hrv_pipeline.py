"""
HRV Analysis Pipeline — MR vs VBT Protocols
Steps 1 (ingestion/cleaning/aggregation) + 2 (EDA/visualization)
Assumes project structure:
  project_root/
    RM_xlsx/     → *_RM.xlsx files
    VPM_xlsx/    → *_VPM.xlsx files
    Lactate.xlsx
    
"""

import re
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from pathlib import Path
from scipy import stats

from config import *
from hrv_utils import resolve_col, trim_exer_post
from hrv_utils import build_master, merge_lactate, fix_metric_scales
from hrv_utils import _ortho_delta, _qq_slot, save_fig, _corr_matrix
from hrv_utils import _violin_slot, _traj_fig, _boxplot_traj_fig, _lactate_scatter


warnings.filterwarnings("ignore")
sns.set_theme(style="whitegrid", palette="Set2", font_scale=1.1)


key_metrics = KEY_METRICS

# ─────────────────────────────────────────────
# STEP 1 — RUN
# ─────────────────────────────────────────────
print("=== Loading HRV files ===")
master = build_master(RM_DIR, VPM_DIR)
print(f"\nRaw shape: {master.shape}")

print("\n=== Merging lactate ===")
master = merge_lactate(master, LACTATE_FILE)

col_map = {k: resolve_col(master, v) for k, v in COL_ALIASES.items()}
print("\n=== Column mapping ===")
for k, v in col_map.items():
    print(f"  {k:15s} → {v}")

# rename to canonical names
rename = {v: k for k, v in col_map.items() if v and v != k}
master.rename(columns=rename, inplace=True)
HRV_METRICS = [k for k, v in col_map.items() if v is not None]

print("\n=== Correcting metric scales ===")
master = fix_metric_scales(master)

# ── Add phase_pos (within-phase position, 1-based) to master ─────────────────
for prefix in ["exer", "post"]:
    mask = master["epoch"].str.startswith(prefix, na=False)
    master.loc[mask, "phase_pos"] = (
        master[mask].groupby(["participant_id","protocol"])["epoch_index"]
        .rank(method="first").astype(int)
    )
# fixed epochs get phase_pos = 1 (each is a single point in its slot)
master["phase_pos"] = master["phase_pos"].fillna(1).astype(int)
master["trimmed"] = False   # updated after master_trimmed is created

# export master (trimmed flag will be correct after next block)
master.to_csv(OUT_DIR / "hrv_long_master.csv", index=False)
print(f"\n✓  Saved hrv_long_master.csv  shape={master.shape}")

# wide delta export
phase_agg = (master.groupby(["participant_id", "participant_name", "protocol", "phase"])[HRV_METRICS]
             .mean().reset_index())
wide = phase_agg.pivot(index=["participant_id","participant_name","protocol"],
                       columns="phase", values=HRV_METRICS)
wide.columns = [f"{m}_{p}" for m, p in wide.columns]
wide.reset_index(inplace=True)

# add lactate columns
lac_cols = master[master["epoch"].isin(["lying_inicial","lying_final"])][
    ["participant_id","protocol","epoch","lactate"]].dropna(subset=["lactate"])
lac_wide = lac_cols.pivot_table(index=["participant_id","protocol"],
                                columns="epoch", values="lactate").reset_index()

# rename explicitly by epoch name 
lac_wide = lac_wide.rename(columns={
    "lying_inicial": "lactate_baseline",
    "lying_final":   "lactate_post"
})
lac_wide["delta_lactate"] = lac_wide["lactate_post"] - lac_wide["lactate_baseline"]
wide = wide.merge(lac_wide, on=["participant_id","protocol"], how="left")

# ── Add orthostatic delta columns to wide ────────────────────────────────────
# delta = stand_up value − lying value, for pre and post exercise separately

ortho_pre_03  = _ortho_delta(master, "lying_inicial",    "stand_up_inicial_03", "pre_03", HRV_METRICS)
ortho_pre_05  = _ortho_delta(master, "lying_inicial",    "stand_up_inicial_05", "pre_05", HRV_METRICS)
ortho_post_03 = _ortho_delta(master, "lying_final",      "stand_up_final_03",   "post_03", HRV_METRICS)
ortho_post_05 = _ortho_delta(master, "lying_final",      "stand_up_final_05",   "post_05", HRV_METRICS)

for ortho_df in [ortho_pre_03, ortho_pre_05, ortho_post_03, ortho_post_05]:
    wide = wide.merge(ortho_df, on=["participant_id","protocol"], how="left")

wide.to_csv(OUT_DIR / "hrv_wide_delta.csv", index=False)
print(f"✓  Saved hrv_wide_delta.csv    shape={wide.shape}")


# ─────────────────────────────────────────────
# STEP 2 — EDA & VISUALIZATION
# ─────────────────────────────────────────────
key_metrics = [m for m in key_metrics if m in master.columns]

traj_metrics = key_metrics
traj_metrics = [m for m in traj_metrics if m in master.columns]


master_trimmed = trim_exer_post(master)

# ── Fill trimmed flag now that master_trimmed exists ─────────────────────────
trimmed_idx = set(zip(master_trimmed["participant_id"],
                      master_trimmed["protocol"],
                      master_trimmed["epoch_index"]))
master["trimmed"] = master.apply(
    lambda r: (r["participant_id"], r["protocol"], r["epoch_index"]) in trimmed_idx,
    axis=1
)
master.to_csv(OUT_DIR / "hrv_long_master.csv", index=False)   # re-export with flag
print(f"✓  Re-saved hrv_long_master.csv with trimmed flag  shape={master.shape}")

# ── 2.1  Missingness heatmap ──────────────────
print("\n=== Plotting ===")
fig, ax = plt.subplots(figsize=(14, 6))
miss = (master.groupby(["epoch","protocol"])[HRV_METRICS]
        .apply(lambda x: x.isna().mean()).reset_index()
        .melt(id_vars=["epoch","protocol"], var_name="metric", value_name="pct_missing"))
pivot_m = miss[miss["protocol"]=="MR"].pivot(index="metric", columns="epoch", values="pct_missing")
sns.heatmap(pivot_m, ax=ax, cmap="Reds", linewidths=0.3, cbar_kws={"label":"Prop. missing"})
ax.set_title("Missingness by Epoch — MR protocol"); ax.set_xlabel(""); ax.set_ylabel("")
save_fig(fig, "01_missingness_MR.png")

# ── 2.2  Distribution + Q-Q — 8 epoch slots + pooled exer/post ──────────────

# slots: 6 fixed epochs + pooled exercise (trimmed) + pooled post (trimmed)
QQ_SLOTS = [(e, e, master[master["epoch"]==e])          for e in FIXED_EPOCHS] + [
            ("exercise_pooled", "Exercise (trimmed)",
             master_trimmed[master_trimmed["phase"]=="exercise"]),
            ("post_pooled",     "Post (trimmed)",
             master_trimmed[master_trimmed["phase"]=="recovery"]
             [master_trimmed[master_trimmed["phase"]=="recovery"]["epoch"].str.startswith("post",na=False)])]


for slot_key, slot_label, sub_df in QQ_SLOTS:
    _qq_slot(slot_key, slot_label, sub_df, key_metrics)
# also trimmed vs untrimmed for pooled phases
for phase_key, phase_label, phase_str in [
        ("exercise_pooled","Exercise","exercise"),
        ("post_pooled","Post","recovery")]:
    sub_un = master[master["phase"]==phase_str]
    if phase_str == "recovery":
        sub_un = sub_un[sub_un["epoch"].str.startswith("post",na=False)]
    _qq_slot(phase_key, phase_label, sub_un, key_metrics, suffix="_untrimmed")

# ── 2.3  Violins — ggplot style, 8 epoch slots + trimmed/untrimmed pooled ─────
VIOLIN_SLOTS = [(e, e, master[master["epoch"]==e]) for e in FIXED_EPOCHS] + [
    ("exercise_trimmed",   "Exercise (trimmed)",   master_trimmed[master_trimmed["phase"]=="exercise"]),
    ("post_trimmed",       "Post (trimmed)",        master_trimmed[master_trimmed["phase"]=="recovery"]
                                                    [master_trimmed[master_trimmed["phase"]=="recovery"]["epoch"].str.startswith("post",na=False)]),
    ("exercise_untrimmed", "Exercise (untrimmed)",  master[master["phase"]=="exercise"]),
    ("post_untrimmed",     "Post (untrimmed)",      master[master["phase"]=="recovery"]
                                                    [master[master["phase"]=="recovery"]["epoch"].str.startswith("post",na=False)]),
]

for slot_key, slot_label, sub_df in VIOLIN_SLOTS:
    _violin_slot(slot_key, slot_label, sub_df, key_metrics)

# ── 2.3b  Lactate violin — pre / post / delta, MR vs VBT ────────────────────
if not lac_wide.empty:
    lac_plot = lac_wide.copy()
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.patch.set_facecolor("white")
    for ax, col, label in zip(axes,
                               ["lactate_baseline", "lactate_post", "delta_lactate"],
                               ["Lactate Pre (mmol/L)", "Lactate Post (mmol/L)",
                                "ΔLactate Post−Pre (mmol/L)"]):
        sub_ok = lac_plot[["protocol", col]].dropna()
        ax.set_facecolor("white")
        sns.violinplot(data=sub_ok, x="protocol", y=col, order=["MR","VBT"],
                       palette=PROTO_COLORS, inner="box", ax=ax, linewidth=1.2)
        sns.stripplot(data=sub_ok, x="protocol", y=col, order=["MR","VBT"],
                      color="white", edgecolor="grey", linewidth=0.5,
                      size=5, alpha=0.7, ax=ax, jitter=True)
        ax.spines[["top","right"]].set_visible(False)
        ax.spines[["left","bottom"]].set_color("#444444")
        ax.yaxis.grid(True, color="#DDDDDD", linewidth=0.7, linestyle="--")
        ax.set_ylim(0,15)
        ax.yaxis.grid(False)
        ax.set_axisbelow(True)
        ax.set_title(label, fontsize=11); ax.set_xlabel("")
        # annotate N per protocol
        for i, proto in enumerate(["MR","VBT"]):
            n = sub_ok[sub_ok["protocol"]==proto][col].notna().sum()
            ax.text(i, ax.get_ylim()[0], f"N={n}", ha="center", va="bottom",
                    fontsize=8, color="grey")
    fig.suptitle("Blood Lactate — MR vs VBT", fontsize=13)
    plt.tight_layout()
    save_fig(fig, "09_violin_lactate.png")

# ── 2.4  Temporal trajectory — ROBUST aggregation ──


_traj_fig(master,         "Temporal Trajectories — Untrimmed (≥60% coverage)",
          "04_temporal_trajectories_untrimmed.png", traj_metrics)

_traj_fig(master_trimmed, "Temporal Trajectories — Trimmed (≥60% coverage)",
          "04_temporal_trajectories_trimmed.png",traj_metrics)

# ── 2.4b  Boxplot trajectories — narrow boxes + median trend line ─────────────

_boxplot_traj_fig(master,
                  "Boxplot Trajectories — Untrimmed (N≥3, ≥60% coverage)",
                  "04b_boxplot_traj_untrimmed.png",traj_metrics)
_boxplot_traj_fig(master_trimmed,
                  "Boxplot Trajectories — Trimmed (N≥3, ≥60% coverage)",
                  "04b_boxplot_traj_trimmed.png",traj_metrics)

# ── 2.5  Pairplots — all 8 epoch slots ───────────────────────────────────────

PAIRPLOT_SLOTS = [(e, e, master[master["epoch"]==e]) for e in FIXED_EPOCHS] + [
    ("exercise_trimmed", "Exercise (trimmed)",
     master_trimmed[master_trimmed["phase"]=="exercise"]),
    ("post_trimmed", "Post (trimmed)",
     master_trimmed[master_trimmed["phase"]=="recovery"]
     [master_trimmed[master_trimmed["phase"]=="recovery"]["epoch"].str.startswith("post",na=False)]),
]
for slot_key, slot_label, sub_df in PAIRPLOT_SLOTS:
    pp_cols = [m for m in key_metrics if m in sub_df.columns]
    sub_pp = sub_df[pp_cols + ["protocol"]].dropna()
    if len(sub_pp) > 3 and len(pp_cols) > 1:
        g = sns.pairplot(sub_pp, hue="protocol", palette=PROTO_COLORS,
                         plot_kws={"alpha":0.6, "s":30}, diag_kind="kde")
        g.fig.suptitle(f"Pairplot — {slot_label}", y=1.02)
        save_fig(g.fig, f"05_pairplot_{slot_key}.png")

# ── 2.6  Correlation matrices — Spearman, upper triangle + dendrogram ──────

# Epoch-specific slots: one obs per participant → straightforward
for ep in FIXED_EPOCHS:
    for proto in ["MR","VBT"]:
        sub = master[(master["epoch"]==ep) & (master["protocol"]==proto)]
        _corr_matrix(sub, f"Spearman Corr — {ep} | {proto}", f"06_corr_{ep}_{proto}.png",key_metrics)

# Full temporal dynamic: participant-median across ALL epochs per protocol
for proto in ["MR","VBT"]:
    sub = master[master["protocol"]==proto]
    _corr_matrix(sub, f"Spearman Corr — All epochs (participant median) | {proto}",
                 f"06_corr_all_epochs_{proto}.png",key_metrics)

# ── 2.7  Lactate vs ΔHRV — all matched pairs + pooled phases ─────────────────

# Build delta dataframes for matched pairs
for pre_ep, post_ep in MATCHED_PAIRS:
    pair_label = f"{pre_ep}_vs_{post_ep}"
    pair_key   = pre_ep.replace("_inicial","").replace("stand_up_","su")
    pre_df  = master[master["epoch"]==pre_ep ][["participant_id","protocol"] + key_metrics].copy()
    post_df = master[master["epoch"]==post_ep][["participant_id","protocol"] + key_metrics].copy()
    merged  = pre_df.merge(post_df, on=["participant_id","protocol"],
                           suffixes=("_pre","_post"))
    delta_df = merged[["participant_id","protocol"]].copy()
    for m in key_metrics:
        if f"{m}_pre" in merged and f"{m}_post" in merged:
            delta_df[f"{m}_delta"] = merged[f"{m}_post"] - merged[f"{m}_pre"]
    delta_df = delta_df.merge(lac_wide[["participant_id","protocol",
                                         "lactate_post","delta_lactate"]],
                               on=["participant_id","protocol"], how="left")
    _lactate_scatter(key_metrics, pair_key, f"{pre_ep} → {post_ep}", delta_df)

# Pooled exercise and post: use participant mean per protocol (no paired delta)
for phase_key, phase_label, phase_df in [
        ("exercise_trimmed", "Exercise (trimmed)",
         master_trimmed[master_trimmed["phase"]=="exercise"]),
        ("post_trimmed", "Post (trimmed)",
         master_trimmed[master_trimmed["phase"]=="recovery"]
         [master_trimmed[master_trimmed["phase"]=="recovery"]["epoch"].str.startswith("post",na=False)])]:
    pmean = (phase_df.groupby(["participant_id","protocol"])[key_metrics]
             .median().reset_index())
    pmean.columns = ["participant_id","protocol"] + [f"{m}_mean" for m in key_metrics
                                                      if m in phase_df.columns]
    pmean = pmean.merge(lac_wide[["participant_id","protocol",
                                   "lactate_post","delta_lactate"]],
                        on=["participant_id","protocol"], how="left")
    _lactate_scatter(key_metrics, phase_key, phase_label, pmean)

# ── 2.8  [removed — superseded by 2.3 per-slot violins] ────────────────────

print("\n✅  All outputs saved to:", OUT_DIR)
print("    Next: open R and load hrv_long_master.csv / hrv_wide_delta.csv")
