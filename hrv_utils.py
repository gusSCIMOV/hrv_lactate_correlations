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

# ─────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────

def clean_col(c: str) -> str:
    """Strip whitespace/special chars from column names."""
    c = c.strip()
    c = re.sub(r"\s+", "_", c)
    c = re.sub(r"[^a-zA-Z0-9_]", "", c)
    c = re.sub(r"_+", "_", c).strip("_")
    return c

def clean_epoch(raw: str) -> str:
    """Normalise epoch label to canonical form; returns None for non-canonical rows."""
    if not isinstance(raw, str):
        return None
    s = raw.strip().lower()
    s = re.sub(r"\(\d+\)", "", s).strip().strip("_ ")
    s = re.sub(r"[-\s]+", "_", s)
    s = s.replace(":", "_")                                        # fix 'lying:final'
    s = re.sub(r"stand_up_inicial[_]?0?3", "stand_up_inicial_03", s)
    s = re.sub(r"stand_up_inicial[_]?0?5", "stand_up_inicial_05", s)
    s = re.sub(r"stand_up_final[_]?0?3",   "stand_up_final_03",   s)
    s = re.sub(r"stand_up_final[_]?0?5",   "stand_up_final_05",   s)
    m = re.match(r"(exer|post)(\d+)$", s)
    if m:
        s = f"{m.group(1)}{int(m.group(2)):02d}"
    # validate against canonical patterns; drop artefacts like 'sample_16'
    canonical = [
        r"^lying_inicial$", r"^lying_final$",
        r"^stand_up_inicial_(03|05)$", r"^stand_up_final_(03|05)$",
        r"^exer\d{2}$", r"^post\d{2}$"
    ]
    return s if any(re.match(p, s) for p in canonical) else None

def parse_filename(path: Path):
    """Return (participant_id, participant_name, protocol) from filename.
    Pattern: {ID}_{Name}_{Surname}_{RM|VPM}.xlsx  (all underscores, homogenised)
    """
    stem = path.stem                                              # e.g. '1_Laura_Guerrero_RM'
    protocol = "VBT" if re.search(r"_VPM$", stem, re.IGNORECASE) else "MR"
    core = re.sub(r"_(VPM|RM)$", "", stem, flags=re.IGNORECASE)  # '1_Laura_Guerrero'
    parts = core.split("_", 1)                                    # ['1', 'Laura_Guerrero']
    pid  = parts[0].strip()
    name = parts[1].strip().replace("_", " ") if len(parts) > 1 else pid
    return pid, name, protocol

def load_hrv_file(path: Path) -> pd.DataFrame:
    """Load one xlsx, drop AR_ and separator rows, clean columns."""
    raw = pd.read_excel(path, dtype=str)

    # deduplicate column names BEFORE cleaning (duplicates cause to_numeric to receive a DataFrame)
    seen = {}
    deduped = []
    for c in raw.columns:
        cc = clean_col(c)
        if cc in seen:
            seen[cc] += 1
            cc = f"{cc}_{seen[cc]}"
        else:
            seen[cc] = 0
        deduped.append(cc)
    raw.columns = deduped

    name_col = raw.columns[0]  # always 'names'

    # drop AR_ rows and blank separator rows
    mask_ar    = raw[name_col].str.strip().str.upper().str.startswith("AR_", na=False)
    mask_blank = raw[name_col].str.strip().isin(["", "nan"]) | raw[name_col].isna()
    raw = raw[~mask_ar & ~mask_blank].copy()

    # clean epoch names
    raw["epoch"] = raw[name_col].apply(clean_epoch)
    raw = raw[raw["epoch"].notna()].copy()

    # convert all HRV columns to numeric — iterate by position to be safe
    hrv_cols = [c for c in raw.columns if c not in [name_col, "epoch"]]
    for c in hrv_cols:
        col_data = raw[c]
        if isinstance(col_data, pd.DataFrame):   # shouldn't happen after dedup, but guard anyway
            col_data = col_data.iloc[:, 0]
        raw[c] = pd.to_numeric(col_data, errors="coerce")

    # drop columns that are entirely NaN (section-header artefacts)
    raw = raw.dropna(axis=1, how="all")

    pid, pname, protocol = parse_filename(path)
    raw.insert(0, "participant_id",   pid)
    raw.insert(1, "participant_name", pname)
    raw.insert(2, "protocol",         protocol)
    raw.drop(columns=[name_col], inplace=True)

    return raw

def resolve_col(df, aliases):
    for a in aliases:
        for c in df.columns:
            if c.lower() == a.lower() or c.lower().startswith(a.lower()):
                return c
    return None

def _ortho_delta(df, lying_ep, standup_ep, suffix, hrv_metrics):
    """Compute stand_up − lying per participant×protocol for each HRV metric."""
    
    HRV_METRICS=hrv_metrics
    ly = (df[df["epoch"]==lying_ep]
          .groupby(["participant_id","protocol"])[HRV_METRICS].mean()
          .reset_index())
    su = (df[df["epoch"]==standup_ep]
          .groupby(["participant_id","protocol"])[HRV_METRICS].mean()
          .reset_index())
    merged = ly.merge(su, on=["participant_id","protocol"], suffixes=("_ly","_su"))
    delta  = merged[["participant_id","protocol"]].copy()
    for m in HRV_METRICS:
        if f"{m}_ly" in merged and f"{m}_su" in merged:
            delta[f"{m}_ortho_{suffix}"] = merged[f"{m}_su"] - merged[f"{m}_ly"]
    return delta

# ─────────────────────────────────────────────
# SCALE CORRECTION
# ─────────────────────────────────────────────

SCALE_RULES = [
    (["RMSSD", "SDNN"],           160),
    (["HF_log", "LF_log"],         10),
    (["HF_nu", "LF_nu"],          100),
    (["LF_HF_ratio"],               4),
    (["SampEn", "alpha1"],          5),
]

def fix_metric_scales(df: pd.DataFrame) -> pd.DataFrame:
    """
    Correct scale/unit inconsistencies introduced by HRV export software.
    For each metric, values exceeding the threshold are divided by 10.
    Operates cell-wise so only affected cells are modified.
    Prints an audit count per metric.
    """
    df = df.copy()
    total_fixed = 0
    for metrics, threshold in SCALE_RULES:
        for m in metrics:
            if m not in df.columns:
                continue
            mask = df[m].notna() & (df[m] > threshold)
            n = mask.sum()
            if n:
                df.loc[mask, m] = df.loc[mask, m] / 10
                print(f"  scale fix: {m:15s} > {threshold:>5} → divided by 10  ({n} cells)")
                total_fixed += n
    if total_fixed == 0:
        print("  scale fix: no corrections needed")
    else:
        print(f"  scale fix: {total_fixed} total cells corrected")
    return df

# ─────────────────────────────────────────────
# STEP 1 — BUILD MASTER LONG-FORMAT DATAFRAME
# ─────────────────────────────────────────────

def build_master(rm_dir: Path, vpm_dir: Path) -> pd.DataFrame:
    frames = []
    for d in [rm_dir, vpm_dir]:
        if not d.exists():
            print(f"  ⚠️  Directory not found: {d}"); continue
        for f in sorted(d.glob("*.xlsx")):
            try:
                frames.append(load_hrv_file(f))
                print(f"  ✓  {f.name}")
            except Exception as e:
                print(f"  ✗  {f.name}: {e}")
    if not frames:
        raise RuntimeError("No files loaded. Check directory paths.")

    df = pd.concat(frames, ignore_index=True)

    # canonical epoch ordering + phase
    df["epoch_index"] = df["epoch"].map(
        {e: i+1 for i, e in enumerate(EPOCH_ORDER)}
    )
    df["phase"] = df["epoch"].map(PHASE_MAP).fillna("other")

    # keep only canonical epochs
    df = df[df["epoch_index"].notna()].copy()
    df["epoch_index"] = df["epoch_index"].astype(int)
    df.sort_values(["participant_id", "protocol", "epoch_index"], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df

def merge_lactate(df: pd.DataFrame, lactate_path: Path) -> pd.DataFrame:
    """Add lactate column; populated only at lying_inicial and lying_final.
    Joins on numeric participant_id (robust to accent/spelling differences in names).
    """
    lac = pd.read_excel(lactate_path)
    lac.columns = [c.strip() for c in lac.columns]
    lac["ID"] = lac["ID"].astype(str).str.strip()

    lac_map = {}
    for _, row in lac.iterrows():
        pid = row["ID"]
        lac_map[(pid, "MR",  "lying_inicial")] = row["RM_baseline_Lactate"]
        lac_map[(pid, "MR",  "lying_final")]   = row["RM_post_Lactate"]
        lac_map[(pid, "VBT", "lying_inicial")] = row["VPM_baseline_Lactate"]
        lac_map[(pid, "VBT", "lying_final")]   = row["VPM_post_Lactate"]

    df["lactate"] = df.apply(
        lambda r: lac_map.get((str(r["participant_id"]), r["protocol"], r["epoch"]), np.nan),
        axis=1
    )

    # audit: report any participant_ids with no lactate match
    matched_ids = {k[0] for k in lac_map}
    file_ids = set(df["participant_id"].astype(str).unique())
    unmatched = file_ids - matched_ids
    if unmatched:
        print(f"  ⚠️  No lactate match for participant_ids: {sorted(unmatched)}")
    return df

# ─────────────────────────────────────────────
# STEP 2 — EDA & VISUALIZATION
# ─────────────────────────────────────────────

def _assign_within_phase_pos(df, phase_prefix):
    """Add within-phase position (1-based, per subject×protocol) for exer or post epochs."""
    mask = df["epoch"].str.startswith(phase_prefix, na=False)
    sub  = df[mask].copy()
    # rank by epoch_index within each subject×protocol group
    sub["phase_pos"] = (sub.groupby(["participant_id","protocol"])["epoch_index"]
                           .rank(method="first").astype(int))
    return sub

def get_phase_nmin(df, phase_prefix):
    """Return the minimum epoch count across all subjects×protocols for this phase."""
    sub    = _assign_within_phase_pos(df, phase_prefix)
    counts = sub.groupby(["participant_id","protocol"])["phase_pos"].max()
    nmin   = int(counts.min()) if len(counts) else 0
    print(f"  Phase '{phase_prefix}*': counts per subject range "
          f"{int(counts.min())}–{int(counts.max())}, N_min={nmin}")
    if nmin < 2 * ANCHOR_N:
        print(f"  ⚠️  N_min={nmin} < 2×ANCHOR_N={2*ANCHOR_N}: "
              f"middle subsampling skipped, only anchors kept where possible.")
    return nmin, sub

def _select_phase_rows(phase_sub_with_pos, nmin):
    """
    Per subject×protocol: keep first ANCHOR_N + evenly-spaced middle + last ANCHOR_N.
    phase_sub_with_pos must have column 'phase_pos' (1-based within-phase rank).
    Returns index labels to keep.
    """
    keep_idx = []
    for (pid, proto), grp in phase_sub_with_pos.groupby(["participant_id","protocol"]):
        grp_sorted = grp.sort_values("phase_pos")
        positions  = grp_sorted["phase_pos"].tolist()
        n          = len(positions)

        if n <= 2 * ANCHOR_N:
            # not enough room for a middle — keep everything
            keep_idx.extend(grp_sorted.index.tolist())
            continue

        first6 = positions[:ANCHOR_N]
        last6  = positions[-ANCHOR_N:]
        middle = positions[ANCHOR_N: n - ANCHOR_N]

        # how many middle slots do we need?
        n_middle_target = nmin - 2 * ANCHOR_N
        if n_middle_target <= 0:
            # N_min too small: just first6 + last6
            selected_middle = []
        elif len(middle) <= n_middle_target:
            selected_middle = middle
        else:
            # evenly-spaced indices into middle list (real observed epochs, no interpolation)
            idx_float  = np.linspace(0, len(middle) - 1, n_middle_target)
            sel_int    = sorted(set(int(round(i)) for i in idx_float))
            selected_middle = [middle[i] for i in sel_int]

        kept_pos = set(first6 + selected_middle + last6)
        keep_idx.extend(grp_sorted[grp_sorted["phase_pos"].isin(kept_pos)].index.tolist())

    return keep_idx

def trim_exer_post(df):
    """
    Return df where exer and post phases are trimmed to N_min per phase:
      - first ANCHOR_N epochs  (per subject, within-phase position)
      - evenly subsampled middle to fill N_min - 2*ANCHOR_N slots
      - last  ANCHOR_N epochs
    Fixed epochs (lying/stand_up) are passed through unchanged.
    """
    keep_rows = []

    # fixed epochs — pass through
    fixed_mask = ~df["epoch"].str.startswith("exer", na=False) &                  ~df["epoch"].str.startswith("post", na=False)
    keep_rows.append(df[fixed_mask])

    for prefix in ["exer", "post"]:
        nmin, phase_sub = get_phase_nmin(df, prefix)
        if nmin == 0:
            continue
        selected = _select_phase_rows(phase_sub, nmin)
        keep_rows.append(df.loc[selected])

    trimmed = pd.concat(keep_rows).sort_values(
        ["participant_id","protocol","epoch_index"]).reset_index(drop=True)

    print(f"  trim_exer_post: {len(df)} → {len(trimmed)} rows  "
          f"(ANCHOR_N={ANCHOR_N}, exer/post trimmed to N_min per phase)")
    return trimmed

def _qq_slot(slot_key, slot_label, sub_df, key_metrics, suffix=""):

    km = [m for m in key_metrics if m in sub_df.columns and sub_df[m].notna().sum() > 3]
    if not km: return
    fig, axes = plt.subplots(2, len(km), figsize=(4*len(km), 8))
    if len(km) == 1: axes = np.array(axes).reshape(2,1)
    for j, m in enumerate(km):
        data = sub_df[m].dropna()
        axes[0,j].hist(data, bins=20, edgecolor="k", color="steelblue", alpha=0.7)
        axes[0,j].set_title(m); axes[0,j].set_xlabel("value")
        stats.probplot(data, plot=axes[1,j])
        axes[1,j].set_title(f"Q-Q {m}")
    fig.suptitle(f"Distribution & Q-Q — {slot_label}{suffix}", fontsize=13, y=1.01)
    plt.tight_layout()
    save_fig(fig, f"02_dist_qq_{slot_key}{suffix}.png")

def save_fig(fig, name):
    p = OUT_DIR / name
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓  {name}")

def _violin_slot(slot_key, slot_label, sub_df, key_metrics):
    km = [m for m in key_metrics if m in sub_df.columns and sub_df[m].notna().sum() > 3]
    if not km: return
    ncols = 4; nrows = int(np.ceil(len(km)/ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows), sharey=False)
    fig.patch.set_facecolor("white")
    axes = np.array(axes).flatten()
    for j, m in enumerate(km):
        ax = axes[j]
        ax.set_facecolor("white")
        sub_ok = sub_df[[m,"protocol"]].dropna()
        if sub_ok.empty: ax.set_visible(False); continue
        sns.violinplot(data=sub_ok, x="protocol", y=m, order=["MR","VBT"],
                       palette=PROTO_COLORS, inner="box", ax=ax, linewidth=1.2)
        sns.stripplot(data=sub_ok, x="protocol", y=m, order=["MR","VBT"],
                      color="white", edgecolor="grey", linewidth=0.5,
                      size=3, alpha=0.6, ax=ax, jitter=True)
        # ggplot-style spines: keep left + bottom only
        ax.spines[["top","right"]].set_visible(False)
        ax.spines[["left","bottom"]].set_color("#444444")
        #ax.yaxis.grid(True, color="#DDDDDD", linewidth=0.7, linestyle="--")
        ax.yaxis.grid(False)
        ax.set_axisbelow(True)
        ax.set_title(m, fontsize=10); ax.set_xlabel("")
    for a in axes[len(km):]: a.set_visible(False)
    fig.suptitle(f"HRV Distributions — {slot_label}", fontsize=13)
    plt.tight_layout()
    save_fig(fig, f"03_violin_{slot_key}.png")

# Robust: use median ± IQR; exclude participants with >40% missing in a metric

def robust_trajectory(df, metric, min_coverage=0.6):
    """
    Per epoch, compute median + IQR across subjects.
    Only include subjects that have data for ≥ min_coverage fraction of epochs.
    """
    # coverage filter
    n_epochs_per_sub = (df.groupby(["participant_id","protocol"])[metric]
                        .apply(lambda x: x.notna().mean()))
    good = n_epochs_per_sub[n_epochs_per_sub >= min_coverage].reset_index()[
        ["participant_id","protocol"]]
    filtered = df.merge(good, on=["participant_id","protocol"])
    
    agg = (filtered.groupby(["protocol","epoch_index","epoch","phase"])[metric]
           .agg(median="median", q25=lambda x: x.quantile(0.25),
                q75=lambda x: x.quantile(0.75), n="count")
           .reset_index())
    return agg

def plot_trajectory(metric, ax, master, min_coverage=0.6):
    agg = robust_trajectory(master, metric, min_coverage)
    if agg.empty: return
    for proto, grp in agg.groupby("protocol"):
        grp = grp.sort_values("epoch_index")
        ax.plot(grp["epoch_index"], grp["median"], label=proto,
                color=PROTO_COLORS[proto], lw=2)
        ax.fill_between(grp["epoch_index"], grp["q25"], grp["q75"],
                        alpha=0.2, color=PROTO_COLORS[proto])
    # phase shading
    for phase, color in [("exercise","#FFFACD"),("recovery","#E8F5E9")]:
        phase_idx = agg[agg["phase"]==phase]["epoch_index"]
        if not phase_idx.empty:
            ax.axvspan(phase_idx.min(), phase_idx.max(), alpha=0.15,
                       color=color, label=f"_{phase}")
    ax.set_title(metric); ax.set_xlabel("Epoch index"); ax.set_ylabel(metric)
    ax.legend(fontsize=8)

def _traj_fig(df, title, fname, traj_metrics, min_coverage=0.6):

    ncols = 2; nrows = int(np.ceil(len(traj_metrics)/ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4*nrows))
    axes = axes.flatten()
    for j, m in enumerate(traj_metrics):
        ax = axes[j]
        agg = robust_trajectory(df, m, min_coverage)
        if agg.empty: ax.set_visible(False); continue
        for proto, grp in agg.groupby("protocol"):
            grp = grp.sort_values("epoch_index")
            ax.plot(grp["epoch_index"], grp["median"], label=proto,
                    color=PROTO_COLORS[proto], lw=2)
            ax.fill_between(grp["epoch_index"], grp["q25"], grp["q75"],
                            alpha=0.2, color=PROTO_COLORS[proto])
        for phase, color in [("exercise","#FFFACD"),("recovery","#E8F5E9")]:
            pi = agg[agg["phase"]==phase]["epoch_index"]
            if not pi.empty:
                ax.axvspan(pi.min(), pi.max(), alpha=0.15, color=color, label=f"_{phase}")
        ax.set_title(m); ax.set_xlabel("Epoch index"); ax.set_ylabel(m)
        ax.legend(fontsize=8)
    for a in axes[len(traj_metrics):]: a.set_visible(False)
    fig.suptitle(title, fontsize=13, y=1.01)
    plt.tight_layout()
    save_fig(fig, fname)

def _boxplot_traj_fig(df, title, fname, traj_metrics, min_coverage=0.6):
    """
    Per metric: narrow side-by-side boxplots at every epoch_index,
    with dotted line connecting medians per protocol.
    Boxes with N < MIN_N_BOX are skipped.
    """
    # apply same coverage filter as robust_trajectory
    ncols = 2; nrows = int(np.ceil(len(traj_metrics)/ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 4*nrows))
    fig.patch.set_facecolor("white")
    axes = axes.flatten()

    offsets = {"MR": -BOX_OFFSET, "VBT": BOX_OFFSET}

    for j, m in enumerate(traj_metrics):
        ax = axes[j]
        ax.set_facecolor("white")
        ax.spines[["top","right"]].set_visible(False)
        ax.spines[["left","bottom"]].set_color("#444444")
        ax.yaxis.grid(True, color="#DDDDDD", linewidth=0.7, linestyle="--")
        ax.set_axisbelow(True)

        # coverage filter (same logic as robust_trajectory)
        coverage = (df.groupby(["participant_id","protocol"])[m]
                    .apply(lambda x: x.notna().mean()))
        good = coverage[coverage >= min_coverage].reset_index()[["participant_id","protocol"]]
        filtered = df.merge(good, on=["participant_id","protocol"])

        epoch_idx_all = sorted(filtered["epoch_index"].unique())
        plotted_any = False

        for proto in ["MR","VBT"]:
            sub_p   = filtered[filtered["protocol"]==proto]
            medians = []
            x_pos   = []

            for ei in epoch_idx_all:
                vals = sub_p[sub_p["epoch_index"]==ei][m].dropna()
                if len(vals) < MIN_N_BOX:
                    continue
                x = ei + offsets[proto]
                x_pos.append(x)
                medians.append(vals.median())

                # boxplot stats manually for narrow box
                q1, med, q3 = vals.quantile([0.25, 0.5, 0.75])
                iqr  = q3 - q1
                wlo  = vals[vals >= q1 - 1.5*iqr].min()
                whi  = vals[vals <= q3 + 1.5*iqr].max()
                col  = PROTO_COLORS[proto]

                # box
                ax.add_patch(plt.Rectangle(
                    (x - BOX_WIDTH/2, q1), BOX_WIDTH, iqr,
                    facecolor=col, alpha=0.35, edgecolor=col, linewidth=0.8))
                # median line
                ax.plot([x - BOX_WIDTH/2, x + BOX_WIDTH/2], [med, med],
                        color=col, lw=1.5, zorder=3)
                # whiskers
                ax.plot([x, x], [wlo, q1], color=col, lw=0.8)
                ax.plot([x, x], [q3, whi], color=col, lw=0.8)

            # dotted trend line connecting medians
            if len(x_pos) > 1:
                ax.plot(x_pos, medians, color=PROTO_COLORS[proto],
                        lw=1.8, ls=":", marker="o", markersize=3,
                        label=proto, zorder=4)
            plotted_any = True

        if not plotted_any:
            ax.set_visible(False); continue

        # phase shading (use full df for reference)
        for phase, color in [("exercise","#FFFACD"),("recovery","#E8F5E9")]:
            pi = df[df["phase"]==phase]["epoch_index"]
            if not pi.empty:
                ax.axvspan(pi.min()-0.5, pi.max()+0.5, alpha=0.12, color=color)

        ax.set_title(m, fontsize=10)
        ax.set_xlabel("Epoch index"); ax.set_ylabel(m)
        ax.legend(fontsize=8)

    for a in axes[len(traj_metrics):]: a.set_visible(False)
    fig.suptitle(title, fontsize=13, y=1.01)
    plt.tight_layout()
    save_fig(fig, fname)

from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

def _corr_matrix(sub_df, title, fname, key_metrics):
    """Participant-median aggregated Spearman correlation, upper triangle + dendrogram."""
    km = [m for m in key_metrics if m in sub_df.columns and sub_df[m].notna().sum() > 5]
    if len(km) < 3: return
    # aggregate: one row per participant × protocol (median across epochs)
    agg = sub_df.groupby(["participant_id","protocol"])[km].median()
    corr = agg[km].corr(method="spearman")
    # hierarchical clustering on distance = 1 - |r|
    dist = squareform(1 - corr.abs().values, checks=False)
    dist = np.clip(dist, 0, None)
    linkage = hierarchy.linkage(dist, method="average")
    order   = hierarchy.leaves_list(linkage)
    corr_ord = corr.iloc[order, order]

    fig = plt.figure(figsize=(max(8, len(km)*0.7+2), max(7, len(km)*0.7)))
    gs  = fig.add_gridspec(2, 2, width_ratios=[1,5], height_ratios=[5,1],
                           hspace=0.02, wspace=0.02)
    ax_dend = fig.add_subplot(gs[0,0])
    ax_heat = fig.add_subplot(gs[0,1])

    with plt.style.context("default"):
        hierarchy.dendrogram(linkage, ax=ax_dend, orientation="left",
                             labels=corr.index[order].tolist(),
                             leaf_font_size=8, color_threshold=0)
    ax_dend.axis("off")

    # upper triangle mask
    mask = np.tril(np.ones_like(corr_ord, dtype=bool), k=-1)
    masked = corr_ord.where(~mask)
    sns.heatmap(masked, ax=ax_heat, cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                annot=True, fmt=".2f", annot_kws={"size":7},
                linewidths=0.3, square=True, cbar_kws={"shrink":0.6},
                xticklabels=corr_ord.columns, yticklabels=corr_ord.index)
    ax_heat.set_xticklabels(ax_heat.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    ax_heat.set_yticklabels(ax_heat.get_yticklabels(), rotation=0, fontsize=8)
    fig.suptitle(title, fontsize=12, y=1.01)
    plt.tight_layout()
    save_fig(fig, fname)

from matplotlib.patches import Ellipse
import matplotlib.transforms as mtransforms

def _confidence_ellipse(x, y, ax, color, n_std=1.5):
    """Draw a covariance ellipse for protocol cluster boundary."""
    if len(x) < 3: return
    cov = np.cov(x, y)
    if cov.shape != (2,2): return
    pearson = cov[0,1] / (np.sqrt(cov[0,0]) * np.sqrt(cov[1,1]) + 1e-9)
    rx, ry = np.sqrt(1 + pearson), np.sqrt(1 - pearson)
    ellipse = Ellipse((0,0), width=rx*2, height=ry*2,
                      facecolor="none", edgecolor=color, linewidth=1.5, linestyle="--", alpha=0.7)
    scale_x = np.sqrt(cov[0,0]) * n_std
    scale_y = np.sqrt(cov[1,1]) * n_std
    transf = (mtransforms.Affine2D()
              .rotate_deg(45)
              .scale(scale_x, scale_y)
              .translate(np.mean(x), np.mean(y)))
    ellipse.set_transform(transf + ax.transData)
    ax.add_patch(ellipse)

def _lactate_scatter(metrics, slot_key, slot_label, delta_df):
    """
    delta_df: columns = participant_id, protocol, {metric}_delta (or mean for pooled),
              lactate_post, delta_lactate
    """
    km = [m for m in metrics if f"{m}_delta" in delta_df.columns or f"{m}_mean" in delta_df.columns]
    if not km or delta_df.empty: return
    for x_col, xlabel in [("lactate_post","Post Lactate (mmol/L)"),
                           ("delta_lactate","ΔLactate (post−pre)")]:
        ncols = 4; nrows = int(np.ceil(len(km)/ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
        axes = np.array(axes).flatten()
        for j, m in enumerate(km):
            ax = axes[j]
            y_col = f"{m}_delta" if f"{m}_delta" in delta_df.columns else f"{m}_mean"
            y_label = f"Δ{m}" if "_delta" in y_col else f"{m} mean"
            sub = delta_df[[x_col, y_col, "protocol"]].dropna()
            for proto, grp in sub.groupby("protocol"):
                ax.scatter(grp[x_col], grp[y_col],
                           color=PROTO_COLORS[proto], label=proto, s=50, alpha=0.8, zorder=3)
                if len(grp) > 3:
                    coefs = np.polyfit(grp[x_col].values, grp[y_col].values, 1)
                    xr = np.linspace(grp[x_col].min(), grp[x_col].max(), 50)
                    ax.plot(xr, np.polyval(coefs, xr),
                            color=PROTO_COLORS[proto], lw=1.5, ls="--")
                    _confidence_ellipse(grp[x_col].values, grp[y_col].values,
                                        ax, PROTO_COLORS[proto])
            ax.axhline(0, color="k", lw=0.7, ls=":")
            ax.set_xlabel(xlabel); ax.set_ylabel(y_label); ax.set_title(m, fontsize=9)
            ax.legend(fontsize=7)
        for a in axes[len(km):]: a.set_visible(False)
        fig.suptitle(f"Lactate vs ΔHRV — {slot_label} | x={xlabel}", fontsize=12)
        plt.tight_layout()
        save_fig(fig, f"07_lactate_{slot_key}_{x_col}.png")
