# hrv_testing.R ──────────────────────────────────────────────────────────────
# Between-Protocol Hypothesis Testing  H1 – H10
# MR vs VBT Resistance Training — HRV & Blood Lactate
# ─────────────────────────────────────────────────────────────────────────────
# Input : outputs/hrv_long_master.csv
# Output: outputs/results_*.csv
# R 4.5.1 (2025-06-13 ucrt)  x86_64-w64-mingw32
# ─────────────────────────────────────────────────────────────────────────────


# ── 0. SETUP ──────────────────────────────────────────────────────────────────

if (!requireNamespace("lmerTest", quietly = TRUE)) install.packages("lmerTest")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(rstatix)
  library(car)
  library(lme4)
  library(lmerTest)
  library(broom)
})

ROOT        <- "C:/Users/dongu/Documents/UNAL_DMCH/HRV/UPN_RM_VPM"
OUT_DIR     <- file.path(ROOT, "outputs")

KEY_METRICS <- c("RMSSD", "SDNN", "HF_log", "LF_log", "LF_nu", "HF_nu",
                 "LF_HF_ratio", "SD1", "SampEn", "alpha1")

ALPHA       <- 0.05   # significance threshold


# ── 1. LOAD & PREPARE ─────────────────────────────────────────────────────────

master <- read.csv(file.path(OUT_DIR, "hrv_long_master.csv"),
                   stringsAsFactors = FALSE)

master$participant_id <- as.character(master$participant_id)
master$protocol       <- factor(master$protocol, levels = c("MR", "VBT"))
master$phase_pos      <- as.integer(master$phase_pos)
master$trimmed        <- master$trimmed %in% c("True", "TRUE", TRUE)

# restrict to key metrics that are actually present
key_metrics <- KEY_METRICS[KEY_METRICS %in% names(master)]

cat("═══════════════════════════════════════════════════════════════\n")
cat("HRV Testing Pipeline — H1 to H10\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat(sprintf("Loaded : %d rows  |  %d participants  |  %d key metrics\n",
            nrow(master),
            length(unique(master$participant_id)),
            length(key_metrics)))
cat("Metrics:", paste(key_metrics, collapse = ", "), "\n\n")


# ── 2. ORTHOSTATIC REACTIVITY  (% change vs lying) ────────────────────────────
#
#   % change = (stand_up − lying) / |lying| × 100
#   computed per participant × protocol before any testing
#   |lying| in denominator handles negative metric values (e.g. PNS_index)
# ─────────────────────────────────────────────────────────────────────────────

compute_ortho_pct <- function(df, lying_ep, standup_ep) {
  ly <- df[df$epoch == lying_ep,
           c("participant_id", "protocol", key_metrics)]
  su <- df[df$epoch == standup_ep,
           c("participant_id", "protocol", key_metrics)]

  merged <- merge(ly, su,
                  by     = c("participant_id", "protocol"),
                  suffixes = c("_ly", "_su"))

  out <- merged[, c("participant_id", "protocol")]
  for (m in key_metrics) {
    ly_col <- paste0(m, "_ly"); su_col <- paste0(m, "_su")
    if (ly_col %in% names(merged) && su_col %in% names(merged)) {
      denom          <- abs(merged[[ly_col]])
      denom[denom == 0] <- NA                        # avoid /0
      out[[m]]       <- (merged[[su_col]] - merged[[ly_col]]) / denom * 100
    }
  }
  out
}

ortho_pre_03  <- compute_ortho_pct(master, "lying_inicial", "stand_up_inicial_03")
ortho_pre_05  <- compute_ortho_pct(master, "lying_inicial", "stand_up_inicial_05")
ortho_post_03 <- compute_ortho_pct(master, "lying_final",   "stand_up_final_03")
ortho_post_05 <- compute_ortho_pct(master, "lying_final",   "stand_up_final_05")

cat(sprintf(
  "Orthostatic %%change tables  pre_03 N=%d | pre_05 N=%d | post_03 N=%d | post_05 N=%d\n\n",
  nrow(ortho_pre_03), nrow(ortho_pre_05),
  nrow(ortho_post_03), nrow(ortho_post_05)))


# ── 3. ASSUMPTION CHECKS ──────────────────────────────────────────────────────
#
#   3a  Shapiro-Wilk on the paired difference (MR − VBT) per metric × slot
#       This is the relevant normality check for a paired test.
#
#   3b  Levene's test on raw values at fixed supine epochs (H1, H7).
# ─────────────────────────────────────────────────────────────────────────────

# 3a — Shapiro-Wilk ────────────────────────────────────────────────────────────

sw_diff <- function(df, hypothesis, metrics = key_metrics) {
  # pivot to wide so each row is one participant with MR and VBT values
  wide <- reshape(
    df[, c("participant_id", "protocol", metrics)],
    idvar     = "participant_id",
    timevar   = "protocol",
    direction = "wide"
  )

  do.call(rbind, lapply(metrics, function(m) {
    mc_mr  <- paste0(m, ".MR");  mc_vbt <- paste0(m, ".VBT")
    if (!all(c(mc_mr, mc_vbt) %in% names(wide))) return(NULL)
    d <- wide[[mc_mr]] - wide[[mc_vbt]]
    d <- d[!is.na(d)]
    if (length(d) < 3) return(NULL)
    sw <- tryCatch(shapiro.test(d), error = function(e) NULL)
    if (is.null(sw)) return(NULL)
    data.frame(
      hypothesis   = hypothesis,
      metric       = m,
      n_pairs      = length(d),
      mean_diff    = round(mean(d), 4),
      SD_diff      = round(sd(d),   4),
      SW_W         = round(sw$statistic, 4),
      SW_p         = round(sw$p.value,   4),
      normal       = sw$p.value > ALPHA,
      stringsAsFactors = FALSE
    )
  }))
}

# convenience: extract fixed epoch as a participant × protocol × metrics table
epoch_df <- function(ep) {
  master[master$epoch == ep, c("participant_id", "protocol", key_metrics)]
}

sw_results <- do.call(rbind, list(
  sw_diff(epoch_df("lying_inicial"),    "H1_lying_inicial"),
  sw_diff(ortho_pre_03,                 "H2_ortho_pre_03_pct"),
  sw_diff(ortho_pre_05,                 "H3_ortho_pre_05_pct"),
  sw_diff(epoch_df("lying_final"),      "H7_lying_final"),
  sw_diff(ortho_post_03,                "H8_ortho_post_03_pct"),
  sw_diff(ortho_post_05,                "H9_ortho_post_05_pct")
))

write.csv(sw_results,
          file.path(OUT_DIR, "results_00_shapiro_normality.csv"),
          row.names = FALSE)
cat("Saved: results_00_shapiro_normality.csv\n")


# 3b — Levene's test (fixed supine epochs H1 and H7) ──────────────────────────

levene_epoch <- function(ep, hypothesis, metrics = key_metrics) {
  sub <- master[master$epoch == ep, c("participant_id", "protocol", metrics)]
  do.call(rbind, lapply(metrics, function(m) {
    vals <- sub[[m]]; grp <- sub$protocol
    ok   <- !is.na(vals)
    if (sum(ok) < 4 || length(unique(grp[ok])) < 2) return(NULL)
    lev <- tryCatch(car::leveneTest(vals[ok] ~ grp[ok]), error = function(e) NULL)
    if (is.null(lev)) return(NULL)
    data.frame(
      hypothesis = hypothesis, metric = m,
      Levene_F   = round(lev[1, "F value"], 4),
      Levene_p   = round(lev[1, "Pr(>F)"],  4),
      equal_var  = lev[1, "Pr(>F)"] > ALPHA,
      stringsAsFactors = FALSE
    )
  }))
}

levene_results <- do.call(rbind, list(
  levene_epoch("lying_inicial", "H1_lying_inicial"),
  levene_epoch("lying_final",   "H7_lying_final")
))

write.csv(levene_results,
          file.path(OUT_DIR, "results_00_levene_variance.csv"),
          row.names = FALSE)
cat("Saved: results_00_levene_variance.csv\n\n")


# ── 4. HELPERS ────────────────────────────────────────────────────────────────

# 4a — Paired t-test or Wilcoxon signed-rank, per metric ──────────────────────
#
#  Decision per metric:
#    SW p > 0.05  →  paired t-test   + Cohen's dz
#    SW p ≤ 0.05  →  Wilcoxon SR     + rank-biserial r  (normal approximation)
#  BH-FDR applied across all metrics within each hypothesis call.
# ─────────────────────────────────────────────────────────────────────────────

paired_test_table <- function(df, hypothesis, sw_tbl,
                               metrics = key_metrics, apply_fdr = TRUE) {
  wide <- reshape(
    df[, c("participant_id", "protocol", metrics)],
    idvar     = "participant_id",
    timevar   = "protocol",
    direction = "wide"
  )

  results <- do.call(rbind, lapply(metrics, function(m) {
    mc_mr  <- paste0(m, ".MR");  mc_vbt <- paste0(m, ".VBT")
    if (!all(c(mc_mr, mc_vbt) %in% names(wide))) return(NULL)

    mr  <- wide[[mc_mr]];  vbt <- wide[[mc_vbt]]
    ok  <- !is.na(mr) & !is.na(vbt)
    mr  <- mr[ok];  vbt <- vbt[ok]
    if (length(mr) < 5) return(NULL)

    n        <- length(mr)
    diff_vec <- mr - vbt

    # normality lookup for this metric
    sw_row   <- sw_tbl[sw_tbl$metric == m, ]
    use_t    <- nrow(sw_row) > 0 && isTRUE(sw_row$normal[1])

    if (use_t) {
      tt  <- t.test(mr, vbt, paired = TRUE)
      dz  <- mean(diff_vec) / sd(diff_vec)
      data.frame(
        hypothesis   = hypothesis,
        metric       = m,
        n            = n,
        test         = "paired_t",
        statistic    = round(tt$statistic, 4),
        df_approx    = round(tt$parameter, 2),
        p_value      = round(tt$p.value,   5),
        effect_size  = round(dz, 3),
        effect_label = "Cohen_dz",
        median_MR    = round(median(mr),  4),
        median_VBT   = round(median(vbt), 4),
        mean_MR      = round(mean(mr),    4),
        mean_VBT     = round(mean(vbt),   4),
        stringsAsFactors = FALSE
      )
    } else {
      wt   <- wilcox.test(mr, vbt, paired = TRUE, exact = FALSE)
      W    <- as.numeric(wt$statistic)
      # rank-biserial r via normal approximation: r = Z / sqrt(N)
      mu_W <- n * (n + 1) / 4
      sd_W <- sqrt(n * (n + 1) * (2 * n + 1) / 24)
      Z    <- (W - mu_W) / sd_W
      r_rb <- Z / sqrt(n)
      data.frame(
        hypothesis   = hypothesis,
        metric       = m,
        n            = n,
        test         = "wilcoxon_SR",
        statistic    = round(W,    4),
        df_approx    = NA_real_,
        p_value      = round(wt$p.value, 5),
        effect_size  = round(r_rb, 3),
        effect_label = "rank_biserial_r",
        median_MR    = round(median(mr),  4),
        median_VBT   = round(median(vbt), 4),
        mean_MR      = round(mean(mr),    4),
        mean_VBT     = round(mean(vbt),   4),
        stringsAsFactors = FALSE
      )
    }
  }))

  if (is.null(results) || nrow(results) == 0) return(results)

  if (apply_fdr) {
    results$p_fdr   <- round(p.adjust(results$p_value, method = "BH"), 5)
    results$sig_raw <- results$p_value < ALPHA
    results$sig_fdr <- results$p_fdr   < ALPHA
  }
  results
}


# 4b — Lactate paired test (single variable, no FDR needed) ───────────────────

lac_paired_test <- function(hypothesis, mr_vals, vbt_vals, label) {
  ok  <- !is.na(mr_vals) & !is.na(vbt_vals)
  mr  <- mr_vals[ok];  vbt <- vbt_vals[ok]
  if (length(mr) < 5) {
    warning(sprintf("%s — %s: fewer than 5 complete pairs, skipping.", hypothesis, label))
    return(NULL)
  }
  n        <- length(mr)
  diff_vec <- mr - vbt
  sw       <- tryCatch(shapiro.test(diff_vec), error = function(e) list(p.value = 0))

  if (sw$p.value > ALPHA) {
    tt  <- t.test(mr, vbt, paired = TRUE)
    dz  <- mean(diff_vec) / sd(diff_vec)
    data.frame(
      hypothesis   = hypothesis, variable = label, n = n,
      SW_p_diff    = round(sw$p.value, 4),
      test         = "paired_t",
      statistic    = round(tt$statistic, 4),
      df_approx    = round(tt$parameter, 2),
      p_value      = round(tt$p.value,   5),
      effect_size  = round(dz, 3),
      effect_label = "Cohen_dz",
      median_MR    = round(median(mr),  4), median_VBT = round(median(vbt), 4),
      mean_MR      = round(mean(mr),    4), mean_VBT   = round(mean(vbt),   4),
      stringsAsFactors = FALSE
    )
  } else {
    wt   <- wilcox.test(mr, vbt, paired = TRUE, exact = FALSE)
    W    <- as.numeric(wt$statistic)
    mu_W <- n * (n + 1) / 4
    sd_W <- sqrt(n * (n + 1) * (2 * n + 1) / 24)
    r_rb <- (W - mu_W) / sd_W / sqrt(n)
    data.frame(
      hypothesis   = hypothesis, variable = label, n = n,
      SW_p_diff    = round(sw$p.value, 4),
      test         = "wilcoxon_SR",
      statistic    = round(W, 4),
      df_approx    = NA_real_,
      p_value      = round(wt$p.value,  5),
      effect_size  = round(r_rb, 3),
      effect_label = "rank_biserial_r",
      median_MR    = round(median(mr),  4), median_VBT = round(median(vbt), 4),
      mean_MR      = round(mean(mr),    4), mean_VBT   = round(mean(vbt),   4),
      stringsAsFactors = FALSE
    )
  }
}


# 4c — Linear Mixed Model (lmerTest) ──────────────────────────────────────────
#
#  Model:  metric ~ protocol * phase_pos + (1 | participant_id)
#    - protocol main effect: overall level difference MR vs VBT
#    - phase_pos main effect: common time trend
#    - interaction: protocol × phase_pos = difference in trajectory slope
#  Fit with ML (REML=FALSE) for LRT of the interaction term.
#  Satterthwaite p-values from lmerTest summary.
# ─────────────────────────────────────────────────────────────────────────────

run_lmm <- function(df_sub, hypothesis, metrics = key_metrics) {

  do.call(rbind, lapply(metrics, function(m) {
    sub <- df_sub[!is.na(df_sub[[m]]),
                  c("participant_id", "protocol", "phase_pos", m)]
    names(sub)[4] <- "value"

    if (nrow(sub) < 10 || length(unique(sub$participant_id)) < 5) return(NULL)

    # full model
    mod_full <- tryCatch(
      lmerTest::lmer(
        value ~ protocol * phase_pos + (1 | participant_id),
        data    = sub,
        REML    = FALSE,
        control = lmerControl(optimizer = "bobyqa")
      ),
      error = function(e) { message("  LMM error for ", m, ": ", e$message); NULL }
    )
    if (is.null(mod_full)) return(NULL)

    # reduced model (no interaction) for LRT
    mod_noint <- tryCatch(
      lmerTest::lmer(
        value ~ protocol + phase_pos + (1 | participant_id),
        data    = sub,
        REML    = FALSE,
        control = lmerControl(optimizer = "bobyqa")
      ),
      error = function(e) NULL
    )

    # coefficient table with Satterthwaite df and p-values
    coef_tbl           <- as.data.frame(coef(summary(mod_full)))
    coef_tbl$term      <- rownames(coef_tbl)
    rownames(coef_tbl) <- NULL
    names(coef_tbl)[names(coef_tbl) == "Estimate"]   <- "estimate"
    names(coef_tbl)[names(coef_tbl) == "Std. Error"] <- "SE"
    names(coef_tbl)[names(coef_tbl) == "df"]         <- "df_satterthwaite"
    names(coef_tbl)[names(coef_tbl) == "t value"]    <- "t_value"
    names(coef_tbl)[names(coef_tbl) == "Pr(>|t|)"]  <- "p_value"

    coef_tbl$hypothesis <- hypothesis
    coef_tbl$metric     <- m
    coef_tbl$singular   <- isSingular(mod_full)
    coef_tbl$n_obs      <- nrow(sub)
    coef_tbl$n_subjects <- length(unique(sub$participant_id))

    # LRT p-value for the interaction term
    if (!is.null(mod_noint)) {
      lrt <- tryCatch(anova(mod_noint, mod_full), error = function(e) NULL)
      if (!is.null(lrt) && nrow(lrt) >= 2) {
        coef_tbl$LRT_interaction_chisq <- round(lrt$Chisq[2],           4)
        coef_tbl$LRT_interaction_p     <- round(lrt$`Pr(>Chisq)`[2],    5)
      }
    }

    coef_tbl[, c("hypothesis", "metric", "term",
                 "estimate", "SE", "df_satterthwaite", "t_value", "p_value",
                 "LRT_interaction_chisq", "LRT_interaction_p",
                 "singular", "n_obs", "n_subjects")] |>
      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 5)))
  }))
}


# ── 5. PREPARE LACTATE WIDE TABLE ─────────────────────────────────────────────

lac_long <- master[master$epoch %in% c("lying_inicial", "lying_final") &
                     !is.na(master$lactate),
                   c("participant_id", "protocol", "epoch", "lactate")]

lac_wide <- reshape(lac_long,
                    idvar     = c("participant_id", "protocol"),
                    timevar   = "epoch",
                    direction = "wide")
names(lac_wide) <- gsub("lactate\\.", "", names(lac_wide))

# rename epoch-derived columns
if ("lying_inicial" %in% names(lac_wide)) {
  names(lac_wide)[names(lac_wide) == "lying_inicial"] <- "lac_base"
}
if ("lying_final" %in% names(lac_wide)) {
  names(lac_wide)[names(lac_wide) == "lying_final"] <- "lac_post"
}

lac_wide$lac_delta <- lac_wide$lac_post - lac_wide$lac_base

# split by protocol for paired helper
lac_MR  <- lac_wide[lac_wide$protocol == "MR",  ]
lac_VBT <- lac_wide[lac_wide$protocol == "VBT", ]
lac_both <- merge(lac_MR, lac_VBT,
                  by     = "participant_id",
                  suffixes = c("_MR", "_VBT"))

cat(sprintf("Lactate pairs available: N=%d\n\n", nrow(lac_both)))


# ── 6. RUN HYPOTHESES ─────────────────────────────────────────────────────────

# H1 — Pre-exercise resting ANS (lying_inicial) ────────────────────────────────
cat("─── H1: Pre-exercise ANS — lying_inicial ───────────────────────\n")
sw_H1  <- sw_results[sw_results$hypothesis == "H1_lying_inicial", ]
res_H1 <- paired_test_table(epoch_df("lying_inicial"), "H1", sw_H1)
print(res_H1[, c("metric","test","statistic","p_value","p_fdr","sig_fdr","effect_size")])
cat("\n")

# H2 — Short-term pre-exercise orthostatic reactivity (pre_03 % change) ────────
cat("─── H2: Pre-exercise ortho — stand_up_inicial_03 vs lying (%%chg) ─\n")
sw_H2  <- sw_results[sw_results$hypothesis == "H2_ortho_pre_03_pct", ]
res_H2 <- paired_test_table(ortho_pre_03, "H2", sw_H2)
print(res_H2[, c("metric","test","statistic","p_value","p_fdr","sig_fdr","effect_size")])
cat("\n")

# H3 — Medium-term pre-exercise orthostatic reactivity (pre_05 % change) ───────
cat("─── H3: Pre-exercise ortho — stand_up_inicial_05 vs lying (%%chg) ─\n")
sw_H3  <- sw_results[sw_results$hypothesis == "H3_ortho_pre_05_pct", ]
res_H3 <- paired_test_table(ortho_pre_05, "H3", sw_H3)
print(res_H3[, c("metric","test","statistic","p_value","p_fdr","sig_fdr","effect_size")])
cat("\n")

# H4 — Pre-exercise blood lactate ──────────────────────────────────────────────
cat("─── H4: Pre-exercise lactate — lying_inicial ────────────────────\n")
res_H4 <- lac_paired_test("H4",
                           lac_both$lac_base_MR,
                           lac_both$lac_base_VBT,
                           "lactate_baseline")
print(res_H4)
cat("\n")

# H5 — Exercise-phase ANS trajectory (LMM, trimmed) ───────────────────────────
cat("─── H5: Exercise trajectory — LMM (trimmed epochs) ─────────────\n")
exer_sub <- master[master$trimmed == TRUE & master$phase == "exercise", ]
cat(sprintf("  Exercise sub-data: %d rows  |  %d participants\n",
            nrow(exer_sub), length(unique(exer_sub$participant_id))))
res_H5 <- run_lmm(exer_sub, "H5_exercise_trajectory")
terms_of_interest <- c("protocolVBT", "phase_pos", "protocolVBT:phase_pos")
print(res_H5[res_H5$term %in% terms_of_interest,
             c("metric","term","estimate","p_value",
               "LRT_interaction_p","singular")])
cat("\n")

# H6 — Post-exercise ANS trajectory (LMM, trimmed) ────────────────────────────
cat("─── H6: Post-exercise trajectory — LMM (trimmed epochs) ────────\n")
post_sub <- master[master$trimmed == TRUE & grepl("^post", master$epoch), ]
cat(sprintf("  Post sub-data: %d rows  |  %d participants\n",
            nrow(post_sub), length(unique(post_sub$participant_id))))
res_H6 <- run_lmm(post_sub, "H6_post_trajectory")
print(res_H6[res_H6$term %in% terms_of_interest,
             c("metric","term","estimate","p_value",
               "LRT_interaction_p","singular")])
cat("\n")

# H7 — Post-exercise resting ANS (lying_final) ────────────────────────────────
cat("─── H7: Post-exercise ANS — lying_final ────────────────────────\n")
sw_H7  <- sw_results[sw_results$hypothesis == "H7_lying_final", ]
res_H7 <- paired_test_table(epoch_df("lying_final"), "H7", sw_H7)
print(res_H7[, c("metric","test","statistic","p_value","p_fdr","sig_fdr","effect_size")])
cat("\n")

# H8 — Short-term post-exercise orthostatic reactivity (post_03 % change) ──────
cat("─── H8: Post-exercise ortho — stand_up_final_03 vs lying (%%chg) ──\n")
sw_H8  <- sw_results[sw_results$hypothesis == "H8_ortho_post_03_pct", ]
res_H8 <- paired_test_table(ortho_post_03, "H8", sw_H8)
print(res_H8[, c("metric","test","statistic","p_value","p_fdr","sig_fdr","effect_size")])
cat("\n")

# H9 — Medium-term post-exercise orthostatic reactivity (post_05 % change) ─────
cat("─── H9: Post-exercise ortho — stand_up_final_05 vs lying (%%chg) ──\n")
sw_H9  <- sw_results[sw_results$hypothesis == "H9_ortho_post_05_pct", ]
res_H9 <- paired_test_table(ortho_post_05, "H9", sw_H9)
print(res_H9[, c("metric","test","statistic","p_value","p_fdr","sig_fdr","effect_size")])
cat("\n")

# H10 — Post-exercise and delta blood lactate ──────────────────────────────────
cat("─── H10: Post-exercise lactate — lying_final + delta ───────────\n")
res_H10 <- do.call(rbind, list(
  lac_paired_test("H10", lac_both$lac_post_MR,  lac_both$lac_post_VBT,  "lactate_post"),
  lac_paired_test("H10", lac_both$lac_delta_MR, lac_both$lac_delta_VBT, "delta_lactate")
))
print(res_H10)
cat("\n")


# ── 7. EXPORT RESULTS ─────────────────────────────────────────────────────────

# Fixed-epoch HRV tests — H1 and H7
write.csv(
  do.call(rbind, list(res_H1, res_H7)),
  file.path(OUT_DIR, "results_H1_H7_fixed_epochs_hrv.csv"),
  row.names = FALSE
)

# Orthostatic reactivity tests — H2, H3, H8, H9
write.csv(
  do.call(rbind, list(res_H2, res_H3, res_H8, res_H9)),
  file.path(OUT_DIR, "results_H2_H3_H8_H9_ortho_pct.csv"),
  row.names = FALSE
)

# Lactate tests — H4 and H10
write.csv(
  do.call(rbind, list(res_H4, res_H10)),
  file.path(OUT_DIR, "results_H4_H10_lactate.csv"),
  row.names = FALSE
)

# LMM trajectory — H5 (exercise)
write.csv(
  res_H5,
  file.path(OUT_DIR, "results_H5_exercise_LMM.csv"),
  row.names = FALSE
)

# LMM trajectory — H6 (post-exercise)
write.csv(
  res_H6,
  file.path(OUT_DIR, "results_H6_post_LMM.csv"),
  row.names = FALSE
)

cat("═══════════════════════════════════════════════════════════════\n")
cat("All outputs saved to:", OUT_DIR, "\n")
cat("  results_00_shapiro_normality.csv\n")
cat("  results_00_levene_variance.csv\n")
cat("  results_H1_H7_fixed_epochs_hrv.csv\n")
cat("  results_H2_H3_H8_H9_ortho_pct.csv\n")
cat("  results_H4_H10_lactate.csv\n")
cat("  results_H5_exercise_LMM.csv\n")
cat("  results_H6_post_LMM.csv\n")
cat("═══════════════════════════════════════════════════════════════\n")
