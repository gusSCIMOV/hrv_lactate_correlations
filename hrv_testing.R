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
cat("Section A outputs saved to:", OUT_DIR, "\n")
cat("  results_00_shapiro_normality.csv\n")
cat("  results_00_levene_variance.csv\n")
cat("  results_H1_H7_fixed_epochs_hrv.csv\n")
cat("  results_H2_H3_H8_H9_ortho_pct.csv\n")
cat("  results_H4_H10_lactate.csv\n")
cat("  results_H5_exercise_LMM.csv\n")
cat("  results_H6_post_LMM.csv\n")
cat("═══════════════════════════════════════════════════════════════\n")


# ═════════════════════════════════════════════════════════════════════════════
# SECTION B — H11 to H15: Δ HRV comparisons & Metabolic–Autonomic Coupling
# Appended to hrv_testing.R
#
# Requires: master, key_metrics, lac_wide, lac_both, OUT_DIR, ALPHA,
#           ortho_pre_03/05, ortho_post_03/05, sw_diff(), paired_test_table()
#           — all defined in Section A above.
# ═════════════════════════════════════════════════════════════════════════════


# B0. Additional packages ─────────────────────────────────────────────────────

for (.pkg in c("ppcor", "boot", "glmnet", "lmtest")) {
  if (!requireNamespace(.pkg, quietly = TRUE)) install.packages(.pkg)
}
suppressPackageStartupMessages({
  library(ppcor)
  library(boot)
  library(glmnet)
  library(lmtest)
})

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("Section B — H11 to H15: Δ HRV & Metabolic–Autonomic Coupling\n")
cat("═══════════════════════════════════════════════════════════════\n\n")


# B1. ΔHRV rest-state table  (lying_final − lying_inicial, per participant × protocol) ─

hrv_lying_base  <- master[master$epoch == "lying_inicial",
                           c("participant_id", "protocol", key_metrics)]
hrv_lying_final <- master[master$epoch == "lying_final",
                           c("participant_id", "protocol", key_metrics)]

delta_hrv_wide <- merge(hrv_lying_base, hrv_lying_final,
                        by       = c("participant_id", "protocol"),
                        suffixes = c("_base", "_final"))

for (.m in key_metrics) {
  delta_hrv_wide[[paste0("delta_", .m)]] <-
    delta_hrv_wide[[paste0(.m, "_final")]] - delta_hrv_wide[[paste0(.m, "_base")]]
}

delta_cols     <- paste0("delta_", key_metrics)
delta_hrv_slim <- delta_hrv_wide[, c("participant_id", "protocol", delta_cols)]


# B2. Δ orthostatic reactivity tables (ortho_post − ortho_pre per metric) ──────
#   ortho_pre_03 / ortho_pre_05 / ortho_post_03 / ortho_post_05 already exist
#   from Section 2 (compute_ortho_pct).
#   Output columns: dortho_{03|05}_{metric}

make_delta_ortho <- function(post_tbl, pre_tbl, suffix) {
  mg  <- merge(pre_tbl, post_tbl,
               by       = c("participant_id", "protocol"),
               suffixes = c("_pre", "_post"))
  out <- mg[, c("participant_id", "protocol")]
  for (.m in key_metrics) {
    pre_c  <- paste0(.m, "_pre");  post_c <- paste0(.m, "_post")
    if (all(c(pre_c, post_c) %in% names(mg)))
      out[[paste0("dortho_", suffix, "_", .m)]] <- mg[[post_c]] - mg[[pre_c]]
  }
  out
}

delta_ortho_03 <- make_delta_ortho(ortho_post_03, ortho_pre_03, "03")
delta_ortho_05 <- make_delta_ortho(ortho_post_05, ortho_pre_05, "05")

cat(sprintf(
  "Δ-tables  |  delta_hrv N=%d  |  dortho_03 N=%d  |  dortho_05 N=%d\n\n",
  nrow(delta_hrv_slim), nrow(delta_ortho_03), nrow(delta_ortho_05)))


# B3. Master coupling table (per participant × protocol) ──────────────────────
#     Joins: Δ HRV + Δ ortho_03 + Δ ortho_05 + lactate (lac_wide from Section 5)

coupling <- Reduce(
  function(a, b) merge(a, b, by = c("participant_id", "protocol"), all = TRUE),
  list(delta_hrv_slim, delta_ortho_03, delta_ortho_05, lac_wide)
)
coupling$protocol_bin <- as.integer(coupling$protocol == "VBT")   # MR=0, VBT=1

cat(sprintf("Coupling table: %d rows (participant × protocol)\n\n", nrow(coupling)))


# ─────────────────────────────────────────────────────────────────────────────
# H11 — ΔHRV rest-state (lying_final − lying_inicial)  MR vs VBT
# ─────────────────────────────────────────────────────────────────────────────
cat("─── H11: ΔHRV rest-state — MR vs VBT ──────────────────────────\n")

# Rename delta_* → bare metric names so sw_diff / paired_test_table can be reused
delta_hrv_test <- delta_hrv_slim
names(delta_hrv_test) <- gsub("^delta_", "", names(delta_hrv_test))

sw_H11  <- sw_diff(delta_hrv_test, "H11_delta_hrv")
res_H11 <- paired_test_table(delta_hrv_test, "H11", sw_H11)

write.csv(res_H11,
          file.path(OUT_DIR, "results_H11_delta_hrv.csv"),
          row.names = FALSE)
cat("Saved: results_H11_delta_hrv.csv\n")
print(res_H11[, c("metric", "test", "statistic", "p_value", "p_fdr", "sig_fdr", "effect_size")])
cat("\n")


# ─────────────────────────────────────────────────────────────────────────────
# H12 — Δ orthostatic reactivity (Δortho_03 and Δortho_05)  MR vs VBT
# ─────────────────────────────────────────────────────────────────────────────
cat("─── H12: Δ orthostatic reactivity — MR vs VBT ──────────────────\n")

delta_ortho_03_test <- delta_ortho_03
names(delta_ortho_03_test) <- gsub("^dortho_03_", "", names(delta_ortho_03_test))

delta_ortho_05_test <- delta_ortho_05
names(delta_ortho_05_test) <- gsub("^dortho_05_", "", names(delta_ortho_05_test))

sw_H12_03  <- sw_diff(delta_ortho_03_test, "H12_ortho_03")
sw_H12_05  <- sw_diff(delta_ortho_05_test, "H12_ortho_05")

res_H12_03 <- paired_test_table(delta_ortho_03_test, "H12_ortho_03", sw_H12_03)
res_H12_05 <- paired_test_table(delta_ortho_05_test, "H12_ortho_05", sw_H12_05)
res_H12    <- do.call(rbind, list(res_H12_03, res_H12_05))

write.csv(res_H12,
          file.path(OUT_DIR, "results_H12_delta_ortho.csv"),
          row.names = FALSE)
cat("Saved: results_H12_delta_ortho.csv\n")
print(res_H12[, c("hypothesis", "metric", "test", "p_value", "p_fdr", "sig_fdr", "effect_size")])
cat("\n")


# ─────────────────────────────────────────────────────────────────────────────
# Helpers for H13 / H14
# ─────────────────────────────────────────────────────────────────────────────

# spearman_boot: Spearman ρ + percentile bootstrap 95% CI (R iterations)
spearman_boot <- function(x, y, R = 2000) {
  ok <- !is.na(x) & !is.na(y)
  x  <- x[ok];  y <- y[ok];  n <- length(x)
  if (n < 5)
    return(list(rho = NA_real_, ci_low = NA_real_,
                ci_high = NA_real_, p = NA_real_, n = n))
  sp      <- cor.test(x, y, method = "spearman", exact = FALSE)
  rho     <- as.numeric(sp$estimate)
  p_val   <- sp$p.value
  boot_fn <- function(dat, idx) cor(dat[idx, 1], dat[idx, 2], method = "spearman")
  b       <- tryCatch(boot::boot(data = cbind(x, y), statistic = boot_fn, R = R),
                      error = function(e) NULL)
  ci      <- if (!is.null(b))
               tryCatch(boot::boot.ci(b, type = "perc", conf = 0.95),
                        error = function(e) NULL)
             else NULL
  list(rho     = rho,
       ci_low  = if (!is.null(ci)) ci$percent[4] else NA_real_,
       ci_high = if (!is.null(ci)) ci$percent[5] else NA_real_,
       p       = p_val,
       n       = n)
}

# fisher_z_test: two-tailed test of ρ1 == ρ2 via Fisher's Z transformation
fisher_z_test <- function(rho1, n1, rho2, n2) {
  if (any(is.na(c(rho1, n1, rho2, n2))) || n1 < 4 || n2 < 4)
    return(list(Z = NA_real_, p = NA_real_))
  z1 <- atanh(rho1);  z2 <- atanh(rho2)
  se <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  Z  <- (z1 - z2) / se
  list(Z = round(Z, 4), p = round(2 * pnorm(-abs(Z)), 5))
}

# run_spearman: per-protocol Spearman ρ + bootstrap CI + inter-protocol Fisher's Z
run_spearman <- function(coupling_tbl, lactate_var, hrv_vars, hypothesis, R = 2000) {
  result <- do.call(rbind, lapply(c("MR", "VBT"), function(prot) {
    sub <- coupling_tbl[coupling_tbl$protocol == prot, ]
    lac <- sub[[lactate_var]]
    do.call(rbind, lapply(hrv_vars, function(dv) {
      if (!dv %in% names(sub)) return(NULL)
      res <- spearman_boot(lac, sub[[dv]], R = R)
      data.frame(hypothesis    = hypothesis,
                 protocol      = prot,
                 lactate_var   = lactate_var,
                 hrv_delta_var = dv,
                 n             = res$n,
                 rho           = round(res$rho,     4),
                 rho_CI_low    = round(res$ci_low,  4),
                 rho_CI_high   = round(res$ci_high, 4),
                 p_value       = round(res$p,       5),
                 stringsAsFactors = FALSE)
    }))
  }))
  if (is.null(result) || nrow(result) == 0) return(result)

  # BH-FDR across all protocols × vars within the hypothesis
  result$p_fdr   <- round(p.adjust(result$p_value, method = "BH"), 5)
  result$sig_fdr <- result$p_fdr < ALPHA

  # Fisher's Z per hrv_delta_var (protocol comparison)
  result$fisher_Z_MR_vs_VBT <- NA_real_
  result$fisher_p_MR_vs_VBT <- NA_real_
  for (.dv in unique(result$hrv_delta_var)) {
    .mr  <- result[result$hrv_delta_var == .dv & result$protocol == "MR",  ]
    .vbt <- result[result$hrv_delta_var == .dv & result$protocol == "VBT", ]
    if (nrow(.mr) > 0 && nrow(.vbt) > 0) {
      .fz <- fisher_z_test(.mr$rho[1], .mr$n[1], .vbt$rho[1], .vbt$n[1])
      result$fisher_Z_MR_vs_VBT[result$hrv_delta_var == .dv] <- .fz$Z
      result$fisher_p_MR_vs_VBT[result$hrv_delta_var == .dv] <- .fz$p
    }
  }
  result
}

# All Δ HRV predictor variable names that are present in the coupling table
hrv_corr_vars <- c(delta_cols,
                   paste0("dortho_03_", key_metrics),
                   paste0("dortho_05_", key_metrics))
hrv_corr_vars <- hrv_corr_vars[hrv_corr_vars %in% names(coupling)]


# ─────────────────────────────────────────────────────────────────────────────
# H13 — lactate_post × Δ HRV  (Spearman ρ, bootstrap CI, Fisher's Z)
# ─────────────────────────────────────────────────────────────────────────────
cat("─── H13: lactate_post × Δ HRV + Δ ortho — Spearman ───────────\n")

res_H13 <- run_spearman(coupling, "lac_post", hrv_corr_vars, "H13", R = 2000)

write.csv(res_H13,
          file.path(OUT_DIR, "results_H13_lactate_post_spearman.csv"),
          row.names = FALSE)
cat("Saved: results_H13_lactate_post_spearman.csv\n")
if (!is.null(res_H13) && nrow(res_H13) > 0) {
  .noteworthy <- res_H13[!is.na(res_H13$p_value) & res_H13$p_value < 0.10, ]
  if (nrow(.noteworthy) > 0)
    print(.noteworthy[, c("protocol", "hrv_delta_var", "rho", "p_value", "p_fdr", "sig_fdr")])
  else cat("  No associations with p < 0.10\n")
}
cat("\n")


# ─────────────────────────────────────────────────────────────────────────────
# H14 — delta_lactate × Δ HRV  (Spearman + partial Spearman, lac_base as covariate)
# ─────────────────────────────────────────────────────────────────────────────
cat("─── H14: delta_lactate × Δ HRV + Δ ortho — Spearman + partial ─\n")

res_H14_biv <- run_spearman(coupling, "lac_delta", hrv_corr_vars, "H14", R = 2000)

# Partial Spearman: control for lac_base (baseline lactate) within each protocol
partial_H14 <- do.call(rbind, lapply(c("MR", "VBT"), function(prot) {
  sub <- coupling[coupling$protocol == prot, ]
  do.call(rbind, lapply(hrv_corr_vars, function(dv) {
    ok  <- !is.na(sub$lac_delta) & !is.na(sub$lac_base) & !is.na(sub[[dv]])
    if (sum(ok) < 5) return(NULL)
    # ppcor::pcor on a 3-column matrix: [lac_delta, hrv, lac_base]
    # element [1,2] gives partial corr of lac_delta ~ hrv | lac_base
    mat <- cbind(sub$lac_delta[ok], sub[[dv]][ok], sub$lac_base[ok])
    pc  <- tryCatch(ppcor::pcor(mat, method = "spearman"), error = function(e) NULL)
    if (is.null(pc)) return(NULL)
    data.frame(protocol      = prot,
               hrv_delta_var = dv,
               partial_rho   = round(pc$estimate[1, 2], 4),
               partial_p     = round(pc$p.value [1, 2], 5),
               n_partial     = sum(ok),
               stringsAsFactors = FALSE)
  }))
}))

res_H14 <- merge(res_H14_biv, partial_H14,
                 by = c("protocol", "hrv_delta_var"), all.x = TRUE)

write.csv(res_H14,
          file.path(OUT_DIR, "results_H14_delta_lactate_spearman.csv"),
          row.names = FALSE)
cat("Saved: results_H14_delta_lactate_spearman.csv\n")
if (!is.null(res_H14) && nrow(res_H14) > 0) {
  .noteworthy <- res_H14[!is.na(res_H14$p_value) & res_H14$p_value < 0.10, ]
  if (nrow(.noteworthy) > 0) {
    .cols <- intersect(c("protocol", "hrv_delta_var", "rho", "p_value",
                         "p_fdr", "partial_rho", "partial_p"), names(.noteworthy))
    print(.noteworthy[, .cols])
  } else cat("  No associations with p < 0.10\n")
}
cat("\n")


# ─────────────────────────────────────────────────────────────────────────────
# H15 — LASSO: lactate_post ~ protocol_bin + Δ HRV predictors
#   Outcome : lac_post (mmol/L)
#   Predictors: protocol_bin (0=MR, 1=VBT) + delta_* + dortho_03_* + dortho_05_*
#   All continuous predictors z-scored before fitting.
#   Lambda selected by LOO-CV (nfolds = N); lambda.1se preferred.
#   OLS diagnostics first: VIF (car::vif), Cook's D, SW residuals, Breusch-Pagan.
#   Bootstrap stability: 500 iterations, select threshold ≥ 50%.
# ─────────────────────────────────────────────────────────────────────────────
cat("─── H15: LASSO — lactate_post as outcome ───────────────────────\n")

predictor_cols_h15 <- c("protocol_bin", hrv_corr_vars)

# Keep rows with non-missing outcome; filter predictors with ≥ 50% data coverage
lasso_df   <- coupling[!is.na(coupling$lac_post),
                        c("participant_id", "protocol", "lac_post", predictor_cols_h15)]
pred_cov   <- sapply(predictor_cols_h15, function(col) mean(!is.na(lasso_df[[col]])))
good_preds <- predictor_cols_h15[pred_cov >= 0.5]
lasso_cc   <- lasso_df[complete.cases(lasso_df[, c("lac_post", good_preds)]), ]

N_lasso <- nrow(lasso_cc)
P_lasso <- length(good_preds)
cat(sprintf("  Dataset: N=%d  |  predictors=%d (≥50%% coverage)\n", N_lasso, P_lasso))

if (N_lasso < 10) {
  cat("  LASSO skipped: insufficient complete cases (N < 10).\n\n")
} else {

  # Z-score continuous predictors; leave protocol_bin as 0/1
  X_mat <- as.matrix(lasso_cc[, good_preds])
  for (.col in setdiff(good_preds, "protocol_bin")) {
    .mu <- mean(X_mat[, .col], na.rm = TRUE)
    .s  <- sd(X_mat[, .col],   na.rm = TRUE)
    if (!is.na(.s) && .s > 0) X_mat[, .col] <- (X_mat[, .col] - .mu) / .s
  }
  y_lasso <- lasso_cc$lac_post

  # ── OLS diagnostics (VIF, Cook's D, SW residuals, Breusch-Pagan) ──────────
  #    Requires N > p + 2; skipped otherwise (LASSO still runs)
  if (N_lasso > P_lasso + 2) {
    ols_data    <- as.data.frame(X_mat)
    ols_data$y  <- y_lasso
    ols_formula <- as.formula(paste("y ~", paste(good_preds, collapse = " + ")))
    ols_fit     <- tryCatch(lm(ols_formula, data = ols_data), error = function(e) NULL)
  } else {
    ols_fit <- NULL
    cat(sprintf("  OLS diagnostics skipped (N=%d ≤ p+2=%d); LASSO proceeds.\n",
                N_lasso, P_lasso + 2))
  }

  infl_obs <- NULL

  if (!is.null(ols_fit)) {
    # VIF — multicollinearity check
    vif_raw <- tryCatch(car::vif(ols_fit), error = function(e) NULL)
    if (!is.null(vif_raw)) {
      diag_vif <- data.frame(predictor = names(vif_raw),
                              VIF       = round(as.numeric(vif_raw), 3),
                              VIF_flag  = as.numeric(vif_raw) > 10,
                              stringsAsFactors = FALSE)
      write.csv(diag_vif,
                file.path(OUT_DIR, "results_H15_lasso_diagnostics.csv"),
                row.names = FALSE)
      cat("  Saved: results_H15_lasso_diagnostics.csv\n")
      .flagged <- diag_vif$predictor[diag_vif$VIF_flag]
      if (length(.flagged) > 0)
        cat("  VIF > 10 flagged:", paste(.flagged, collapse = ", "), "\n")
      else
        cat("  VIF: no predictors exceed threshold (>10)\n")
    }

    # Cook's distance — influential observation detection
    .cooks     <- cooks.distance(ols_fit)
    .infl_flag <- .cooks > 4 / N_lasso
    infl_obs   <- data.frame(
      participant_id       = lasso_cc$participant_id,
      protocol             = lasso_cc$protocol,
      cooks_d              = round(.cooks, 5),
      influential          = .infl_flag,
      sensitivity_excluded = NA_character_,
      stringsAsFactors     = FALSE
    )
    cat(sprintf("  Influential obs (Cook's D > 4/N=%.3f): %d\n",
                4 / N_lasso, sum(.infl_flag)))

    # Shapiro-Wilk on OLS residuals
    .sw_res <- tryCatch(shapiro.test(residuals(ols_fit)), error = function(e) NULL)
    if (!is.null(.sw_res))
      cat(sprintf("  OLS residual SW: W=%.4f, p=%.4f  → %s\n",
                  as.numeric(.sw_res$statistic), .sw_res$p.value,
                  ifelse(.sw_res$p.value > 0.05,
                         "residual normality OK",
                         "violation — interpret LASSO; report as caveat")))

    # Breusch-Pagan — heteroscedasticity
    .bp <- tryCatch(lmtest::bptest(ols_fit), error = function(e) NULL)
    if (!is.null(.bp))
      cat(sprintf("  Breusch-Pagan: BP=%.4f, p=%.4f  → %s\n",
                  as.numeric(.bp$statistic), .bp$p.value,
                  ifelse(.bp$p.value > 0.05,
                         "homoscedastic OK",
                         "heteroscedastic — use robust SE on OLS interpretation")))
  }

  # ── LASSO cv.glmnet with LOO-CV ───────────────────────────────────────────
  set.seed(2025)
  cv_fit <- tryCatch(
    glmnet::cv.glmnet(X_mat, y_lasso, alpha = 1, nfolds = N_lasso),
    error = function(e) { message("  cv.glmnet error: ", e$message); NULL }
  )

  if (!is.null(cv_fit)) {
    lam_1se  <- cv_fit$lambda.1se
    lam_min  <- cv_fit$lambda.min

    coef_1se <- as.matrix(coef(cv_fit, s = "lambda.1se"))
    coef_min <- as.matrix(coef(cv_fit, s = "lambda.min"))

    coef_df <- data.frame(
      predictor  = rownames(coef_1se),
      lambda_min = round(as.numeric(coef_min), 5),
      lambda_1se = round(as.numeric(coef_1se), 5),
      selected   = as.numeric(coef_1se) != 0,
      stringsAsFactors = FALSE
    )

    # Cross-validation performance table
    cv_tbl <- data.frame(
      lambda        = round(cv_fit$lambda, 6),
      cvm           = round(cv_fit$cvm,    5),
      cvsd          = round(cv_fit$cvsd,   5),
      lambda_is_min = cv_fit$lambda == lam_min,
      lambda_is_1se = cv_fit$lambda == lam_1se,
      stringsAsFactors = FALSE
    )
    .cvm_1se              <- cv_fit$cvm[which.min(abs(cv_fit$lambda - lam_1se))]
    .cv_R2                <- round(1 - .cvm_1se / var(y_lasso), 4)
    cv_tbl$cv_R2_lambda_1se <- ifelse(cv_tbl$lambda_is_1se, .cv_R2, NA_real_)

    .sel_preds <- coef_df$predictor[coef_df$selected & coef_df$predictor != "(Intercept)"]
    cat(sprintf("  lambda.1se=%.4f  |  cv-R²=%.4f\n", lam_1se, .cv_R2))
    cat(sprintf("  Selected at lambda.1se (%d): %s\n",
                length(.sel_preds),
                if (length(.sel_preds) > 0) paste(.sel_preds, collapse = ", ") else "none"))

    # ── Bootstrap stability (500 iterations, threshold ≥ 50%) ─────────────
    cat("  Bootstrap stability (500 iter) ...\n")
    n_coef       <- nrow(coef_1se)
    boot_sel_mat <- matrix(0L, nrow = n_coef, ncol = 500,
                           dimnames = list(rownames(coef_1se), NULL))
    set.seed(42)
    for (.b in seq_len(500)) {
      .idx <- sample(N_lasso, N_lasso, replace = TRUE)
      .fb  <- tryCatch(
        glmnet::glmnet(X_mat[.idx, , drop = FALSE], y_lasso[.idx],
                       alpha = 1, lambda = lam_1se),
        error = function(e) NULL
      )
      if (!is.null(.fb))
        boot_sel_mat[, .b] <- as.integer(as.numeric(coef(.fb, s = lam_1se)) != 0)
    }
    .boot_pct                  <- rowMeans(boot_sel_mat) * 100
    coef_df$boot_selection_pct <- round(.boot_pct[coef_df$predictor], 1)
    coef_df$stable             <- coef_df$boot_selection_pct >= 50
    .stable_preds <- coef_df$predictor[coef_df$stable & coef_df$predictor != "(Intercept)"]
    cat(sprintf("  Stable predictors (≥50%% bootstrap): %s\n",
                if (length(.stable_preds) > 0) paste(.stable_preds, collapse = ", ")
                else "none"))

    write.csv(coef_df, file.path(OUT_DIR, "results_H15_lasso_coefficients.csv"),
              row.names = FALSE)
    write.csv(cv_tbl,  file.path(OUT_DIR, "results_H15_lasso_cv.csv"),
              row.names = FALSE)
    cat("  Saved: results_H15_lasso_coefficients.csv\n")
    cat("  Saved: results_H15_lasso_cv.csv\n")

    # ── Sensitivity: refit LASSO excluding each influential observation ────
    if (!is.null(infl_obs) && any(infl_obs$influential)) {
      for (.i in which(infl_obs$influential)) {
        .fb_s <- tryCatch(
          glmnet::glmnet(X_mat[-.i, , drop = FALSE], y_lasso[-.i],
                         alpha = 1, lambda = lam_1se),
          error = function(e) NULL
        )
        if (!is.null(.fb_s)) {
          .c_s <- as.numeric(coef(.fb_s, s = lam_1se))
          .p_s <- rownames(coef_1se)[.c_s != 0 & rownames(coef_1se) != "(Intercept)"]
          infl_obs$sensitivity_excluded[.i] <-
            ifelse(!setequal(.p_s, .sel_preds), "selection_changed", "stable")
        }
      }
      write.csv(infl_obs, file.path(OUT_DIR, "results_H15_influential_obs.csv"),
                row.names = FALSE)
      cat("  Saved: results_H15_influential_obs.csv\n")
    }

  }  # end if cv_fit
}  # end if N_lasso >= 10


# ─────────────────────────────────────────────────────────────────────────────
# Section B — final summary
# ─────────────────────────────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("Section B complete. Outputs written to:", OUT_DIR, "\n")
cat("  results_H11_delta_hrv.csv\n")
cat("  results_H12_delta_ortho.csv\n")
cat("  results_H13_lactate_post_spearman.csv\n")
cat("  results_H14_delta_lactate_spearman.csv\n")
cat("  results_H15_lasso_coefficients.csv      (if N ≥ 10)\n")
cat("  results_H15_lasso_cv.csv                (if N ≥ 10)\n")
cat("  results_H15_lasso_diagnostics.csv       (if N > p+2)\n")
cat("  results_H15_influential_obs.csv         (if OLS feasible + influential obs found)\n")
cat("═══════════════════════════════════════════════════════════════\n")

