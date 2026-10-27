# Supplementary sensitivity analysis: does including sub-0011 (flagged
# elsewhere in this pipeline as a whole-brain total-P outlier) change Fig.
# 1's primary current dipole moment (I_dip) and primary dipole power (P)
# results? Main-text Results exclude sub-0011 (n=30); this script quantifies
# the effect of including it (n=31) as a sensitivity check.
#
# Panels:
#   A) Histogram of whole-brain total P across all n=31 subjects, sub-0011
#      highlighted.
#   B) Primary I_dip, parcellated (Schaefer-600), averaged across all n=31
#      subjects (WITH sub-0011) -- compare against the main pipeline's n=30
#      (excluding sub-0011) version, e.g. rest_s600_mean_i_out_filt.png.
#   C) Primary dipole P, parcellated, averaged across all n=31 subjects --
#      compare against rms_2min_total_p_out_filt.png.
#
# Also prints: the >5 SD outlier statistic, the n=31-vs-n=30 mean P
# comparison, and which regions show the largest sub-0011-driven excess
# (with a check for whether they're disproportionately frontal).
#
# Inputs (produced by parcellate_central_surface_s600.m; per-subject files
# are NOT deleted when a subject is dropped from that script's SubjectNames,
# so all 31 subjects' files -- including sub-0011's -- are still on disk
# even though the main-pipeline subavg CSVs are now n=30):
# - s600_sub-<ID>_snr3_p_pri_2min_rest.csv (per-subject, per-region P)
# - s600_sub-<ID>_snr3_i_pri_2min_rest.csv (per-subject, per-region I_dip)

library(ggsegSchaefer)
library(ggseg)
library(dplyr)
library(ggplot2)
library(patchwork)

schaefer7_600 <- schaefer7_600()

outputs_dir <- "/export02/data/vikramn/hbm_manuscript_code/outputs"
primary_p_dir <- file.path(outputs_dir, "primary_p", "cent_surf_fem_fwd")
fig_dir <- file.path(outputs_dir, "primary_p", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

OUTLIER_SUBJECT <- 'sub-0011'
SNR <- 3

# Same region-name helpers as schaeffer_rest_power_map.R / snr_sensitivity_supp_fig.R.
transform_s600_colnames <- function(column_names) {
  transformed_names <- vector("character", length(column_names))
  for (i in 1:length(column_names)) {
    original_name <- column_names[i]
    hemi_letter <- substr(original_name, nchar(original_name), nchar(original_name))
    hemi_code <- if (hemi_letter == 'L') 'LH' else 'RH'
    base_name <- substr(original_name, 1, nchar(original_name) - 2)
    transformed_names[i] <- paste0('7Networks_', hemi_code, '_', base_name)
  }
  return(transformed_names)
}

fix_region_names <- function(region) {
  region <- trimws(region)
  region <- sub("PFCI_", "PFCl_", region, fixed = TRUE)
  region <- sub("^7Networks_LH_Default_PFCdPFCm_8$", "7Networks_LH_Default_PFC_8", region)
  region <- sub("^7Networks_RH_Vis_45$", "7Networks_RH_Vis_44", region)
  region
}

# Same per-region outlier clipping used throughout schaeffer_rest_power_map.R
# (its "*_out_filt" columns/figures) -- Tukey's IQR rule (values beyond
# Q3+1.5*IQR clipped to the max non-outlier value), NOT a fixed SD cutoff.
# Applied here to Panels B/C below so a single sub-0011-driven region
# doesn't just saturate the whole color scale and wash out every other
# region -- matches how the main pipeline already handles this, rather than
# a one-off threshold invented just for this supplementary figure.
replace_outliers_with_max_non_outlier <- function(vector) {
  q <- quantile(vector, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  lower_bound <- q[1] - 1.5 * iqr
  upper_bound <- q[2] + 1.5 * iqr

  non_outliers <- vector[!is.na(vector) & vector >= lower_bound & vector <= upper_bound]
  if (length(non_outliers) == 0) return(vector)
  max_non_outlier <- max(non_outliers, na.rm = TRUE)

  vector[vector > upper_bound] <- max_non_outlier
  return(vector)
}

atlas_names <- unique(ggseg:::prepare_polygon_atlas(schaefer7_600, hemi = NULL, view = NULL,
                                                      position = position_brain(), context = TRUE,
                                                      focus = NULL)$region)

# =============================================================================
# Load all 31 subjects' per-region P and I_dip (SnrFixed=3)
# =============================================================================
load_long <- function(pattern) {
  files <- list.files(primary_p_dir, pattern = pattern, full.names = TRUE)
  rows <- list()
  for (f in files) {
    subject <- sub(".*(sub-[0-9]+)_snr.*", "\\1", basename(f))
    d <- read.csv(f, stringsAsFactors = FALSE)
    d$subject <- subject
    rows[[length(rows) + 1]] <- d
  }
  do.call(rbind, rows)
}

p_long <- load_long(sprintf('^s600_sub-[0-9]+_snr%d_p_pri_2min_rest\\.csv$', SNR))
i_long <- load_long(sprintf('^s600_sub-[0-9]+_snr%d_i_pri_2min_rest\\.csv$', SNR))

subjects_found <- sort(unique(p_long$subject))
cat(sprintf("Found %d subjects: %s\n", length(subjects_found), paste(subjects_found, collapse = ", ")))
if (!(OUTLIER_SUBJECT %in% subjects_found)) {
  stop(OUTLIER_SUBJECT, " not found in ", primary_p_dir,
       " -- its per-subject files may have been deleted; this sensitivity analysis needs them.")
}

# =============================================================================
# Outlier statistics: whole-brain total P per subject, sub-0011 vs. the
# remaining n=30 (mean/SD computed EXCLUDING sub-0011, per the >5 SD test
# definition -- an outlier can't be allowed to inflate its own threshold).
# =============================================================================
subject_totals_p <- p_long %>% group_by(subject) %>% summarise(total_p = sum(value, na.rm = TRUE), .groups = "drop")
write.csv(subject_totals_p, file.path(outputs_dir, "primary_p", "outlier_sensitivity_subject_totals_p.csv"), row.names = FALSE)

outlier_val <- subject_totals_p$total_p[subject_totals_p$subject == OUTLIER_SUBJECT]
rest_vals <- subject_totals_p$total_p[subject_totals_p$subject != OUTLIER_SUBJECT]
mean_excl <- mean(rest_vals)
sd_excl <- sd(rest_vals)
z_score <- (outlier_val - mean_excl) / sd_excl
mean_incl <- mean(subject_totals_p$total_p)
pct_increase <- 100 * (mean_incl - mean_excl) / mean_excl

cat(sprintf("\n=== Outlier test: %s whole-brain total P vs. remaining n=%d ===\n", OUTLIER_SUBJECT, length(rest_vals)))
cat(sprintf("  %s total P        = %.4e W\n", OUTLIER_SUBJECT, outlier_val))
cat(sprintf("  mean (excl., n=%d) = %.4e W\n", length(rest_vals), mean_excl))
cat(sprintf("  SD (excl., n=%d)   = %.4e W\n", length(rest_vals), sd_excl))
cat(sprintf("  z-score            = %.2f SD %s\n", abs(z_score), if (z_score > 0) "above" else "below"))
cat(sprintf("  exceeds 5 SD test? %s\n", if (abs(z_score) > 5) "YES" else "no"))
cat(sprintf("\n=== Main-text summary: including vs. excluding %s ===\n", OUTLIER_SUBJECT))
cat(sprintf("  mean total P, n=31 (including) = %.4e W\n", mean_incl))
cat(sprintf("  mean total P, n=30 (excluding) = %.4e W\n", mean_excl))
cat(sprintf("  including %s inflates the mean by %.1f%%\n", OUTLIER_SUBJECT, pct_increase))

# =============================================================================
# Panel A: histogram of whole-brain total P (n=31), sub-0011 highlighted.
# Linear binning is uninformative here -- with one point 16 SD out, a linear
# scale either squashes the n=30 cluster into 1-2 bars (few wide bins) or
# produces unreadable, arbitrarily-spaced tick labels (many narrow bins).
# Binning and plotting in log10 space instead spreads the cluster out
# properly while still accommodating the outlier, with standard "1-2-5"
# round-number tick labels rather than arbitrary bin-edge values.
# =============================================================================
subject_totals_p$is_outlier <- subject_totals_p$subject == OUTLIER_SUBJECT
subject_totals_p$log_p <- log10(subject_totals_p$total_p)
log_binwidth <- diff(range(log10(rest_vals))) / 12
threshold_5sd <- mean_excl + 5 * sd_excl

# Standard log-scale "1-2-5" ticks spanning the data range, in W.
log_breaks_125 <- function(lo, hi) {
  exps <- floor(log10(lo)):ceiling(log10(hi))
  cand <- as.vector(outer(c(1, 2, 5), 10^exps))
  sort(cand[cand >= lo / 1.5 & cand <= hi * 1.5])
}
x_breaks_w <- log_breaks_125(min(subject_totals_p$total_p), max(subject_totals_p$total_p))
x_breaks_log <- log10(x_breaks_w)
x_labels <- sub("e\\+?0*(-?[0-9]+)$", "e\\1", formatC(x_breaks_w, format = "e", digits = 0))

hist_panel <- ggplot(subject_totals_p, aes(x = log_p, fill = is_outlier)) +
  geom_histogram(binwidth = log_binwidth, color = "white", boundary = 0) +
  geom_rug(aes(color = is_outlier), sides = "b", linewidth = 0.7, length = unit(0.04, "npc"), show.legend = FALSE) +
  geom_vline(xintercept = log10(mean_excl), linetype = "dashed", color = "gray30", linewidth = 0.6) +
  geom_vline(xintercept = log10(threshold_5sd), linetype = "dashed", color = "#C44E52", linewidth = 0.6) +
  stat_bin(binwidth = log_binwidth, boundary = 0, geom = "text",
           aes(label = ifelse(after_stat(count) > 0, after_stat(count), "")),
           vjust = -0.4, size = 3.5, color = "gray30") +
  annotate("text", x = log10(mean_excl), y = Inf, label = "mean (n=30)", color = "gray30",
           angle = 90, hjust = 1.1, vjust = -0.5, size = 3.8) +
  annotate("text", x = log10(threshold_5sd), y = Inf, label = "5 SD threshold", color = "#C44E52",
           angle = 90, hjust = 1.1, vjust = -0.5, size = 3.8) +
  annotate("text", x = log10(outlier_val), y = 1, label = OUTLIER_SUBJECT, color = "#C44E52",
           angle = 90, hjust = -0.2, vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c(`FALSE` = "#4C72B0", `TRUE` = "#C44E52"), guide = "none") +
  scale_color_manual(values = c(`FALSE` = "#4C72B0", `TRUE` = "#C44E52")) +
  scale_x_continuous(breaks = x_breaks_log, labels = x_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Whole-brain total primary dipole P (W, log scale)", y = "Number of subjects",
       title = sprintf("A  Total primary dipole P across all subjects (n=%d, SnrFixed=%d)", nrow(subject_totals_p), SNR)) +
  theme_minimal(base_size = 14)
ggsave(file.path(fig_dir, "outlier_sensitivity_p_histogram.png"), plot = hist_panel, width = 9, height = 6, bg = "white")

# =============================================================================
# Panels B/C: parcellated I_dip and P, averaged across ALL n=31 subjects
# (with sub-0011), for visual comparison against the main n=30 figures.
# =============================================================================
mean_by_region <- function(long_df, value_col) {
  d <- long_df %>% group_by(region) %>% summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")
  d$region <- transform_s600_colnames(d$region)
  d$region <- fix_region_names(d$region)
  colnames(d)[colnames(d) == "mean_val"] <- value_col
  d
}

i_mean31 <- mean_by_region(i_long, "mean_i")
p_mean31 <- mean_by_region(p_long, "mean_p")
i_mean31$mean_i_out_filt <- replace_outliers_with_max_non_outlier(i_mean31$mean_i)
p_mean31$mean_p_out_filt <- replace_outliers_with_max_non_outlier(p_mean31$mean_p)

# Same colorbar limits as Fig. 1B's n=30 (excl. sub-0011) maps
# (rms_2min_mean_i_out_filt.png / rms_2min_total_p_out_filt.png in
# schaeffer_rest_power_map.R) -- those figures never set explicit
# scale_fill_distiller limits, so their color range is implicitly their own
# data's min/max; reproduced here from the same underlying broad CSV and
# the same thresholding function so both figures are on an identical scale.
# The point: a region that's MORE saturated here than it ever gets in Fig.
# 1 is a region that visibly "gained" power/current from including
# sub-0011 -- not just an artifact of two independently-autoscaled legends.
fig1_broad <- read.csv(file.path(primary_p_dir, "rms_2min_cent_pri_power_met_broad.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
fig1_p_out_filt <- replace_outliers_with_max_non_outlier(as.numeric(fig1_broad$mean_p))
fig1_i_out_filt <- replace_outliers_with_max_non_outlier(as.numeric(fig1_broad$mean_i))
fig1_p_limits <- range(fig1_p_out_filt, na.rm = TRUE)
fig1_i_limits <- range(fig1_i_out_filt, na.rm = TRUE)
cat(sprintf("\nUsing Fig. 1 (n=30) colorbar limits: mean_p_out_filt = [%.3e, %.3e] W, mean_i_out_filt = [%.3e, %.3e] A.m\n",
            fig1_p_limits[1], fig1_p_limits[2], fig1_i_limits[1], fig1_i_limits[2]))

brain_row <- function(data, varname, fill_label, title, limits) {
  ggplot(data) +
    geom_brain(data = data, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping = aes(fill = !!sym(varname))) +
    labs(fill = fill_label) + ggtitle(title) + theme_void() +
    scale_fill_distiller(palette = "Blues", direction = 1, limits = limits, oob = scales::squish)
}

fig_i31 <- brain_row(i_mean31, "mean_i_out_filt", "A.m",
                      "B  Primary I_dip, n=31, sub-0011 included (Thresholded, Fig. 1 color scale)", fig1_i_limits)
ggsave(file.path(fig_dir, "outlier_sensitivity_i_pri_n31.png"), plot = fig_i31, width = 10, height = 3.5, bg = "white")

fig_p31 <- brain_row(p_mean31, "mean_p_out_filt", "W",
                      "C  Primary dipole P, n=31, sub-0011 included (Thresholded, Fig. 1 color scale)", fig1_p_limits)
ggsave(file.path(fig_dir, "outlier_sensitivity_p_pri_n31.png"), plot = fig_p31, width = 10, height = 3.5, bg = "white")

combined <- hist_panel / fig_i31 / fig_p31 +
  plot_layout(heights = c(1.2, 1, 1)) +
  plot_annotation(title = sprintf("Sensitivity analysis: effect of including %s (outlier) on primary dipole results", OUTLIER_SUBJECT))
ggsave(file.path(fig_dir, "supp_fig4_outlier_sensitivity.png"), plot = combined, width = 10, height = 10, bg = "white")

# =============================================================================
# Which regions show the largest sub-0011-driven excess? (sub-0011's value
# minus the n=30-excluding-outlier mean, per region), and are they
# disproportionately frontal? Frontal markers based on the actual Schaefer-600
# 7-network region abbreviations observed in this data: PFC* (all prefrontal
# variants), OFC (orbitofrontal), FEF (frontal eye fields), FrOper (frontal
# operculum) -- verified against this atlas's real region names, not guessed.
# =============================================================================
p_wide30 <- p_long %>% filter(subject != OUTLIER_SUBJECT) %>%
  group_by(region) %>% summarise(mean_p_excl = mean(value, na.rm = TRUE), .groups = "drop")
p_sub11 <- p_long %>% filter(subject == OUTLIER_SUBJECT) %>% select(region, sub11_p = value)

region_excess <- merge(p_wide30, p_sub11, by = "region")
region_excess$excess <- region_excess$sub11_p - region_excess$mean_p_excl
region_excess$ratio <- region_excess$sub11_p / region_excess$mean_p_excl
region_excess <- region_excess[order(-region_excess$excess), ]

FRONTAL_PATTERN <- "PFC|OFC|FEF|FrOper"
region_excess$is_frontal <- grepl(FRONTAL_PATTERN, region_excess$region)

top_n <- 30 # top ~5% of 600 regions
top_regions <- head(region_excess, top_n)
pct_frontal_top <- 100 * mean(top_regions$is_frontal)
pct_frontal_all <- 100 * mean(region_excess$is_frontal)

write.csv(region_excess, file.path(outputs_dir, "primary_p", "outlier_sensitivity_region_excess.csv"), row.names = FALSE)

cat(sprintf("\n=== Regions with the largest %s-driven excess primary P (top %d of %d) ===\n",
            OUTLIER_SUBJECT, top_n, nrow(region_excess)))
print(head(top_regions[, c("region", "mean_p_excl", "sub11_p", "excess", "ratio")], 15))
cat(sprintf("\n%% of top-%d excess regions that are frontal (PFC/OFC/FEF/FrOper): %.1f%%\n", top_n, pct_frontal_top))
cat(sprintf("%% of all %d regions that are frontal (baseline rate): %.1f%%\n", nrow(region_excess), pct_frontal_all))
cat(sprintf("enrichment ratio: %.2fx\n", pct_frontal_top / pct_frontal_all))

cat("\nSaved:\n")
cat(" -", file.path(fig_dir, "outlier_sensitivity_p_histogram.png"), "\n")
cat(" -", file.path(fig_dir, "outlier_sensitivity_i_pri_n31.png"), "\n")
cat(" -", file.path(fig_dir, "outlier_sensitivity_p_pri_n31.png"), "\n")
cat(" -", file.path(fig_dir, "supp_fig4_outlier_sensitivity.png"), "\n")
cat(" -", file.path(outputs_dir, "primary_p", "outlier_sensitivity_subject_totals_p.csv"), "\n")
cat(" -", file.path(outputs_dir, "primary_p", "outlier_sensitivity_region_excess.csv"), "\n")
