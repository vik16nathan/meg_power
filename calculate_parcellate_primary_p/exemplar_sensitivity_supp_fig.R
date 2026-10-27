# Supplementary sensitivity analysis: are the 4 exemplar subjects used for
# the secondary-current FEM comparison (sub-0002, sub-0008, sub-0019,
# sub-0022 -- secondary-current FEM power is only computed for these 4, see
# secondary_p/README.md) representative of the remaining n=26 subjects in
# terms of primary current dipole moment (I_dip) and primary dipole power
# (P)? If not, the secondary-vs-primary comparison in Fig. 3 could be
# distorted by an unrepresentative subset.
#
# Panels:
#   A) Barplot: each of the 4 exemplar subjects' whole-brain total primary
#      dipole P, alongside the remaining n=26 subjects' mean +/- SEM.
#   B) Scatterplot: parcellated (Schaefer-600) average total primary P,
#      exemplar-4 mean per region (y) vs. remaining-26 mean per region (x),
#      with a line of unity.
#
# Also prints the same magnitude/spatial-correlation comparison for I_dip
# (not just P), since the main-text sentence this supports claims both
# quantities are similar between the two groups.
#
# Input (produced by parcellate_central_surface_s600.m):
# - s600_sub-<ID>_snr3_p_pri_2min_rest.csv (per-subject, per-region P)
# - s600_sub-<ID>_snr3_i_pri_2min_rest.csv (per-subject, per-region I_dip)

library(dplyr)
library(ggplot2)

outputs_dir <- "/export02/data/vikramn/hbm_manuscript_code/outputs"
primary_p_dir <- file.path(outputs_dir, "primary_p", "cent_surf_fem_fwd")
fig_dir <- file.path(outputs_dir, "primary_p", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

SNR <- 3

# sub-0013 excluded throughout this pipeline (DUNEuro forward-model failure);
# sub-0011 excluded as a confirmed outlier (snr_sensitivity_idip_rms.R). Same
# n=30 list used by every figure/statistic in this pipeline.
SubjectNames <- c(
  'sub-0002', 'sub-0008', 'sub-0009', 'sub-0010',
  'sub-0012', 'sub-0014', 'sub-0015',
  'sub-0016', 'sub-0018', 'sub-0019', 'sub-0020',
  'sub-0021', 'sub-0022', 'sub-0023', 'sub-0024', 'sub-0025',
  'sub-0026', 'sub-0028', 'sub-0029', 'sub-0030', 'sub-0031',
  'sub-0032', 'sub-0033', 'sub-0034', 'sub-0035', 'sub-0036',
  'sub-0037', 'sub-0039', 'sub-0040', 'sub-0041'
)

# The 4 exemplar subjects secondary-current FEM power was computed for.
EXEMPLAR_SUBJECTS <- c('sub-0002', 'sub-0008', 'sub-0019', 'sub-0022')
stopifnot(all(EXEMPLAR_SUBJECTS %in% SubjectNames))
REMAINING_SUBJECTS <- setdiff(SubjectNames, EXEMPLAR_SUBJECTS)
cat(sprintf("Exemplar subjects (n=%d): %s\n", length(EXEMPLAR_SUBJECTS), paste(EXEMPLAR_SUBJECTS, collapse = ", ")))
cat(sprintf("Remaining subjects (n=%d)\n", length(REMAINING_SUBJECTS)))

# =============================================================================
# Load all 30 subjects' per-region P and I_dip (SnrFixed=3)
# =============================================================================
load_long <- function(pattern) {
  files <- list.files(primary_p_dir, pattern = pattern, full.names = TRUE)
  files <- files[grepl(paste(sprintf("s600_%s_", SubjectNames), collapse = "|"), basename(files))]
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

found_subjects <- sort(unique(p_long$subject))
if (!setequal(found_subjects, SubjectNames)) {
  stop("Subject mismatch -- found ", length(found_subjects), " but expected ", length(SubjectNames),
       ". Missing: ", paste(setdiff(SubjectNames, found_subjects), collapse = ", "))
}

# =============================================================================
# Whole-brain totals per subject (for Panel A and the magnitude comparison)
# =============================================================================
subject_totals_p <- p_long %>% group_by(subject) %>% summarise(total_p = sum(value, na.rm = TRUE), .groups = "drop")
subject_totals_i <- i_long %>% group_by(subject) %>% summarise(total_i = sum(value, na.rm = TRUE), .groups = "drop")
subject_totals_p$group <- ifelse(subject_totals_p$subject %in% EXEMPLAR_SUBJECTS, "exemplar", "remaining")
subject_totals_i$group <- ifelse(subject_totals_i$subject %in% EXEMPLAR_SUBJECTS, "exemplar", "remaining")
write.csv(subject_totals_p, file.path(outputs_dir, "primary_p", "exemplar_sensitivity_subject_totals_p.csv"), row.names = FALSE)
write.csv(subject_totals_i, file.path(outputs_dir, "primary_p", "exemplar_sensitivity_subject_totals_i.csv"), row.names = FALSE)

# Magnitude comparison: exemplar-4 vs. remaining-26 means (+ Wilcoxon rank-
# sum -- n=4 is too small for a meaningful t-test, and rank-sum doesn't
# assume normality).
summarize_groups <- function(df, value_col) {
  exemplar_vals <- df[[value_col]][df$group == "exemplar"]
  remaining_vals <- df[[value_col]][df$group == "remaining"]
  wt <- wilcox.test(exemplar_vals, remaining_vals)
  list(
    exemplar_mean = mean(exemplar_vals), exemplar_sd = sd(exemplar_vals),
    remaining_mean = mean(remaining_vals), remaining_sem = sd(remaining_vals) / sqrt(length(remaining_vals)),
    pct_diff = 100 * (mean(exemplar_vals) - mean(remaining_vals)) / mean(remaining_vals),
    wilcox_W = wt$statistic, wilcox_p = wt$p.value
  )
}
p_group_stats <- summarize_groups(subject_totals_p, "total_p")
i_group_stats <- summarize_groups(subject_totals_i, "total_i")

cat(sprintf("\n=== Magnitude comparison: exemplar (n=4) vs. remaining (n=26), total P ===\n"))
cat(sprintf("  exemplar mean = %.3e W (SD = %.3e)\n", p_group_stats$exemplar_mean, p_group_stats$exemplar_sd))
cat(sprintf("  remaining mean = %.3e W (SEM = %.3e)\n", p_group_stats$remaining_mean, p_group_stats$remaining_sem))
cat(sprintf("  exemplar is %.1f%% %s than remaining\n", abs(p_group_stats$pct_diff),
            if (p_group_stats$pct_diff >= 0) "larger" else "smaller"))
cat(sprintf("  Wilcoxon rank-sum: W = %.1f, p = %.3g\n", p_group_stats$wilcox_W, p_group_stats$wilcox_p))

cat(sprintf("\n=== Magnitude comparison: exemplar (n=4) vs. remaining (n=26), total I_dip ===\n"))
cat(sprintf("  exemplar mean = %.3e A.m (SD = %.3e)\n", i_group_stats$exemplar_mean, i_group_stats$exemplar_sd))
cat(sprintf("  remaining mean = %.3e A.m (SEM = %.3e)\n", i_group_stats$remaining_mean, i_group_stats$remaining_sem))
cat(sprintf("  exemplar is %.1f%% %s than remaining\n", abs(i_group_stats$pct_diff),
            if (i_group_stats$pct_diff >= 0) "larger" else "smaller"))
cat(sprintf("  Wilcoxon rank-sum: W = %.1f, p = %.3g\n", i_group_stats$wilcox_W, i_group_stats$wilcox_p))

# =============================================================================
# Panel A: barplot -- each exemplar subject individually, plus the
# remaining-26 group mean +/- SEM, so it's visible whether the 4 specific
# exemplar subjects are typical or extreme, not just their group average.
# =============================================================================
remaining_bar <- data.frame(
  label = "Remaining\n(n=26 mean)", value = p_group_stats$remaining_mean,
  se = p_group_stats$remaining_sem, group = "remaining"
)
exemplar_bars <- subject_totals_p[subject_totals_p$group == "exemplar", c("subject", "total_p")]
exemplar_bars <- data.frame(label = exemplar_bars$subject, value = exemplar_bars$total_p, se = 0, group = "exemplar")
bar_df <- rbind(exemplar_bars, remaining_bar)
bar_df$label <- factor(bar_df$label, levels = c(sort(EXEMPLAR_SUBJECTS), "Remaining\n(n=26 mean)"))

barplot_A <- ggplot(bar_df, aes(x = label, y = value, fill = group)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.15) +
  scale_fill_manual(values = c(exemplar = "#4C72B0", remaining = "#DD8452")) +
  labs(x = NULL, y = "Whole-brain summed\nprimary dipole P (W)") +
  theme_minimal(base_size = 20) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = rel(0.75)),
    axis.text.y = element_text(size = rel(0.85))
  )
ggsave(file.path(fig_dir, "exemplar_sensitivity_p_barplot.png"), plot = barplot_A, width = 8, height = 7, bg = "white")

# =============================================================================
# Panel B: scatterplot -- parcellated (600-region) average total primary P,
# exemplar-4 mean per region vs. remaining-26 mean per region.
# =============================================================================
region_means <- function(long_df) {
  exemplar_means <- long_df[long_df$subject %in% EXEMPLAR_SUBJECTS, ] %>%
    group_by(region) %>% summarise(exemplar_mean = mean(value, na.rm = TRUE), .groups = "drop")
  remaining_means <- long_df[long_df$subject %in% REMAINING_SUBJECTS, ] %>%
    group_by(region) %>% summarise(remaining_mean = mean(value, na.rm = TRUE), .groups = "drop")
  merge(exemplar_means, remaining_means, by = "region")
}
p_region_df <- region_means(p_long)
i_region_df <- region_means(i_long)

p_spatial_r <- cor(p_region_df$remaining_mean, p_region_df$exemplar_mean, method = "spearman")
i_spatial_r <- cor(i_region_df$remaining_mean, i_region_df$exemplar_mean, method = "spearman")
cat(sprintf("\n=== Spatial correlation (Spearman r across %d regions), exemplar-4 mean vs. remaining-26 mean ===\n", nrow(p_region_df)))
cat(sprintf("  Primary dipole P:   r = %.3f\n", p_spatial_r))
cat(sprintf("  Primary I_dip:      r = %.3f\n", i_spatial_r))

line_range <- range(c(p_region_df$remaining_mean, p_region_df$exemplar_mean))
line_df <- data.frame(x = line_range, y = line_range)

scatter_B <- ggplot(p_region_df, aes(x = remaining_mean, y = exemplar_mean)) +
  geom_point(alpha = 0.4, size = 1.5, color = "#4C72B0") +
  geom_line(data = line_df, aes(x = x, y = y), color = "black", linewidth = 1, linetype = "dashed", inherit.aes = FALSE) +
  annotate("text", x = min(p_region_df$remaining_mean), y = max(p_region_df$exemplar_mean),
           label = sprintf("Spearman r = %.2f\n(n=%d regions)", p_spatial_r, nrow(p_region_df)),
           hjust = -0.05, vjust = 1.1, size = 6, color = "darkred", fontface = "bold") +
  labs(x = "Remaining (n=26) mean\nprimary dipole P (W)",
       y = "Exemplar (n=4) mean\nprimary dipole P (W)") +
  theme_minimal(base_size = 22) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = rel(0.8))
  )
ggsave(file.path(fig_dir, "exemplar_sensitivity_p_scatter.png"), plot = scatter_B, width = 7, height = 7, bg = "white")

cat("\nSaved:\n")
cat(" -", file.path(outputs_dir, "primary_p", "exemplar_sensitivity_subject_totals_p.csv"), "\n")
cat(" -", file.path(outputs_dir, "primary_p", "exemplar_sensitivity_subject_totals_i.csv"), "\n")
cat(" -", file.path(fig_dir, "exemplar_sensitivity_p_barplot.png"), "\n")
cat(" -", file.path(fig_dir, "exemplar_sensitivity_p_scatter.png"), "\n")
