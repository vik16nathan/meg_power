# Sensitivity analysis: effect of the MNE regularization parameter
# (SnrFixed, swept by rms_idip_rest.m: 1, 3, 5) directly on the reconstructed
# per-vertex RMS primary current dipole moment (I_dip) -- rms_idip_rest.m's
# raw output, upstream of cortical-column power and Schaefer parcellation.
#
# Power P = I_dip^2 * R (R = cortical column resistance, SnrFixed-independent
# -- see process_brainstorm_data/parcellate_cortical_column_res_vol.m), so any
# SnrFixed effect on P follows directly from its effect on I_dip; a full
# sensitivity analysis does not require re-running calculate_fem_column_power.m
# / parcellate_central_surface_s600.m at every SnrFixed value.
#
# Input (produced by calculate_parcellate_primary_p/rms_idip_rest.m):
# - <subject>_snr<N>_idip_rms.csv (per-vertex [x,y,z] RMS I_dip, A.m)

library(dplyr)
library(ggplot2)
library(R.matlab) # for reading calculate_fem_column_power.m's dip_p_i_rms.mat output
library(patchwork) # for the combined 4-panel figure at the end of this script

outputs_dir <- "/export02/data/vikramn/hbm_manuscript_code/outputs"
rms_dir <- file.path(outputs_dir, "primary_p", "idip_rms") # rms_idip_rest.m output
power_dir <- file.path(outputs_dir, "primary_p", "fem_column_power") # calculate_fem_column_power.m output
primary_p_dir <- file.path(outputs_dir, "primary_p")
fig_dir <- file.path(primary_p_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# calculate_fem_column_power.m saves 'dip_p'/'dip_i'; R.matlab::readMat()
# mangles underscores to dots when returning MAT variables as list names
# (verified directly against these files: names(readMat(...)) == c("dip.i",
# "dip.p")), so every access below uses d$dip.p / d$dip.i, not d$dip_p.
load_dip_p <- function(matfile) as.numeric(readMat(matfile)$dip.p)

# sub-0013 excluded throughout this pipeline (DUNEuro forward-model failed
# even after widening the BEM brain-inner margin to 3/6/9mm -- see
# SUPPLEMENTARY_METHODS.md Section 1). sub-0011 excluded here too: confirmed
# outlier (~11x the remaining-cohort mean whole-brain current dipole moment
# at SnrFixed=3, R and cortical thickness both normal -- flagged but not
# root-caused, likely subject-specific data quality rather than a pipeline
# bug). Final analyzed sample for every figure/statistic in this script is
# n=30 -- matches parcellate_central_surface_s600.m's SubjectNames.
SubjectNames <- c(
  'sub-0002', 'sub-0008', 'sub-0009', 'sub-0010',
  'sub-0012', 'sub-0014', 'sub-0015',
  'sub-0016', 'sub-0018', 'sub-0019', 'sub-0020',
  'sub-0021', 'sub-0022', 'sub-0023', 'sub-0024', 'sub-0025',
  'sub-0026', 'sub-0028', 'sub-0029', 'sub-0030', 'sub-0031',
  'sub-0032', 'sub-0033', 'sub-0034', 'sub-0035', 'sub-0036',
  'sub-0037', 'sub-0039', 'sub-0040', 'sub-0041'
)
SnrValues <- c(1, 3, 5) # must match the sweep in rms_idip_rest.m

# Kept as a named constant purely for documentation/traceability (also
# referenced in plot_total_p_histogram.py) -- SubjectNames above already
# excludes it, so nothing below needs to filter by this.
OUTLIER_SUBJECT <- 'sub-0011'

# =============================================================================
# Load per-subject, per-SnrFixed whole-brain I_dip summaries
# =============================================================================
rows <- list()
for (subject in SubjectNames) {
  for (snr in SnrValues) {
    f <- file.path(rms_dir, sprintf('%s_snr%d_idip_rms.csv', subject, snr))
    if (!file.exists(f)) {
      stop("Missing ", f, " -- rms_idip_rest.m has not finished this subject/SnrFixed yet.")
    }
    d <- as.matrix(read.csv(f, header = FALSE))
    # Each row is one cortical vertex's RMS current dipole moment [x, y, z]
    # (A.m); per-vertex magnitude, matching
    # calculate_fem_column_power.m's dip_i = vecnorm(curr_xyz, 2, 3).
    vertex_mag <- sqrt(rowSums(d^2))
    rows[[length(rows) + 1]] <- data.frame(
      subject = subject, snr = snr,
      total_i = sum(vertex_mag), mean_i = mean(vertex_mag), n_vertices = nrow(d))
  }
}
idip_df <- do.call(rbind, rows)
write.csv(idip_df, file.path(primary_p_dir, 'snr_sensitivity_idip_rms_per_subject.csv'), row.names = FALSE)

# =============================================================================
# Summary across subjects, per SnrFixed
# =============================================================================
summary_df <- idip_df %>%
  group_by(snr) %>%
  summarise(mean_total_i = mean(total_i), sd_total_i = sd(total_i), n = n(), .groups = "drop") %>%
  mutate(se_total_i = sd_total_i / sqrt(n))
write.csv(summary_df, file.path(primary_p_dir, 'snr_sensitivity_idip_rms_summary.csv'), row.names = FALSE)

# Main-text figure: max/min ratio across SnrFixed, in orders of magnitude
# (log10 of the ratio) -- reported directly rather than left for the reader
# to compute from the table above.
oom_ratio <- max(summary_df$mean_total_i) / min(summary_df$mean_total_i)
cat(sprintf("\n=== Main-text summary: SnrFixed effect on total primary current (n=%d subjects, %s excluded) ===\n",
            length(SubjectNames), OUTLIER_SUBJECT))
for (i in seq_len(nrow(summary_df))) {
  cat(sprintf("  SnrFixed=%d: mean total I_dip = %.2e A.m (SE = %.2e)\n",
              summary_df$snr[i], summary_df$mean_total_i[i], summary_df$se_total_i[i]))
}
cat(sprintf("  max/min ratio = %.2fx (%.2f orders of magnitude, log10 scale)\n", oom_ratio, log10(oom_ratio)))

# Same panel sizing/styling as p_barplot_panel below (fig2d_res_barplot_subjavg.png
# convention): no title, bold wrapped axis text sized for a composite figure.
barplot_panel <- ggplot(summary_df, aes(x = factor(snr), y = mean_total_i)) +
  geom_col(fill = "#4C72B0", width = 0.6) +
  geom_errorbar(aes(ymin = mean_total_i - se_total_i, ymax = mean_total_i + se_total_i), width = 0.15) +
  labs(x = "SnrFixed (MNE\nregularization parameter)",
       y = "Whole-brain summed\nRMS I_dip (A.m)") +
  theme_minimal(base_size = 24) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = rel(0.85))
  )
ggsave(file.path(fig_dir, 'snr_sensitivity_idip_rms_barplot.png'), plot = barplot_panel, width = 6, height = 7)

# =============================================================================
# Repeated-measures ANOVA: same n=30 subjects (sub-0011 excluded as an
# outlier) measured at all 3 SnrFixed levels, testing whether SnrFixed
# significantly affects whole-brain summed I_dip.
# =============================================================================
idip_df$subject <- factor(idip_df$subject)
idip_df$snr <- factor(idip_df$snr)
rm_anova <- aov(total_i ~ snr + Error(subject/snr), data = idip_df)
anova_summary <- summary(rm_anova)

sink(file.path(primary_p_dir, 'snr_sensitivity_idip_rms_anova.txt'))
cat("Repeated-measures ANOVA: SnrFixed effect on whole-brain summed RMS I_dip\n")
cat(sprintf("n = %d subjects, SnrFixed levels: %s\n\n", length(SubjectNames), paste(SnrValues, collapse = ", ")))
print(anova_summary)
sink()

cat("=== Summary across SnrFixed (n=", length(SubjectNames), " subjects) ===\n", sep = "")
print(summary_df)
cat("\n=== Repeated-measures ANOVA ===\n")
print(anova_summary)
cat("\nSaved:\n")
cat(" -", file.path(primary_p_dir, 'snr_sensitivity_idip_rms_per_subject.csv'), "\n")
cat(" -", file.path(primary_p_dir, 'snr_sensitivity_idip_rms_summary.csv'), "\n")
cat(" -", file.path(primary_p_dir, 'snr_sensitivity_idip_rms_anova.txt'), "\n")
cat(" -", file.path(fig_dir, 'snr_sensitivity_idip_rms_barplot.png'), "\n")

# =============================================================================
# Same SnrFixed sensitivity, but for total primary dipole POWER (P) instead
# of I_dip -- read directly from calculate_fem_column_power.m's dip_p output
# (P = (I_dip/thickness)^2 * R). Already computed for all three SnrFixed
# values (that script's SnrValues was widened to [1,3,5] earlier in this
# pipeline), so this is a read of existing files, not a new computation.
# =============================================================================
p_rows <- list()
for (subject in SubjectNames) {
  for (snr in SnrValues) {
    f <- file.path(power_dir, sprintf('%s_snr%d_dip_p_i_rms.mat', subject, snr))
    if (!file.exists(f)) {
      stop("Missing ", f, " -- calculate_fem_column_power.m has not finished this subject/SnrFixed yet.")
    }
    dip_p <- load_dip_p(f)
    p_rows[[length(p_rows) + 1]] <- data.frame(
      subject = subject, snr = snr,
      total_p = sum(dip_p), mean_p = mean(dip_p), n_vertices = length(dip_p))
  }
}
p_df <- do.call(rbind, p_rows)
write.csv(p_df, file.path(primary_p_dir, 'snr_sensitivity_p_per_subject.csv'), row.names = FALSE)

p_summary_df <- p_df %>%
  group_by(snr) %>%
  summarise(mean_total_p = mean(total_p), sd_total_p = sd(total_p), n = n(), .groups = "drop") %>%
  mutate(se_total_p = sd_total_p / sqrt(n))
write.csv(p_summary_df, file.path(primary_p_dir, 'snr_sensitivity_p_summary.csv'), row.names = FALSE)

oom_ratio_p <- max(p_summary_df$mean_total_p) / min(p_summary_df$mean_total_p)
cat(sprintf("\n=== Main-text summary: SnrFixed effect on total primary dipole POWER (n=%d subjects, %s excluded) ===\n",
            length(SubjectNames), OUTLIER_SUBJECT))
for (i in seq_len(nrow(p_summary_df))) {
  cat(sprintf("  SnrFixed=%d: mean total P = %.2e W (SE = %.2e)\n",
              p_summary_df$snr[i], p_summary_df$mean_total_p[i], p_summary_df$se_total_p[i]))
}
cat(sprintf("  max/min ratio = %.2fx (%.2f orders of magnitude, log10 scale)\n", oom_ratio_p, log10(oom_ratio_p)))

# Sized/styled to match the fig2d_res_barplot_subjavg.png convention
# (res_correlations.R): no title (a composite figure supplies its own panel
# caption; a title this long was getting clipped when rendered small), bold
# larger axis text so it stays legible once shrunk into a multi-panel figure.
p_barplot_panel <- ggplot(p_summary_df, aes(x = factor(snr), y = mean_total_p)) +
  geom_col(fill = "#DD8452", width = 0.6) +
  geom_errorbar(aes(ymin = mean_total_p - se_total_p, ymax = mean_total_p + se_total_p), width = 0.15) +
  labs(x = "SnrFixed (MNE\nregularization parameter)",
       y = "Whole-brain summed\nprimary dipole P (W)") +
  theme_minimal(base_size = 24) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = rel(0.85))
  )
ggsave(file.path(fig_dir, 'snr_sensitivity_p_barplot.png'), plot = p_barplot_panel, width = 6, height = 7)

p_df$subject <- factor(p_df$subject)
p_df$snr <- factor(p_df$snr)
rm_anova_p <- aov(total_p ~ snr + Error(subject/snr), data = p_df)
anova_summary_p <- summary(rm_anova_p)

sink(file.path(primary_p_dir, 'snr_sensitivity_p_anova.txt'))
cat("Repeated-measures ANOVA: SnrFixed effect on whole-brain summed primary dipole P\n")
cat(sprintf("n = %d subjects, SnrFixed levels: %s\n\n", length(SubjectNames), paste(SnrValues, collapse = ", ")))
print(anova_summary_p)
sink()

cat("=== Summary across SnrFixed, total P (n=", length(SubjectNames), " subjects) ===\n", sep = "")
print(p_summary_df)
cat("\n=== Repeated-measures ANOVA (P) ===\n")
print(anova_summary_p)
cat("\nSaved:\n")
cat(" -", file.path(primary_p_dir, 'snr_sensitivity_p_per_subject.csv'), "\n")
cat(" -", file.path(primary_p_dir, 'snr_sensitivity_p_summary.csv'), "\n")
cat(" -", file.path(primary_p_dir, 'snr_sensitivity_p_anova.txt'), "\n")
cat(" -", file.path(fig_dir, 'snr_sensitivity_p_barplot.png'), "\n")

# =============================================================================
# Targeted sensitivity analysis: CONSTRAINED (fixed, normal-to-cortex) vs.
# UNCONSTRAINED (free, 3D) dipole orientation, at SnrFixed=3 only -- NOT a
# SnrFixed sweep (see process_brainstorm_data/process_meg_compute_mn_source_fem_constrained.m
# and rms_idip_rest_constrained.m, which compute the constrained kernel/RMS).
#
# Hypothesis: unconstrained I_dip RMS should be >= constrained I_dip RMS,
# since unconstrained is the RMS-norm of the full 3D dipole vector while
# constrained is a single (normal-direction) component of essentially the
# same underlying current -- a vector's norm can only be >= any one of its
# own components' magnitude.
#
# Input (produced by rms_idip_rest_constrained.m):
# - <subject>_snr3_idip_rms_constrained.csv (per-vertex signed RMS I_dip
#   along the cortical normal, A.m -- 1 column, not [x,y,z])
# =============================================================================
orientation_rows <- list()
for (subject in SubjectNames) {
  f <- file.path(rms_dir, sprintf('%s_snr3_idip_rms_constrained.csv', subject))
  if (!file.exists(f)) {
    stop("Missing ", f, " -- rms_idip_rest_constrained.m has not finished this subject yet.")
  }
  d <- as.matrix(read.csv(f, header = FALSE))
  vertex_mag <- abs(d[, 1]) # RMS is already >=0 by construction; abs() is a no-op safeguard
  orientation_rows[[length(orientation_rows) + 1]] <- data.frame(
    subject = subject, orientation = "constrained",
    total_i = sum(vertex_mag), mean_i = mean(vertex_mag), n_vertices = length(vertex_mag))
}
orientation_df <- do.call(rbind, orientation_rows)

# Reuse the unconstrained SnrFixed=3 rows already loaded above (idip_df) --
# same subjects, same SnrFixed, only the orientation constraint differs.
unconstrained_snr3 <- idip_df[idip_df$snr == 3, ]
unconstrained_snr3$orientation <- "unconstrained"
col_order <- c("subject", "orientation", "total_i", "mean_i", "n_vertices")
orientation_df <- rbind(orientation_df[, col_order], unconstrained_snr3[, col_order])
write.csv(orientation_df, file.path(primary_p_dir, 'orientation_sensitivity_idip_rms_per_subject.csv'), row.names = FALSE)

orientation_summary <- orientation_df %>%
  group_by(orientation) %>%
  summarise(mean_total_i = mean(total_i), sd_total_i = sd(total_i), n = n(), .groups = "drop") %>%
  mutate(se_total_i = sd_total_i / sqrt(n))
write.csv(orientation_summary, file.path(primary_p_dir, 'orientation_sensitivity_idip_rms_summary.csv'), row.names = FALSE)

# Same panel sizing/styling as p_orientation_barplot below.
orientation_barplot <- ggplot(orientation_summary, aes(x = orientation, y = mean_total_i, fill = orientation)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_total_i - se_total_i, ymax = mean_total_i + se_total_i), width = 0.15) +
  labs(x = "Dipole orientation\nconstraint (SnrFixed=3)",
       y = "Whole-brain summed\nRMS I_dip (A.m)") +
  theme_minimal(base_size = 24) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = rel(0.85))
  )
ggsave(file.path(fig_dir, 'orientation_sensitivity_idip_rms_barplot.png'), plot = orientation_barplot, width = 6, height = 7)

# Paired test (same subjects measured both ways at SnrFixed=3) -- avoids
# tidyr (not used elsewhere in this pipeline) via a plain merge() instead of
# pivot_wider().
wide <- merge(
  orientation_df[orientation_df$orientation == "unconstrained", c("subject", "total_i")],
  orientation_df[orientation_df$orientation == "constrained", c("subject", "total_i")],
  by = "subject", suffixes = c("_unconstrained", "_constrained")
)
paired_test <- t.test(wide$total_i_unconstrained, wide$total_i_constrained, paired = TRUE)
pct_diff <- 100 * (mean(wide$total_i_unconstrained) - mean(wide$total_i_constrained)) / mean(wide$total_i_constrained)

sink(file.path(primary_p_dir, 'orientation_sensitivity_idip_rms_ttest.txt'))
cat("Paired t-test: orientation constraint effect on whole-brain summed RMS I_dip (SnrFixed=3)\n")
cat(sprintf("n = %d subjects\n", nrow(wide)))
cat(sprintf("unconstrained is %.1f%% %s than constrained on average\n\n",
            abs(pct_diff), if (pct_diff >= 0) "larger" else "smaller"))
print(paired_test)
sink()

cat("\n=== Orientation sensitivity (constrained vs. unconstrained, SnrFixed=3, n=", nrow(wide), " subjects) ===\n", sep = "")
print(orientation_summary)
cat(sprintf("\nunconstrained is %.1f%% %s than constrained on average\n",
            abs(pct_diff), if (pct_diff >= 0) "larger" else "smaller"))
cat("\n=== Paired t-test ===\n")
print(paired_test)
cat("\nSaved:\n")
cat(" -", file.path(primary_p_dir, 'orientation_sensitivity_idip_rms_per_subject.csv'), "\n")
cat(" -", file.path(primary_p_dir, 'orientation_sensitivity_idip_rms_summary.csv'), "\n")
cat(" -", file.path(primary_p_dir, 'orientation_sensitivity_idip_rms_ttest.txt'), "\n")
cat(" -", file.path(fig_dir, 'orientation_sensitivity_idip_rms_barplot.png'), "\n")

# =============================================================================
# Same orientation comparison, but for total primary dipole POWER (P)
# instead of I_dip -- read from calculate_fem_column_power_constrained.m's
# dip_p output (same P = (I_dip/thickness)^2 * R formula, computed directly
# from the constrained I_dip data plus the SAME thickness/R used for the
# unconstrained case -- both are orientation-invariant properties of the
# cortical column geometry, not the source estimate).
# =============================================================================
p_orientation_rows <- list()
for (subject in SubjectNames) {
  f <- file.path(power_dir, sprintf('%s_snr3_dip_p_i_rms_constrained.mat', subject))
  if (!file.exists(f)) {
    stop("Missing ", f, " -- calculate_fem_column_power_constrained.m has not finished this subject yet.")
  }
  dip_p <- load_dip_p(f)
  p_orientation_rows[[length(p_orientation_rows) + 1]] <- data.frame(
    subject = subject, orientation = "constrained",
    total_p = sum(dip_p), mean_p = mean(dip_p), n_vertices = length(dip_p))
}
p_orientation_df <- do.call(rbind, p_orientation_rows)

unconstrained_p_snr3 <- p_df[p_df$snr == 3, ]
unconstrained_p_snr3$orientation <- "unconstrained"
p_col_order <- c("subject", "orientation", "total_p", "mean_p", "n_vertices")
p_orientation_df <- rbind(p_orientation_df[, p_col_order], unconstrained_p_snr3[, p_col_order])
write.csv(p_orientation_df, file.path(primary_p_dir, 'orientation_sensitivity_p_per_subject.csv'), row.names = FALSE)

p_orientation_summary <- p_orientation_df %>%
  group_by(orientation) %>%
  summarise(mean_total_p = mean(total_p), sd_total_p = sd(total_p), n = n(), .groups = "drop") %>%
  mutate(se_total_p = sd_total_p / sqrt(n))
write.csv(p_orientation_summary, file.path(primary_p_dir, 'orientation_sensitivity_p_summary.csv'), row.names = FALSE)

# Same panel sizing/styling as p_barplot_panel above.
p_orientation_barplot <- ggplot(p_orientation_summary, aes(x = orientation, y = mean_total_p, fill = orientation)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_total_p - se_total_p, ymax = mean_total_p + se_total_p), width = 0.15) +
  labs(x = "Dipole orientation\nconstraint (SnrFixed=3)",
       y = "Whole-brain summed\nprimary dipole P (W)") +
  theme_minimal(base_size = 24) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = rel(0.85))
  )
ggsave(file.path(fig_dir, 'orientation_sensitivity_p_barplot.png'), plot = p_orientation_barplot, width = 6, height = 7)

p_wide <- merge(
  p_orientation_df[p_orientation_df$orientation == "unconstrained", c("subject", "total_p")],
  p_orientation_df[p_orientation_df$orientation == "constrained", c("subject", "total_p")],
  by = "subject", suffixes = c("_unconstrained", "_constrained")
)
p_paired_test <- t.test(p_wide$total_p_unconstrained, p_wide$total_p_constrained, paired = TRUE)
p_pct_diff <- 100 * (mean(p_wide$total_p_unconstrained) - mean(p_wide$total_p_constrained)) / mean(p_wide$total_p_constrained)

sink(file.path(primary_p_dir, 'orientation_sensitivity_p_ttest.txt'))
cat("Paired t-test: orientation constraint effect on whole-brain summed primary dipole P (SnrFixed=3)\n")
cat(sprintf("n = %d subjects\n", nrow(p_wide)))
cat(sprintf("unconstrained is %.1f%% %s than constrained on average\n\n",
            abs(p_pct_diff), if (p_pct_diff >= 0) "larger" else "smaller"))
print(p_paired_test)
sink()

cat("\n=== Orientation sensitivity, total P (constrained vs. unconstrained, SnrFixed=3, n=", nrow(p_wide), " subjects) ===\n", sep = "")
print(p_orientation_summary)
cat(sprintf("\nunconstrained is %.1f%% %s than constrained on average\n",
            abs(p_pct_diff), if (p_pct_diff >= 0) "larger" else "smaller"))
cat("\n=== Paired t-test (P) ===\n")
print(p_paired_test)
cat("\nSaved:\n")
cat(" -", file.path(primary_p_dir, 'orientation_sensitivity_p_per_subject.csv'), "\n")
cat(" -", file.path(primary_p_dir, 'orientation_sensitivity_p_summary.csv'), "\n")
cat(" -", file.path(primary_p_dir, 'orientation_sensitivity_p_ttest.txt'), "\n")
cat(" -", file.path(fig_dir, 'orientation_sensitivity_p_barplot.png'), "\n")

# =============================================================================
# Combined 4-panel figure: A = I_dip (i: SnrFixed sweep, ii: orientation
# constraint), B = primary dipole P (same two analyses) -- lets a reader
# compare the SnrFixed and orientation sensitivity results for both
# quantities at a glance instead of as four separate standalone files. The
# four standalone PNGs above are still saved independently too.
# =============================================================================
tag_style <- theme(plot.title = element_text(face = "bold", size = 20, hjust = 0))

panel_A_i   <- barplot_panel         + ggtitle("A(i)")  + tag_style
panel_A_ii  <- orientation_barplot   + ggtitle("A(ii)") + tag_style
panel_B_i   <- p_barplot_panel       + ggtitle("B(i)")  + tag_style
panel_B_ii  <- p_orientation_barplot + ggtitle("B(ii)") + tag_style

combined_4panel <- (panel_A_i | panel_A_ii) / (panel_B_i | panel_B_ii)
ggsave(file.path(fig_dir, 'snr_orientation_sensitivity_4panel.png'), plot = combined_4panel,
       width = 12, height = 14, bg = "white")
cat(" -", file.path(fig_dir, 'snr_orientation_sensitivity_4panel.png'), "\n")
