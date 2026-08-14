# Dendritic/neuron-density-derived resistance vs. uniform grey-matter
# columnar resistance (Fig. 3D). One point per subject (subject-mean R under
# each formulation), not one point per cortical column pooled across
# subjects -- pooling ~15,000 non-independent columns per subject into a
# single correlation treats each subject's columns as independent samples
# and inflates the apparent correlation (pseudoreplication). Subject-wise
# correlations (each subject's own columns) are computed separately and
# summarized (median/range) as a check on whether the relationship holds
# within-subject, not just across subject averages.
#
# IMPORTANT CONFOUND: both formulations divide by the same per-vertex support
# area A (cell_res_avg.m: vertex_res_map_old = gm_res*thickness/A;
# vertex_resistance_map = dend_r/(neuronDensityValue*1e6*A)) -- so part of
# their correlation is mechanical (both get large/small together wherever A
# is small/large), not evidence that thickness and neuron density themselves
# agree. The "A-removed" analysis below multiplies each R by its own vertex's
# A to cancel this shared 1/A term (R_old*A = gm_res*thickness;
# R_new*A = dend_r/(neuronDensityValue*1e6), a region-level constant) and
# correlates what's left.
library(R.matlab)
library(ggplot2)

RESISTANCES_DIR <- "/export02/data/vikramn/hbm_manuscript_code/outputs/resistances"
OUT_DIR <- "/export02/data/vikramn/hbm_manuscript_code/outputs/primary_p/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

subjects <- sub("_region_res_cell_column\\.mat$", "",
                 basename(Sys.glob(file.path(RESISTANCES_DIR, "*_region_res_cell_column.mat"))))

# sub-0011 excluded: confirmed outlier (~11x the remaining-cohort mean
# whole-brain current dipole moment at SnrFixed=3) flagged in
# calculate_parcellate_primary_p/snr_sensitivity_idip_rms.R and
# plot_total_p_histogram.py -- excluded here too for consistency across every
# figure/statistic that reports subject-level summaries. Final n=30.
OUTLIER_SUBJECT <- "sub-0011"
subjects <- setdiff(subjects, OUTLIER_SUBJECT)
cat(sprintf("Found %d subjects (excluding %s): %s\n", length(subjects), OUTLIER_SUBJECT, paste(subjects, collapse = ", ")))

# Per-subject: subject-mean R under each formulation, this subject's own
# within-subject Spearman correlation across its columns (raw R, still
# confounded by shared A), the A-removed within-subject correlation (R*A on
# both sides, confound removed), and subject-mean A-removed values (for the
# subject-averaged A-removed scatter below).
subj_rows <- list()
for (subject in subjects) {
  mat_data <- readMat(file.path(RESISTANCES_DIR, paste0(subject, "_region_res_cell_column.mat")))
  vertex_res_map_old <- as.numeric(mat_data[[1]])
  vertex_resistance_map <- as.numeric(mat_data[[2]])
  support_areas <- as.numeric(as.matrix(read.csv(
    file.path(RESISTANCES_DIR, paste0(subject, "_support_areas.csv")), header = FALSE)))

  d <- data.frame(vertex_res_map_old = vertex_res_map_old,
                   vertex_resistance_map = vertex_resistance_map,
                   support_areas = support_areas)
  d <- d[which(d$vertex_resistance_map != 0), ] # filter zeroes, same as original
  d$a_removed_old <- d$vertex_res_map_old * d$support_areas       # = gm_res * thickness
  d$a_removed_new <- d$vertex_resistance_map * d$support_areas    # = dend_r / (neuronDensityValue * 1e6)

  subj_rows[[subject]] <- data.frame(
    subject = subject,
    mean_old = mean(d$vertex_res_map_old, na.rm = TRUE),
    mean_new = mean(d$vertex_resistance_map, na.rm = TRUE),
    within_subj_r = cor(d$vertex_res_map_old, d$vertex_resistance_map, method = "spearman"),
    within_subj_r_a_removed = cor(d$a_removed_old, d$a_removed_new, method = "spearman"),
    mean_a_removed_old = mean(d$a_removed_old, na.rm = TRUE),
    mean_a_removed_new = mean(d$a_removed_new, na.rm = TRUE))
}
subj_df <- do.call(rbind, subj_rows)
rownames(subj_df) <- NULL
write.csv(subj_df, file.path(OUT_DIR, "fig3d_res_correlation_per_subject.csv"), row.names = FALSE)

cat("\nPer-subject mean resistances, within-subject correlation, and A-removed within-subject correlation:\n")
print(subj_df)

### Subject-level summary stats (raw R, shared-A confound still present) ###
subj_level_r <- cor(subj_df$mean_old, subj_df$mean_new, method = "spearman")
within_subj_median <- median(subj_df$within_subj_r)
within_subj_range <- range(subj_df$within_subj_r)

### Paired comparison of subject-level mean resistances across formulations ###
paired_test <- wilcox.test(subj_df$mean_old, subj_df$mean_new, paired = TRUE)
paired_diff <- subj_df$mean_new - subj_df$mean_old

cat(sprintf("\nSubject-level (n=%d) Spearman r of subject-mean R = %.3f\n", nrow(subj_df), subj_level_r))
cat(sprintf("Within-subject Spearman r: median = %.3f (range %.3f-%.3f)\n",
            within_subj_median, within_subj_range[1], within_subj_range[2]))
cat(sprintf("Paired Wilcoxon signed-rank (mean_old vs. mean_new): V = %.1f, p = %.3g\n",
            paired_test$statistic, paired_test$p.value))
cat(sprintf("Median paired difference (new - old) = %.3f Ohms\n", median(paired_diff)))

# =============================================================================
# COMMENTED OUT: subject-mean spatial-comparison scatter (superseded here by
# the A-removed analysis below; code kept for reference, not run).
# =============================================================================
# annotation_label <- paste(
#   sprintf("Subject-level r = %.2f (n=%d)", subj_level_r, nrow(subj_df)),
#   sprintf("Within-subject r: median = %.2f (range %.2f-%.2f)", within_subj_median, within_subj_range[1], within_subj_range[2]),
#   sprintf("Paired Wilcoxon p = %.3g", paired_test$p.value),
#   sep = "\n"
# )
#
# p <- ggplot(subj_df, aes(x = mean_old, y = mean_new)) +
#   geom_point(size = 3, color = "#4C72B0") +
#   labs(
#     x = "Subject-mean R, GM resistivity in\ncortical column (Ohms)",
#     y = "Subject-mean R, dendrite and\nneuron density (Ohms)"
#   ) +
#   annotate(
#     "text",
#     x = min(subj_df$mean_old, na.rm = TRUE),
#     y = max(subj_df$mean_new, na.rm = TRUE),
#     label = annotation_label,
#     hjust = -0.05, vjust = 1.1, size = 5.5, color = "darkred", fontface = "bold"
#   ) +
#   theme_minimal(base_size = 22) +
#   theme(
#     legend.position = "none",
#     axis.title = element_text(face = "bold"),
#     axis.text = element_text(size = rel(0.8))
#   )
#
# ggsave(file.path(OUT_DIR, "fig3d_res_correlation_scatter_subjavg.png"), plot = p, width = 7, height = 7)

#### barplot: mean +/- SEM of subject-level means, across n=30 subjects ####
# SEM here is across subjects (independent samples), not across pooled
# columns -- pooling ~15,000 non-independent columns per subject (the
# original approach) treats them as independent observations and understates
# the true SEM by roughly sqrt(n_columns/n_subjects).
n_subj <- nrow(subj_df)
mean_old <- mean(subj_df$mean_old, na.rm = TRUE)
mean_new <- mean(subj_df$mean_new, na.rm = TRUE)
sem_old <- sd(subj_df$mean_old, na.rm = TRUE) / sqrt(n_subj)
sem_new <- sd(subj_df$mean_new, na.rm = TRUE) / sqrt(n_subj)

bar_df <- data.frame(
  Axis = c("Cortical column R", "Cell resistance R"),
  Mean_R = c(mean_old, mean_new),
  SEM = c(sem_old, sem_new)
)

# Same panel sizing/styling as the other Fig 3D panels: no title, larger bold
# axis text so it holds up when shrunk into the composite figure.
p_bar <- ggplot(bar_df, aes(x = Axis, y = Mean_R, fill = Axis)) +
  geom_bar(stat = "identity", width = 0.6, show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = Mean_R - SEM, ymax = Mean_R + SEM),
    width = 0.15,
    linewidth = 1
  ) +
  labs(
    y = "Mean Resistance R (Ohms)",
    x = NULL
  ) +
  theme_minimal(base_size = 24) +
  theme(
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold", size = rel(0.85)),
    axis.text.y = element_text(size = rel(0.85))
  )

ggsave(file.path(OUT_DIR, "fig3d_res_barplot_subjavg.png"), plot = p_bar, width = 6, height = 7)

# =============================================================================
# NEW: subject-averaged A-removed scatter -- one point per subject (n=30),
# each subject's own vertex-wise A-removed values (a_removed_old,
# a_removed_new) averaged first, THEN correlated across subjects. Tests
# whether thickness (a_removed_old) and region-level neuron density
# (a_removed_new) agree at the subject level once the shared 1/A term is no
# longer inflating the raw-R correlation above.
# =============================================================================
a_removed_subjavg_r <- cor(subj_df$mean_a_removed_old, subj_df$mean_a_removed_new, method = "spearman")
cat(sprintf("\nSubject-averaged A-removed Spearman r = %.3f (n=%d subjects)\n",
            a_removed_subjavg_r, nrow(subj_df)))

p_a_removed_subjavg <- ggplot(subj_df, aes(x = mean_a_removed_old, y = mean_a_removed_new)) +
  geom_point(size = 3, color = "#4C72B0") +
  labs(
    x = "Subject-mean GM resistivity term\nwith A removed (gm_res × thickness, Ω·m²)",
    y = "Subject-mean cell-derived term\nwith A removed (dend_r / neuron density, Ω·m²)"
  ) +
  annotate(
    "text",
    x = min(subj_df$mean_a_removed_old, na.rm = TRUE),
    y = max(subj_df$mean_a_removed_new, na.rm = TRUE),
    label = sprintf("Subject-averaged\nSpearman r = %.2f (n=%d)", a_removed_subjavg_r, nrow(subj_df)),
    hjust = -0.05, vjust = 1.1, size = 5.5, color = "darkred", fontface = "bold"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold", size = rel(0.6)),
    axis.text = element_text(size = rel(0.8)),
    plot.margin = margin(t = 10, r = 16, b = 10, l = 10)
  )

ggsave(file.path(OUT_DIR, "fig3d_res_a_removed_scatter_subjavg.png"), plot = p_a_removed_subjavg, width = 7, height = 7)

cat(sprintf("\nSaved:\n  %s\n  %s\n  %s\n",
            file.path(OUT_DIR, "fig3d_res_correlation_per_subject.csv"),
            file.path(OUT_DIR, "fig3d_res_barplot_subjavg.png"),
            file.path(OUT_DIR, "fig3d_res_a_removed_scatter_subjavg.png")))
