# Supplementary figure: effect of the MNE regularization parameter
# (SnrFixed, swept by rms_idip_rest.m: 1, 3, 5) on the reconstructed primary
# current dipole moment (I_dip). Two panels:
#   A) Whole-brain summed I_dip, averaged across subjects, one bar per
#      SnrFixed value (+/- SEM across subjects).
#   B) Schaeffer-600 parcellated I_dip spatial maps, one per SnrFixed value,
#      on a shared color scale so they're visually comparable.
#
# Inputs (produced by parcellate_central_surface_s600.m, this directory):
# - s600_<subject>_snr<N>_i_pri_2min_rest.csv (per-subject, per-region I_dip) -> panel A
# - s600_i_pri_2min_rest_subavg_snr<N>.csv (subject-averaged, per-region I_dip) -> panel B

library(ggsegSchaefer)
library(ggseg)
library(dplyr)
library(ggplot2)
library(patchwork)

# ggsegSchaefer >= 2.x exports atlases as zero-arg constructor functions
# (schaefer7_600() -> the actual ggseg_atlas list-with-$data), not
# pre-loaded data objects -- see schaeffer_rest_power_map.R's matching note.
schaefer7_600 <- schaefer7_600()

outputs_dir <- "/export02/data/vikramn/hbm_manuscript_code/outputs"
primary_p_dir <- file.path(outputs_dir, "primary_p", "cent_surf_fem_fwd") # parcellate_central_surface_s600.m output
fig_dir <- file.path(outputs_dir, "primary_p", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

SnrValues <- c(1, 3, 5) # must match the sweep in rms_idip_rest.m

# Confirmed outlier (see snr_sensitivity_idip_rms.R's matching comment) --
# excluded explicitly by subject ID rather than relying on its per-subject
# CSVs being absent, since stale sub-0011 files from before the exclusion
# may still be sitting in primary_p_dir.
OUTLIER_SUBJECT <- 'sub-0011'

# Same region-name convention/patches as schaeffer_rest_power_map.R (see its
# matching comment): ggsegSchaefer's atlas embeds hemisphere in the region
# string itself (e.g. "7Networks_LH_SomMot_15"), not a separate column --
# without the LH_/RH_ marker every row silently fails geom_brain's
# atlas-to-data join (renders as uniform-grey/no-legend, no error).
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
  # Both patches are hemisphere-specific -- see schaeffer_rest_power_map.R's
  # matching comment (verified against the real atlas region list).
  region <- sub("^7Networks_LH_Default_PFCdPFCm_8$", "7Networks_LH_Default_PFC_8", region)
  region <- sub("^7Networks_RH_Vis_45$", "7Networks_RH_Vis_44", region)
  region
}

# =============================================================================
# Panel A: whole-brain summed I_dip per subject per SnrFixed
# =============================================================================
whole_brain_rows <- list()
for (snr in SnrValues) {
  subj_files <- list.files(primary_p_dir,
                            pattern = sprintf('^s600_sub-[0-9]+_snr%d_i_pri_2min_rest\\.csv$', snr),
                            full.names = TRUE)
  subj_files <- subj_files[!grepl(OUTLIER_SUBJECT, basename(subj_files), fixed = TRUE)]
  for (f in subj_files) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    whole_brain_rows[[length(whole_brain_rows) + 1]] <- data.frame(
      file = basename(f), snr = snr, total_i = sum(d$value, na.rm = TRUE))
  }
}
if (length(whole_brain_rows) == 0) {
  stop("No s600_sub-*_snr<N>_i_pri_2min_rest.csv files found in ", primary_p_dir,
       " -- run parcellate_central_surface_s600.m first.")
}
whole_brain_df <- do.call(rbind, whole_brain_rows)

summary_df <- whole_brain_df %>%
  group_by(snr) %>%
  summarise(mean_i = mean(total_i), sd_i = sd(total_i), n = n(), .groups = "drop") %>%
  mutate(se_i = sd_i / sqrt(n))

barplot_panel <- ggplot(summary_df, aes(x = factor(snr), y = mean_i)) +
  geom_col(fill = "#4C72B0", width = 0.6) +
  geom_errorbar(aes(ymin = mean_i - se_i, ymax = mean_i + se_i), width = 0.15) +
  labs(x = "SnrFixed (MNE regularization parameter)",
       y = "Whole-brain summed I_dip (A.m)",
       title = sprintf("A  Average total primary current dipole moment (n=%d subjects)", max(summary_df$n))) +
  theme_minimal(base_size = 13)

ggsave(file.path(fig_dir, 'snr_sensitivity_barplot.png'), plot = barplot_panel, width = 6, height = 5, bg = "white")

# =============================================================================
# Panel B: parcellated I_dip maps, one per SnrFixed, shared color scale
# =============================================================================
region_maps <- list()
for (snr in SnrValues) {
  f <- file.path(primary_p_dir, sprintf('s600_i_pri_2min_rest_subavg_snr%d.csv', snr))
  if (!file.exists(f)) {
    warning("Missing ", f, " -- skipping SnrFixed=", snr, " in the parcellated maps.")
    next
  }
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(d) <- c('region', 'mean_i')
  d$region <- transform_s600_colnames(d$region)
  d$region <- fix_region_names(d$region)
  d$snr <- snr
  region_maps[[length(region_maps) + 1]] <- d
}
if (length(region_maps) == 0) {
  stop("No s600_i_pri_2min_rest_subavg_snr<N>.csv files found in ", primary_p_dir,
       " -- run parcellate_central_surface_s600.m first.")
}
region_maps_df <- do.call(rbind, region_maps)
shared_limit <- max(region_maps_df$mean_i, na.rm = TRUE)
shared_min <- min(region_maps_df$mean_i, na.rm = TRUE)

# SnrFixed=5's max is ~4x SnrFixed=1's (see rms_idip_rest.m's SnrFixed effect
# on magnitude), so a shared LINEAR scale (limits = c(0, shared_limit))
# collapses SnrFixed=1's whole map into near-white -- visually indistinguishable,
# looking like only SnrFixed=3/5 "have data" even though every region has a
# real nonzero value. A shared log10 scale keeps all three rows' spatial
# texture visible while still using one consistent scale for comparison
# across SnrFixed (verified no zero/NA values in any of the three CSVs, so
# log10 is safe here).
brain_row <- function(data, title) {
  ggplot(data) +
    geom_brain(data = data, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping = aes(fill = mean_i)) +
    labs(fill = 'A.m') + ggtitle(title) + theme_void() +
    scale_fill_distiller(palette = "Blues", direction = 1, trans = "log10", limits = c(shared_min, shared_limit))
}

brain_plots <- lapply(sort(unique(region_maps_df$snr)), function(snr) {
  brain_row(region_maps_df[region_maps_df$snr == snr, ], sprintf("SnrFixed = %d", snr))
})

maps_panel <- wrap_plots(brain_plots, ncol = 1) +
  plot_annotation(title = "B  Parcellated primary current dipole moment (I_dip)")
# geom_brain enforces a fixed x/y aspect ratio (real brain shapes, not
# stretched), so each row only actually needs ~1.8in of height at width=10in
# -- the old `4 * length(brain_plots)` requested far more canvas than that,
# and since patchwork can't stretch fixed-aspect content to fill it, the
# excess just became blank margin above and below the whole 3-row block
# (measured: ~29% blank top, ~29% blank bottom of the old 12in-tall canvas).
ggsave(file.path(fig_dir, 'snr_sensitivity_maps.png'), plot = maps_panel,
       width = 10, height = 2 * length(brain_plots), bg = "white")

# =============================================================================
# Combined supplementary figure (barplot + stacked parcellated maps)
# =============================================================================
combined <- wrap_plots(c(list(barplot_panel), brain_plots), ncol = 1) +
  plot_layout(heights = c(1, rep(1.3, length(brain_plots)))) +
  plot_annotation(title = "Effect of the MNE regularization parameter (SnrFixed) on primary current dipole moment")

ggsave(file.path(fig_dir, 'supp_fig_snr_sensitivity.png'), plot = combined,
       width = 10, height = 5 + 2 * length(brain_plots), bg = "white")
