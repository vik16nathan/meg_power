# Enable this universe
options(repos = c(
  ggseg = 'https://ggseg.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))

# Install some packages
#install.packages('ggsegSchaefer')
library(ggsegSchaefer)
library(ggseg)
library(dplyr)
library(ggplot2)
library(patchwork) # for the combined Fig. 1B/1C-style panels below

# ggsegSchaefer >= 2.x exports atlases as zero-arg constructor functions
# (schaefer7_600() -> the actual ggseg_atlas list-with-$data), not
# pre-loaded data objects like older versions -- evaluate it once here so
# every later bare `schaefer7_600` reference below (all 29 of them) keeps
# working unchanged against the real atlas object instead of the function.
schaefer7_600 <- schaefer7_600()

# ---- Paths (headless: no setwd() to a machine-specific working directory) --
# outputs_dir matches process_brainstorm_data/calculate_parcellate_primary_p's
# shared OUTPUTS_DIR convention.
outputs_dir <- "/export02/data/vikramn/hbm_manuscript_code/outputs"
primary_p_dir <- file.path(outputs_dir, "primary_p", "cent_surf_fem_fwd") # parcellate_central_surface_s600.m output
fig_dir <- file.path(outputs_dir, "primary_p", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

custom_sort <- function(schaefer7_sort) {
  lh_strings <- grep("LH", schaefer7_sort, value = TRUE)
  rh_strings <- grep("RH", schaefer7_sort, value = TRUE)

  lh_strings <- lh_strings[order(sub("LH", "", lh_strings))]
  rh_strings <- rh_strings[order(sub("RH", "", rh_strings))]

  min_len <- min(length(lh_strings), length(rh_strings))

  result <- c()
  for (i in 1:min_len) {
    result <- c(result, lh_strings[i], rh_strings[i])
  }
  result <- c(result, lh_strings[(min_len + 1):length(lh_strings)], rh_strings[(min_len + 1):length(rh_strings)])
  return(result)
}

transform_s600_colnames <- function(column_names) {
  # ggsegSchaefer's atlas polygon data (flattened by geom_brain internally)
  # embeds the hemisphere right in the region string, e.g.
  # "7Networks_LH_SomMot_15" -- not "7Networks_SomMot_15" with hemi kept in
  # a separate column. Without the LH_/RH_ marker every row fails
  # geom_brain's atlas-to-data join (verified: 0/600 matched before this
  # fix, 600/600 after), which doesn't error -- it just renders every
  # parcel with NA fill (uniform grey, no legend). Confirmed against the
  # installed ggsegSchaefer/ggseg versions on 2026-08-05.
  transformed_names <- vector("character", length(column_names))
  for (i in 1:length(column_names)) {
    original_name <- column_names[i]
    hemi_letter <- substr(original_name, nchar(original_name), nchar(original_name)) # 'L' or 'R'
    hemi_code <- if (hemi_letter == 'L') 'LH' else 'RH'
    base_name <- substr(original_name, 1, nchar(original_name) - 2) # strip trailing " L"/" R"
    transformed_names[i] <- paste0('7Networks_', hemi_code, '_', base_name)
  }
  return(transformed_names)
}

replace_outliers_with_max_non_outlier <- function(vector) {
  q <- quantile(vector, probs = c(0.25, 0.75), na.rm = TRUE)   # <<< added na.rm
  iqr <- q[2] - q[1]
  lower_bound <- q[1] - 1.5 * iqr
  upper_bound <- q[2] + 1.5 * iqr

  non_outliers <- vector[!is.na(vector) & vector >= lower_bound & vector <= upper_bound]  # <<< NA-safe
  if (length(non_outliers) == 0) return(vector)                                           # <<< guard
  max_non_outlier <- max(non_outliers, na.rm = TRUE)                                      # <<< NA-safe

  vector[vector > upper_bound] <- max_non_outlier
  return(vector)
}

# path: full path to a merged primary-P + metabolism CSV with columns
# region, mean_p, cmro2, cmrglu, mean_i. Produced by build_pri_power_met_csv.py
# (this directory), which joins parcellate_central_surface_s600.m's
# s600_p_pri_2min_rest_subavg_snr3.csv / s600_i_pri_2min_rest_subavg_snr3.csv
# with a neuromaps CMRO2/CMRGlu reference CSV.
processSchaeferInput <- function(path) {
  # keep original read.csv but ensure names + strings; do NOT fetch atlas                      # <<< comment
  data <- as.data.frame(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))        # <<< added args

  ## --- minimal column-name sanitization (prevents dplyr/ggplot errors) ---                  # <<< added block
  bad <- which(is.na(names(data)) | names(data) == "")
  if (length(bad) > 0) names(data)[bad] <- paste0("V", bad)

  rgn_vec <- as.vector(data['region'])$region
  #Rename columns to align with ggseg Schaeffer 600 names
  data['region'] <- transform_s600_colnames(rgn_vec)

  mean_p_out_filt <- replace_outliers_with_max_non_outlier(as.numeric(as.matrix((data['mean_p']))))
  mean_i_out_filt <- replace_outliers_with_max_non_outlier(as.numeric(as.matrix((data['mean_i']))))

  hemi_string <- c()
  for (i in c(1:length(rgn_vec))) {
    rgn <- rgn_vec[i]
    if(substr(rgn,nchar(rgn),nchar(rgn)) == 'L') {
      hemi_string <- c(hemi_string, 'left')
    } else {
      hemi_string <- c(hemi_string, 'right')
    }
  }

  # data[,2:ncol(data)] (not 3:) -- keeps ALL 4 original raw columns
  # (mean_p, cmro2, cmrglu, mean_i) alongside the 2 derived *_out_filt ones,
  # giving exactly 8 columns matching new8 below one-for-one. The previous
  # data[,3:ncol(data)] dropped the raw 'mean_p' column entirely, which
  # silently shifted every later column's data out of alignment with its
  # name (the column labeled 'mean_p' actually held cmro2's values, 'cmro2'
  # held cmrglu's, 'cmrglu' held mean_i's, and no column was named 'mean_i'
  # at all) -- not a crash, just wrong data in every affected plot.
  data <- cbind(data['region'], hemi_string,
                mean_p_out_filt, mean_i_out_filt, data[,2:ncol(data)])

  new8 <- c('region', 'hemi', 'mean_p_out_filt', 'mean_i_out_filt',
            'mean_p', 'cmro2', 'cmrglu','mean_i')
  oldn <- colnames(data)
  len <- min(length(oldn), length(new8))
  oldn[1:len] <- new8[1:len]
  colnames(data) <- oldn

  return(data)
}


# =============================================================================
# Read in primary P data (PRIMARY ONLY -- see bottom of script for the
# secondary-P comparison, which is disabled/commented out)
# =============================================================================
s600_pri_p_rms <- processSchaeferInput(file.path(primary_p_dir, 'rms_2min_cent_pri_power_met_broad.csv'))

rest_data_schaefer <- s600_pri_p_rms
###FIX REGION NAMES###
## --- minimal join hygiene + two name patches ---
# schaefer7_600$data is a nested vertices/geom structure in this ggsegSchaefer
# version, not a flat table with its own $region column -- the flattened,
# joinable region list is what geom_brain's internal prepare_polygon_atlas()
# produces, so pull the sanity-check list from the same place.
atlas_names <- unique(ggseg:::prepare_polygon_atlas(schaefer7_600, hemi = NULL, view = NULL,
                                                      position = position_brain(), context = TRUE,
                                                      focus = NULL)$region)

# ensure clean keys
rest_data_schaefer$region <- trimws(rest_data_schaefer$region)
rest_data_schaefer$hemi   <- tolower(trimws(rest_data_schaefer$hemi))

# revert any accidental PFCl->PFCI changes (atlas uses PFCl)
rest_data_schaefer$region <- sub("PFCI_", "PFCl_", rest_data_schaefer$region, fixed = TRUE)

# fix the two labels that don't exist in your atlas build -- both are
# hemisphere-specific (verified against the real atlas region list): LH's
# atlas name is "PFC_8" (not "PFCdPFCm_8", unlike every other PFCdPFCm
# parcel and unlike RH's own "PFCdPFCm_8" which IS correct as-is); RH has no
# "Vis_45" (LH genuinely has both Vis_44 and Vis_45, so this patch must NOT
# touch LH).
rest_data_schaefer$region <- sub("^7Networks_LH_Default_PFCdPFCm_8$", "7Networks_LH_Default_PFC_8",
                                 rest_data_schaefer$region)
rest_data_schaefer$region <- sub("^7Networks_RH_Vis_45$", "7Networks_RH_Vis_44",
                                 rest_data_schaefer$region)

# sanity check: any leftovers?
setdiff(unique(rest_data_schaefer$region), atlas_names)

# =============================================================================
# PRIMARY P FIGURES
# =============================================================================

#Start with p
varname <- 'mean_p_out_filt'
ggplot(rest_data_schaefer) +
  geom_brain(data = rest_data_schaefer, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping = aes(fill = !!sym(varname))) +
  labs(fill='W') +
  ggtitle('Total FEM-MNE Primary Dipole P (Thresholded)') +
  theme_void()
ggsave(file.path(fig_dir, 'rms_2min_total_p_out_filt.png'), width = 10, height = 3.5, bg = "white")

varname <- 'mean_p'
ggplot(rest_data_schaefer) + geom_brain(data = rest_data_schaefer, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='W') + ggtitle('FEM-MNE Primary P, 120 s RMS')
ggsave(file.path(fig_dir, 'rms_2min_total_p_pri.png'), width = 10, height = 3.5, bg = "white")

varname <- 'mean_i_out_filt'
ggplot(rest_data_schaefer) + geom_brain(data = rest_data_schaefer, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='A.m') + ggtitle('Total Primary Current Dipole (Thresholded)') + theme_void()
ggsave(file.path(fig_dir, 'rms_2min_mean_i_out_filt.png'), width = 10, height = 3.5, bg = "white")

varname <- 'mean_i'
ggplot(rest_data_schaefer) + geom_brain(data = rest_data_schaefer, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='A.m') + ggtitle('FEM-MNE Primary Current Dipole, 120 s RMS')
ggsave(file.path(fig_dir, 'rms_2min_mean_i.png'), width = 10, height = 3.5, bg = "white")

varname <- 'cmrglu'
ggplot(rest_data_schaefer) + geom_brain(data = rest_data_schaefer, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='cmrglu (a.u.)') + ggtitle('CMrGlu (Neuromaps)')
ggsave(file.path(fig_dir, 'rest_s600_mean_cmrglu.png'), width = 10, height = 3.5, bg = "white")

varname <- 'cmro2'
ggplot(rest_data_schaefer) + geom_brain(data = rest_data_schaefer, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='cmro2') + ggtitle('CMrO2 (Neuromaps)')
ggsave(file.path(fig_dir, 'rest_s600_mean_cmro2.png'), width = 10, height = 3.5, bg = "white")


#Repeat for z scores
varname <- 'mean_p_out_filt_z'
rest_data_schaefer[varname] <- scale(rest_data_schaefer['mean_p_out_filt'])[1:600]
max_abs_value <- max(abs(rest_data_schaefer[varname]))
ggplot(rest_data_schaefer) + geom_brain(data = rest_data_schaefer, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled P') + ggtitle('Mean Rest P (Thresholded)') + theme_void() +
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))
ggsave(file.path(fig_dir, 'rest_s600_mean_p_out_filt_z.png'), width = 10, height = 3.5, bg = "white")

varname <- 'mean_p_z'
rest_data_schaefer[varname] <- scale(rest_data_schaefer['mean_p'])
max_abs_value <- max(abs(rest_data_schaefer[varname]))
ggplot(rest_data_schaefer) + geom_brain(data = rest_data_schaefer, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled P') + ggtitle('Mean Rest P') +  theme_void() +
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))
ggsave(file.path(fig_dir, 'rest_s600_mean_p_z.png'), width = 10, height = 3.5, bg = "white")

varname <- 'mean_i_out_filt_z'
rest_data_schaefer[varname] <- scale(rest_data_schaefer['mean_i_out_filt'])[1:600]
max_abs_value <- max(abs(rest_data_schaefer[varname]))
ggplot(rest_data_schaefer) + geom_brain(data = rest_data_schaefer, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled Thresholded I') + ggtitle('Mean Rest I (No Outliers)') + theme_void() +
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))
ggsave(file.path(fig_dir, 'rest_s600_mean_i_no_out_z.png'), width = 10, height = 3.5, bg = "white")

varname <- 'mean_i_z'
rest_data_schaefer[varname] <- scale(rest_data_schaefer['mean_i'])
max_abs_value <- max(abs(rest_data_schaefer[varname]))
ggplot(rest_data_schaefer) + geom_brain(data = rest_data_schaefer, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) + theme_void() +
  labs(fill='Scaled I') + ggtitle('Mean Rest I') +
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))
ggsave(file.path(fig_dir, 'rest_s600_mean_i_z.png'), width = 10, height = 3.5, bg = "white")

varname <- 'cmrglu_z'
rest_data_schaefer[varname] <- scale(rest_data_schaefer['cmrglu'])
max_abs_value <- max(abs(rest_data_schaefer[varname]))
ggplot(rest_data_schaefer) + geom_brain(data = rest_data_schaefer, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled CMrGlu') + ggtitle('Mean Rest CMrGlu') + theme_void() +
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))
ggsave(file.path(fig_dir, 'rest_s600_mean_cmrglu_z.png'), width = 10, height = 3.5, bg = "white")

varname <- 'cmro2_z'
rest_data_schaefer[varname] <- scale(rest_data_schaefer['cmro2'])
max_abs_value <- max(abs(rest_data_schaefer[varname]))
ggplot(rest_data_schaefer) + geom_brain(data = rest_data_schaefer, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled CMrO2') + ggtitle('Mean Rest CMrO2') + theme_void() +
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))
ggsave(file.path(fig_dir, 'rest_s600_mean_cmro2_z.png'), width = 10, height = 3.5, bg = "white")


#Plot cortical column volume (m^3) and resistance (Ohms)
#Power = squared current dipole / cortical column volume
#Power = (current (dipole/thickness))^2 * resistance
# vol/res/thick/sa CSVs are process_brainstorm_data/parcellate_cortical_column_res_vol.m's
# subject-averaged outputs, written to the shared outputs_dir root.
vol_s600 <- as.data.frame(read.csv(file.path(outputs_dir, 's600_subavg_column_vol.csv'))) ###repeat for n=4 subjects
colnames(vol_s600) <- c('region', 'vol')

res_s600 <- as.data.frame(read.csv(file.path(outputs_dir, 's600_subavg_column_res.csv')))
colnames(res_s600) <- c('region', 'res')

#NEW (Fig. 2D dendrite-derived resistance -- see meg_power/cell_res; not produced by this pipeline)
# s600_rgn_r_map_default.csv is never produced by any script in this
# pipeline (it's Fig. 2D's separate cell_res_avg.m/res_correlations.R
# analysis) -- commented out rather than loaded.
# res_new_path <- file.path(outputs_dir, 's600_rgn_r_map_default.csv')
# res_s600_new <- as.data.frame(read.csv(res_new_path))
# colnames(res_s600_new) <- c('region', 'res_new')

thick_s600 <- as.data.frame(read.csv(file.path(outputs_dir, 's600_subavg_column_thick.csv')))
colnames(thick_s600) <- c('region', 'thick')

thick_s600['thick_no_out'] <- replace_outliers_with_max_non_outlier(thick_s600['thick']$thick)


sa_s600 <- as.data.frame(read.csv(file.path(outputs_dir, 's600_subavg_column_sa.csv')))
colnames(sa_s600) <- c('region', 'sa')

vol_res_s600 <- cbind(vol_s600, res_s600[,2],
                      thick_s600[,2:3], sa_s600[,2])

colnames(vol_res_s600) <- c('region', 'vol', 'res', 'thick', 'thick_no_out','sa')

rgn_vec <- as.vector(vol_res_s600['region'])$region
vol_res_s600['region'] <- transform_s600_colnames(rgn_vec)

#Visualize resistance, thickness, and volume per cortical column
varname <- 'thick_no_out'
vol_res_s600[,'thick_no_out'] <- vol_res_s600[,'thick_no_out']*1000
ggplot(vol_res_s600) + geom_brain(data = vol_res_s600, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) + theme_void() +
  labs(fill='mm') + ggtitle('Mean Cortical Thickness')
ggsave(file.path(fig_dir, 'rest_s600_mean_cort_thick_no_out.png'), width = 10, height = 3.5, bg = "white")

varname <- 'res' #uniform grey-matter resistivity, R=rho*L/A (Fig. 1A/1B(iii))
ggplot(vol_res_s600) + geom_brain(data = vol_res_s600, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='Ohms') + ggtitle('Mean Cortical Columnar Resistance')
ggsave(file.path(fig_dir, 'rest_s600_mean_res.png'), width = 10, height = 3.5, bg = "white")

# varname <- 'res_new' #dendrite/neuron-density-derived resistance (Fig. 2D)
# ggplot(vol_res_s600) + geom_brain(data = vol_res_s600, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
#   labs(fill='Ohms') + ggtitle('Mean Neuron-Derived Resistance')
# ggsave(file.path(fig_dir, 'rest_s600_mean_res_new.png'), bg = "white")
#
# varname <- 'res_new_no_out' #single-cell, number of neurons
# vol_res_s600['res_new_no_out'] <- replace_outliers_with_max_non_outlier(vol_res_s600['res_new']$res_new)
#
# ggplot(vol_res_s600) + geom_brain(data = vol_res_s600, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
#   labs(fill='Ohms') + ggtitle('Mean Neuron-Derived Resistance (Thresholded)')
# ggsave(file.path(fig_dir, 'rest_s600_mean_res_new_no_out.png'), bg = "white")

varname <- 'res_no_out'
vol_res_s600['res_no_out'] <- replace_outliers_with_max_non_outlier(vol_res_s600['res']$res)
ggplot(vol_res_s600) + geom_brain(data = vol_res_s600, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='Ohms') + ggtitle('Mean Cortical Column Resistance (Thresholded)') + theme_void()
ggsave(file.path(fig_dir, 'rest_s600_mean_r_out_filt.png'), width = 10, height = 3.5, bg = "white")

varname <- 'vol'
ggplot(vol_res_s600) + geom_brain(data = vol_res_s600, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname)))  +
  labs(fill='m^3') + ggtitle('Mean Cortical Column Volume')
ggsave(file.path(fig_dir, 'rest_s600_mean_v.png'), width = 10, height = 3.5, bg = "white")

varname <- 'sa'
ggplot(vol_res_s600) + geom_brain(data = vol_res_s600, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='m^2') + ggtitle('Mean Cortical Column Support Area')
ggsave(file.path(fig_dir, 'rest_s600_mean_sa.png'), width = 10, height = 3.5, bg = "white")

#Repeat for z scores
varname <- 'thick_z'
vol_res_s600[varname] <- scale(vol_res_s600['thick'])[1:600]
max_abs_value <- max(abs(vol_res_s600[varname]))
ggplot(vol_res_s600) + geom_brain(data = vol_res_s600, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled Cortical Thickness') + ggtitle('Mean Thickness (Standardized)') + theme_void()+
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))
ggsave(file.path(fig_dir, 'rest_s600_mean_cort_thick_z.png'), width = 10, height = 3.5, bg = "white")

varname <- 'thick_no_out_z'
vol_res_s600[varname] <- scale(vol_res_s600['thick_no_out'])[1:600]
max_abs_value <- max(abs(vol_res_s600[varname]))
ggplot(vol_res_s600) + geom_brain(data = vol_res_s600, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled Cortical Thickness') + ggtitle('Mean thickness (Standardized)') + theme_void()+
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))
ggsave(file.path(fig_dir, 'rest_s600_mean_cort_thick_no_out_z.png'), width = 10, height = 3.5, bg = "white")

varname <- 'res_z'
vol_res_s600[varname] <- scale(vol_res_s600['res'])[1:600]
max_abs_value <- max(abs(vol_res_s600[varname]))
ggplot(vol_res_s600) + geom_brain(data = vol_res_s600, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled Resistance') + ggtitle('Mean Resistance (Standardized)') + theme_void() +
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))
ggsave(file.path(fig_dir, 'rest_s600_mean_cort_r_z.png'), width = 10, height = 3.5, bg = "white")

varname <- 'res_no_out_z'
vol_res_s600[varname] <- scale(vol_res_s600['res_no_out'])[1:600]
max_abs_value <- max(abs(vol_res_s600[varname]))
ggplot(vol_res_s600) + geom_brain(data = vol_res_s600, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled Resistance') + ggtitle('Mean Resistance (Standardized)') + theme_void() +
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))
ggsave(file.path(fig_dir, 'rest_s600_mean_cort_r_no_out_z.png'), width = 10, height = 3.5, bg = "white")

varname <- 'vol_z'
vol_res_s600[varname] <- scale(vol_res_s600['vol'])[1:600]
max_abs_value <- max(abs(vol_res_s600[varname]))
ggplot(vol_res_s600) + geom_brain(data = vol_res_s600, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled Column Volume') + ggtitle('Mean Column Volume (Standardized)') + theme_void() +
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))
ggsave(file.path(fig_dir, 'rest_s600_mean_cort_v_z.png'), width = 10, height = 3.5, bg = "white")

varname <- 'sa_z'
vol_res_s600[varname] <- scale(vol_res_s600['sa'])[1:600]
max_abs_value <- max(abs(vol_res_s600[varname]))
ggplot(vol_res_s600) + geom_brain(data = vol_res_s600, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled Column Area') + ggtitle('Mean Column Area (Standardized)') + theme_void() +
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))
ggsave(file.path(fig_dir, 'rest_s600_mean_cort_sa_z.png'), width = 10, height = 3.5, bg = "white")


#Correlate power with 1/volume - all residual variation is due to currents
vol <- vol_res_s600['vol']$vol
res <- vol_res_s600['res']$res

power <- as.vector(rest_data_schaefer['mean_p'])$mean_p
curr_dip <- as.vector(rest_data_schaefer['mean_i'])$mean_i

png(file.path(fig_dir, 'rest_s600_power_vs_vol.png'))
plot(x=vol, y=power)
dev.off()

#Correlate dipole current with volume
png(file.path(fig_dir, 'rest_s600_curr_dip_vs_vol.png'))
plot(x=vol, y=curr_dip, ylab='Avg. Current Dipole in Parcel (A.m)',
     xlab='Avg. Cortical Column Volume (m^3)', main='Schaefer 600 Current Dipole vs. Column Volume')
dev.off()

#Correlate dipole current with resistance
png(file.path(fig_dir, 'rest_s600_curr_dip_vs_res.png'))
plot(x=res, y=curr_dip, xlab='Avg. Column Resistance in Parcel (Ohms)',
     ylab='Avg. Current Dipole in Parcel (A.m)', main='Schaefer 600 Current Dipole vs. Column Resistances')
dev.off()


# =============================================================================
# Combined multi-panel figures resembling Fig. 1B/1C of Nathan et al.,
# "Lower Bound Estimates for Electrophysiological Power Dissipation in Human
# Gray Matter". Each row is one brain-map panel (LH lateral / LH medial /
# RH medial / RH lateral, via ggseg's default hemi x view layout), stacked
# vertically with patchwork, one independently-scaled legend per row -- same
# arrangement as the preprint figure.
# =============================================================================

brain_row <- function(data, varname, fill_label, title, palette = "Blues") {
  ggplot(data) +
    geom_brain(data = data, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping = aes(fill = !!sym(varname))) +
    labs(fill = fill_label) + ggtitle(title) + theme_void() +
    scale_fill_distiller(palette = palette, direction = 1)
}

diverging_row <- function(data, varname, fill_label, title, limit) {
  ggplot(data) +
    geom_brain(data = data, atlas = schaefer7_600, view = c('lateral', 'medial'), mapping = aes(fill = !!sym(varname))) +
    labs(fill = fill_label) + ggtitle(title) + theme_void() +
    scale_fill_gradientn(colours = c('blue', 'white', 'red'), limits = c(-limit, limit))
}

# Fig. 1B: input features and primary dipole P
fig1b_i   <- brain_row(rest_data_schaefer, 'mean_i_out_filt', 'A.m',
                        '(i) Average current dipole moment; resting state (120 s)')
fig1b_ii  <- brain_row(vol_res_s600, 'thick_no_out', 'mm',
                        '(ii) CAT12 estimated cortical thickness')
fig1b_iii <- brain_row(vol_res_s600, 'res_no_out', 'Ohms',
                        '(iii) Cortical column resistance (R)')
fig1b_iv  <- brain_row(rest_data_schaefer, 'mean_p_out_filt', 'W',
                        '(iv) Total primary dipole power (averaged over 120 s)')

fig1b <- fig1b_i / fig1b_ii / fig1b_iii / fig1b_iv +
  plot_annotation(title = 'B  Input features and primary dipole P')
ggsave(file.path(fig_dir, 'fig1b_primary_p_inputs.png'), plot = fig1b, width = 10, height = 14, bg = "white")

# Fig. 1C: neuromaps glucose and oxygen metabolism (shared diverging scale,
# matching the preprint's shared -5/+5 z-score range)
metab_limit <- max(abs(rest_data_schaefer$cmrglu_z), abs(rest_data_schaefer$cmro2_z), na.rm = TRUE)
fig1c_i  <- diverging_row(rest_data_schaefer, 'cmrglu_z', 'Scaled CMrGlu',
                           '(i) Z-scored glucose metabolism', metab_limit)
fig1c_ii <- diverging_row(rest_data_schaefer, 'cmro2_z', 'Scaled CMrO2',
                           '(ii) Standardized oxygen metabolism', metab_limit)

fig1c <- fig1c_i / fig1c_ii +
  plot_annotation(title = 'C  neuromaps glucose and oxygen metabolism')
ggsave(file.path(fig_dir, 'fig1c_metabolism.png'), plot = fig1c, width = 10, height = 7, bg = "white")


# =============================================================================
# Fig. 3D: primary vs. secondary dipole P, parcellated (Schaefer-600),
# averaged over the same 4 exemplar subjects used throughout Fig. 3 (sub-0002/
# 8/19/22 -- secondary-current FEM power is only computed for these 4, not
# the full cohort; see secondary_p/README.md). Inputs are
# calculate_parcellate_primary_p/parcellate_secondary_dipole_p.py's output:
# that script parcellates secondary_p/eeg_dipole_power_computation.py's
# per-vertex fem_p_rms.csv into the same Schaefer-600 7-network regions
# parcellate_central_surface_s600.m uses for primary P (region value = sum of
# the per-vertex quantity over that region's vertices), then averages both
# primary and secondary P over exactly these 4 subjects -- a matched
# comparison, not primary's full-cohort subavg vs. secondary's n=4.
# =============================================================================
sec_dir <- file.path(outputs_dir, "secondary_p", "cent_surf_fem_fwd") # parcellate_secondary_dipole_p.py output

processSimpleRegionValueInput <- function(path, value_colname) {
  data <- as.data.frame(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  data$region <- transform_s600_colnames(data$region)

  # same two subject-anatomy-vs-atlas-template name mismatches
  # processSchaeferInput patches above (verified against the real atlas
  # region list, not specific to primary vs. secondary P)
  data$region <- sub("PFCI_", "PFCl_", data$region, fixed = TRUE)
  data$region <- sub("^7Networks_LH_Default_PFCdPFCm_8$", "7Networks_LH_Default_PFC_8", data$region)
  data$region <- sub("^7Networks_RH_Vis_45$", "7Networks_RH_Vis_44", data$region)

  data[[paste0(value_colname, "_out_filt")]] <- replace_outliers_with_max_non_outlier(data$value)
  colnames(data)[colnames(data) == "value"] <- value_colname
  return(data)
}

s600_pri_p_exemplars <- processSimpleRegionValueInput(
  file.path(sec_dir, sprintf("s600_p_pri_2min_rest_subavg_snr%d_exemplars.csv", 3)), "mean_p")
s600_sec_p_exemplars <- processSimpleRegionValueInput(
  file.path(sec_dir, sprintf("s600_p_sec_dipole_2min_rest_subavg_snr%d_exemplars.csv", 3)), "mean_p_sec")

# sanity check: any leftovers not in the real atlas region list? (same check
# rest_data_schaefer gets above)
setdiff(unique(s600_pri_p_exemplars$region), atlas_names)
setdiff(unique(s600_sec_p_exemplars$region), atlas_names)

# Same "Blues" single-hue sequential colorscheme as Fig. 1B's brain_row()
# (scale_fill_distiller, direction=1) -- each panel independently scaled
# (primary and secondary P differ by roughly an order of magnitude), matching
# how Fig. 1B's rows are each independently scaled too.
fig3d_primary <- brain_row(s600_pri_p_exemplars, "mean_p_out_filt", "W",
                            "Total FEM-MNE Primary Dipole P (Thresholded)")
fig3d_secondary <- brain_row(s600_sec_p_exemplars, "mean_p_sec_out_filt", "W",
                              "FEM-MNE Secondary Dipole P (Thresholded)")

fig3d <- fig3d_primary / fig3d_secondary +
  plot_annotation(title = "D  FEM primary dipole P and secondary dipole P show spatial differences.")
ggsave(file.path(fig_dir, "fig3d_primary_vs_secondary_dipole_p.png"), plot = fig3d,
       width = 10, height = 7, bg = "white")
