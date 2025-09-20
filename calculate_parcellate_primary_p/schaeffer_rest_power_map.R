# Enable this universe
options(repos = c(
  ggseg = 'https://ggseg.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))

# Install some packages
#install.packages('ggsegSchaefer')
library(ggsegSchaefer)
library(ggseg3d)
library(ggseg)
library(dplyr)
library(ggplot2)

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

#Load in the data we're working with
setwd("C:/Users/vik16/Documents/Baillet Lab Manuscript/Baillet Lab/Baillet Lab")

transform_s600_colnames <- function(column_names) {
  transformed_names <- vector("character", length(column_names))
  for (i in 1:length(column_names)) {
    original_name <- column_names[i]
    final_name <- paste0('7Networks_', substr(original_name, 1, nchar(original_name) - 2))
    transformed_names[i] <- final_name
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
  
  data <- cbind(data['region'], hemi_string,
                mean_p_out_filt, mean_i_out_filt, data[,3:ncol(data)])
  
  ## assign names to the FIRST 8 columns only, keep the rest untouched                         # <<< changed
  new8 <- c('region', 'hemi', 'mean_p_out_filt', 'mean_i_out_filt',
            'mean_p', 'cmro2', 'cmrglu','mean_i')
  oldn <- colnames(data)
  len <- min(length(oldn), length(new8))
  oldn[1:len] <- new8[1:len]
  colnames(data) <- oldn
  
  return(data)
}



#Read in data
s600_sec_p_rms <- processSchaeferInput('./rest_p_rms/sec/rms_2min_cent_sec_power_met_broad.csv')
s600_pri_p_rms <- processSchaeferInput('./rest_p_rms/pri/rms_2min_cent_pri_power_met_broad.csv')

rest_data_schaefer <- s600_pri_p_rms
###FIX REGION NAMES###
## --- minimal join hygiene + two name patches ---
atlas_names <- unique(schaefer7_600$data$region)

# ensure clean keys
rest_data_schaefer$region <- trimws(rest_data_schaefer$region)
rest_data_schaefer$hemi   <- tolower(trimws(rest_data_schaefer$hemi))

# revert any accidental PFCl->PFCI changes (atlas uses PFCl)
rest_data_schaefer$region <- sub("PFCI_", "PFCl_", rest_data_schaefer$region, fixed = TRUE)

# fix the two labels that don't exist in your atlas build
rest_data_schaefer$region <- sub("^7Networks_Default_PFCdPFCm_8$", "7Networks_Default_PFC_8",
                                 rest_data_schaefer$region)
rest_data_schaefer$region <- sub("^7Networks_Vis_45$", "7Networks_Vis_44",
                                 rest_data_schaefer$region)

# sanity check: any leftovers?
setdiff(unique(rest_data_schaefer$region), atlas_names)

#Start with p
varname <- 'mean_p_out_filt'
ggplot(rest_data_schaefer) +
  geom_brain(atlas = schaefer7_600, mapping = aes(fill = !!sym(varname))) +
  labs(fill='W') +
  ggtitle('Total FEM-MNE Primary Dipole P (Thresholded)') +
  theme_void()

ggsave('./rms_2min_sec_p_out_filt.png')

# Example input
#schaefer7_sort <- sort(unique(schaefer7_600$data$region))
#sorted_result <- unique(na.omit(custom_sort(schaefer7_sort)))
#write.csv(sorted_result, 's600_ggseg_colnames.csv',row.names=FALSE)

#Read in data
#schaefer7_600_names <- as.data.frame(read.csv('s600_ggseg_colnames.csv'))
s600_sec_p_rms <- processSchaeferInput('./rest_p_rms/sec/rms_2min_cent_sec_power_met_broad.csv')
s600_pri_p_rms <- processSchaeferInput('./rest_p_rms/pri/rms_2min_cent_pri_power_met_broad.csv')
#s600_100s_subavg <- processSchaeferInput('./rest_power_redo/omega_s600_power_met_broad.csv')

rest_data_schaefer <- s600_pri_p_rms

#Start with p
varname <- 'mean_p_out_filt'
ggplot(rest_data_schaefer) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='W') + ggtitle('Total FEM-MNE Primary Dipole P (Thresholded)') + theme_void()
#ggsave('./rest_s600_total_p_out_filt.png')
ggsave('./rms_2min_sec_p_out_filt.png')

###create new variable for sec - pri
rest_data_schaefer[,"pri_min_sec"] <- s600_pri_p_rms$mean_p_out_filt - s600_sec_p_rms$mean_p_out_filt
varname <- 'pri_min_sec'
ggplot(rest_data_schaefer) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='W') + ggtitle('Pri. - Sec. Dipole P (Thresholded)') + theme_void()
#ggsave('./rest_s600_total_p_out_filt.png')
ggsave('./rms_2min_pri_min_sec_out_filt.png')

varname <- 'mean_p'
ggplot(rest_data_schaefer) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='W') + ggtitle('FEM-MNE Primary P, 120 s RMS')
ggsave('./rms_2min_total_p_pri.png')

varname <- 'mean_i_out_filt'
ggplot(rest_data_schaefer) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='A.m') + ggtitle('Total Primary Current Dipole (Thresholded)') + theme_void()
ggsave('./rms_2min_mean_i_out_filt.png')

varname <- 'mean_i'
ggplot(rest_data_schaefer) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='A.m') + ggtitle('FEM-MNE Primary Current Dipole, 120 s RMS')
ggsave('./rms_2min_mean_i.png')

varname <- 'cmrglu'
ggplot(rest_data_schaefer) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='cmrglu (a.u.)') + ggtitle('CMrGlu (Neuromaps)')
ggsave('./rest_s600_mean_cmrglu.png')

varname <- 'cmro2'
ggplot(rest_data_schaefer) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='cmro2') + ggtitle('CMrO2 (Neuromaps)')
ggsave('./rest_s600_mean_cmro2.png')


#Repeat for z scores
varname <- 'mean_p_out_filt_z'
rest_data_schaefer[varname] <- scale(rest_data_schaefer['mean_p_out_filt'])[1:600]
max_abs_value <- max(abs(rest_data_schaefer[varname]))
ggplot(rest_data_schaefer) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled P') + ggtitle('Mean Rest P (Thresholded)') + theme_void() + 
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))
                       
ggsave('./rest_s600_mean_p_out_filt_z.png')

varname <- 'mean_p_z'
rest_data_schaefer[varname] <- scale(rest_data_schaefer['mean_p'])
max_abs_value <- max(abs(rest_data_schaefer[varname]))
ggplot(rest_data_schaefer) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled P') + ggtitle('Mean Rest P') +  theme_void() +
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))

ggsave('./rest_s600_mean_p_z.png')

varname <- 'mean_i_out_filt_z'
rest_data_schaefer[varname] <- scale(rest_data_schaefer['mean_i_out_filt'])[1:600]
max_abs_value <- max(abs(rest_data_schaefer[varname]))
ggplot(rest_data_schaefer) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled Thresholded I') + ggtitle('Mean Rest I (No Outliers)') + theme_void() + 
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))

ggsave('./rest_s600_mean_i_no_out_z.png')

varname <- 'mean_i_z'
rest_data_schaefer[varname] <- scale(rest_data_schaefer['mean_i'])
max_abs_value <- max(abs(rest_data_schaefer[varname]))
ggplot(rest_data_schaefer) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) + theme_void() +
  labs(fill='Scaled I') + ggtitle('Mean Rest I') + 
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))

ggsave('./rest_s600_mean_i_z.png')

varname <- 'cmrglu_z'
rest_data_schaefer[varname] <- scale(rest_data_schaefer['cmrglu'])
max_abs_value <- max(abs(rest_data_schaefer[varname]))
ggplot(rest_data_schaefer) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled CMrGlu') + ggtitle('Mean Rest CMrGlu') + theme_void() +
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))

ggsave('./rest_s600_mean_cmrglu_z.png')

varname <- 'cmro2_z'
rest_data_schaefer[varname] <- scale(rest_data_schaefer['cmro2'])
max_abs_value <- max(abs(rest_data_schaefer[varname]))
ggplot(rest_data_schaefer) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled CMrO2') + ggtitle('Mean Rest CMrO2') + theme_void() +
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))

ggsave('./rest_s600_mean_cmro2_z.png')


#Plot cortical column volume (m^3) and resistance (Ohms)
#Power = squared current dipole / cortical column volume
#Power = (current (dipole/thickness))^2 * resistance

vol_s600 <- as.data.frame(read.csv('./s600_subavg_column_vol.csv')) ###repeat for n=4 subjects
colnames(vol_s600) <- c('region', 'vol')

res_s600 <- as.data.frame(read.csv('./s600_subavg_column_res.csv'))
colnames(res_s600) <- c('region', 'res')

#NEW
res_s600_new <- as.data.frame(read.csv('./s600_rgn_r_map_default.csv'))
colnames(res_s600_new) <- c('region', 'res_new')

thick_s600 <- as.data.frame(read.csv('./s600_subavg_column_thick.csv'))
colnames(thick_s600) <- c('region', 'thick')

thick_s600['thick_no_out'] <- replace_outliers_with_max_non_outlier(thick_s600['thick']$thick)


sa_s600 <- as.data.frame(read.csv('./s600_subavg_column_sa.csv'))
colnames(sa_s600) <- c('region', 'sa')

vol_res_s600 <- cbind(vol_s600, res_s600[,2],res_s600_new[,2],
                      thick_s600[,2:3], sa_s600[,2])

colnames(vol_res_s600) <- c('region', 'vol', 'res', 'res_new', 'thick', 'thick_no_out','sa')

rgn_vec <- as.vector(vol_res_s600['region'])$region
vol_res_s600['region'] <- transform_s600_colnames(rgn_vec)

#Visualize resistance, thickness, and volume per cortical column
varname <- 'thick_no_out'
vol_res_s600[,'thick_no_out'] <- vol_res_s600[,'thick_no_out']*1000
ggplot(vol_res_s600) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) + theme_void() +
  labs(fill='mm') + ggtitle('Mean Cortical Thickness')
ggsave('./rest_s600_mean_cort_thick_no_out.png')

varname <- 'res' #single-cell, number of neurons
ggplot(vol_res_s600) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='Ohms') + ggtitle('Mean Cortical Columnar Resistance')
ggsave('./rest_s600_mean_res_new.png')

varname <- 'res_new' #single-cell, number of neurons
ggplot(vol_res_s600) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='Ohms') + ggtitle('Mean Neuron-Derived Resistance')
ggsave('./rest_s600_mean_res_new.png')

varname <- 'res_new_no_out' #single-cell, number of neurons
vol_res_s600['res_new_no_out'] <- replace_outliers_with_max_non_outlier(vol_res_s600['res_new']$res_new)

ggplot(vol_res_s600) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='Ohms') + ggtitle('Mean Neuron-Derived Resistance (Thresholded)')
ggsave('./rest_s600_mean_res_new_no_out.png')

varname <- 'res_no_out'
vol_res_s600['res_no_out'] <- replace_outliers_with_max_non_outlier(vol_res_s600['res']$res)
ggplot(vol_res_s600) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='Ohms') + ggtitle('Mean Cortical Column Resistance (Thresholded)') + theme_void()
ggsave('./rest_s600_mean_r_out_filt.png')

varname <- 'vol'
ggplot(vol_res_s600) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname)))  +
  labs(fill='m^3') + ggtitle('Mean Cortical Column Volume')
ggsave('./rest_s600_mean_v.png')

varname <- 'sa'
ggplot(vol_res_s600) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='m^2') + ggtitle('Mean Cortical Column Support Area')
ggsave('./rest_s600_mean_sa.png')

#Repeat for z scores
varname <- 'thick_z'
vol_res_s600[varname] <- scale(vol_res_s600['thick'])[1:600]
max_abs_value <- max(abs(vol_res_s600[varname]))
ggplot(vol_res_s600) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled Cortical Thickness') + ggtitle('Mean Thickness (Standardized)') + theme_void()+
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value)) 

ggsave('./rest_s600_mean_cort_thick_z.png')

varname <- 'thick_no_out_z'
vol_res_s600[varname] <- scale(vol_res_s600['thick_no_out'])[1:600]
max_abs_value <- max(abs(vol_res_s600[varname]))
ggplot(vol_res_s600) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled Cortical Thickness') + ggtitle('Mean thickness (Standardized)') + theme_void()+
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))  
  
  ggsave('./rest_s600_mean_cort_thick_no_out_z.png')

varname <- 'res_z'
vol_res_s600[varname] <- scale(vol_res_s600['res'])[1:600]
max_abs_value <- max(abs(vol_res_s600[varname]))
ggplot(vol_res_s600) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled Resistance') + ggtitle('Mean Resistance (Standardized)') + theme_void() +
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))

ggsave('./rest_s600_mean_cort_r_z.png')

varname <- 'res_no_out_z'
vol_res_s600[varname] <- scale(vol_res_s600['res_no_out'])[1:600]
max_abs_value <- max(abs(vol_res_s600[varname]))
ggplot(vol_res_s600) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled Resistance') + ggtitle('Mean Resistance (Standardized)') + theme_void() +
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))

ggsave('./rest_s600_mean_cort_r_no_out_z.png')

varname <- 'vol_z'
vol_res_s600[varname] <- scale(vol_res_s600['vol'])[1:600]
max_abs_value <- max(abs(vol_res_s600[varname]))
ggplot(vol_res_s600) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled Column Volume') + ggtitle('Mean Column Volume (Standardized)') + theme_void() + 
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))

ggsave('./rest_s600_mean_cort_v_z.png')

varname <- 'sa_z'
vol_res_s600[varname] <- scale(vol_res_s600['sa'])[1:600]
max_abs_value <- max(abs(vol_res_s600[varname]))
ggplot(vol_res_s600) + geom_brain(atlas=schaefer7_600, mapping=aes(fill=!!sym(varname))) +
  labs(fill='Scaled Column Area') + ggtitle('Mean Column Area (Standardized)') + theme_void() +
  scale_fill_gradientn(colours=c('blue','white','red'), limits=c(-max_abs_value, max_abs_value))

ggsave('./rest_s600_mean_cort_sa_z.png')


#Correlate power with 1/volume - all residual variation is due to currents
vol <- vol_res_s600['vol']$vol
res <- vol_res_s600['res']$res

power <- as.vector(rest_data_schaefer['mean_p'])$mean_p
curr_dip <- as.vector(rest_data_schaefer['mean_i'])$mean_i

plot(x=vol, y=power)
#Correlate dipole current with volume 
plot(x=vol, y=curr_dip, ylab='Avg. Current Dipole in Parcel (A.m)',
     xlab='Avg. Cortical Column Volume (m^3)', main='Schaefer 600 Current Dipole vs. Column Volume')

#Correlate dipole current with resistance 
plot(x=res, y=curr_dip, xlab='Avg. Column Resistance in Parcel (Ohms)',
     ylab='Avg. Current Dipole in Parcel (A.m)', main='Schaefer 600 Current Dipole vs. Column Resistances')

