## Pelvis Statistical Appearance Model
# Last Updated: 2023-01-20

# =============================================================================

# Import functions
source("BoneDensityFunctions.R")


# =============================================================================

# info
## subject info
template_subject <- "PBQ002"
subjects <- c("PBQ005", "PBQ007", "PBQ008", "PBQ009", "PBQ011", "PBQ012", 
              "PBQ013", "PBQ014", "PBQ015", "PBQ017", "PBQ018", "PBQ019", 
              "PBQ020", "PBQ021", "PBQ023", "PBQ024", "PBQ027", "PBQ028", 
              "PBQ029", "PBQ031", "PBQ032", "PBQ052", "PBQ054", "PBQ055", 
              "PBQ056", "PBQ059", "PBQ060", #"PBQ063", "PBQ064", 
              "PBQ065", 
              "PBQ066", "PBQ067", "PBQ068", "PBQ070", "PBQ071", "PBQ072", 
              "PBQ073", "PBQ074", "PBQ075", "PBQ076", "PBQ077", "PBQ078", 
              "PBQ079", "PBQ080", "PBQ081", "PBQ082", "PBQ083", "PBQ084", 
              "PBQ085", "PBQ086", "PBQ087", "PBQ088", "PBQ089", "PBQ090", 
              "PBQ091", "PBQ092", "PBQ093", "PBQ102", "PBQ110", "PBQ124")
all_subjects <- c(template_subject, subjects)
base_dir <- "Z:/Projects/058_Pelvis Bone Quality/data/"
demographics_df <- read.csv("C:/Users/telfe/Dropbox/My_Projects/Pelvic_Bone_Quality/demographics.csv")
demographics_df <- demographics_df %>% filter(subject_id %in% all_subjects)

## check files are present
for (i in seq_along(all_subjects)) {
  if (file.exists(paste0(base_dir, all_subjects[i],
                         "/pelvis_resampled.nii")) == FALSE) {
    stop(paste0(".nii file missing ", all_subjects[i]))
  }
  if (file.exists(paste0(base_dir, all_subjects[i],
                         "/surface_models/pelvisRight_surfaceModel1-5.ply")) == FALSE) {
    stop(paste0(".ply file missing ", all_subjects[i]))
  }
  if (file.exists(paste0(base_dir, all_subjects[i],
                         "/surface_models/", all_subjects[i], "_landmarks_right.fcsv")) == FALSE) {
    stop(paste0("landmark file missing ", all_subjects[i]))
  }
}

## surface mesh paths
surface_meshes <- c()
for (subject in seq_along(subjects)) {
  surface_meshes[subject] <- paste0(base_dir, subjects[subject], 
                                    "/surface_models/pelvisRight_surfaceModel1-5.ply")
}
template_surface_path <- paste0(base_dir, template_subject, 
                                "/surface_models/pelvisRight_surfaceModel1-5.ply")

## landmark paths
landmark_paths <- c()
for (subject in seq_along(subjects)) {
  landmark_paths[subject] <- paste0(base_dir, subjects[subject], "/surface_models/",
                                    subjects[subject], "_landmarks_right.fcsv")
}
template_landmarks_path <- paste0(base_dir, template_subject, "/surface_models/",
                                  template_subject, "_landmarks_right.fcsv")

## check landmarks are on bone model surfaces
landmark_check(template_surface_path, template_landmarks_path)
for (subject in seq_along(subjects)) {  
  print(as.character(subjects[subject]))
  landmark_check(surface_meshes[subject], landmark_paths[subject], threshold = 2)
}


# =============================================================================

# material coefficients
calibration_curves <- read.csv(paste0("C:/Users/telfe/Dropbox/My_Projects/", 
                                      "Pelvic_Bone_Quality/PelvicBoneQuality/ct_calibrations.csv"))
demographics_df <- demographics_df %>% rowwise() %>% 
  mutate(sigma = ct_coefficients(table_height, calibration_curves, scanner, "sigma")) %>%
  mutate(beta = ct_coefficients(table_height, calibration_curves, scanner, "beta"))
mat_coeff <- demographics_df %>% select(subject_id, sigma, beta)


# =============================================================================

# calculate density values for template mesh
## make template
template_surface <- import_ply(template_surface_path)
template_landmarks <- read.csv(template_landmarks_path, skip = 2)[, 2:4]
template <- surface_points_template(template_surface, template_landmarks, 10000)

## generate surface density values
### empty matrix to store all density values
density_mat <- matrix(NA, nrow = 10015, ncol = (length(subjects) + 1))

### calculate density
template_nifti <- readNIfTI(paste0(base_dir, template_subject, "/pelvis_resampled.nii"))
template_beta <- mat_coeff %>% filter(subject_id == template_subject) %>% pull(beta)
template_sigma <- mat_coeff %>% filter(subject_id == template_subject) %>% pull(sigma)
density_mat[, 1] <- surface_normal_intersect(template_surface, template, 3.0, template_nifti, 
                                             template_beta, template_sigma, 0.5)
rm(template_nifti)

# calculate density values for all meshes
## generate new surface points
new_pts <- list()
for (subject in seq_along(subjects)) {
  surface_mesh <- import_ply(surface_meshes[subject])
  lnmks <- read.csv(landmark_paths[subject], skip = 2)[2:4]
  new_pts[[subject]] <- surface_points_new(surface_mesh, lnmks, template)
  
}
## calculate density
for (subject in seq_along(subjects)) {
  surface_mesh <- import_ply(surface_meshes[subject])
  nifti <- readNIfTI(paste0(base_dir, subjects[subject], "/pelvis_resampled.nii"))
  beta_val <- mat_coeff %>% filter(subject_id == subjects[subject]) %>% pull(beta)
  sigma_val <- mat_coeff %>% filter(subject_id == subjects[subject]) %>% pull(sigma)
  density_mat[, subject + 1] <- surface_normal_intersect(surface_mesh, new_pts[[subject]], 
                                                         3.0, nifti, beta_val, sigma_val, 0.5)
  print(subject)
}



# make density_mat into data frame
density_df <- as.data.frame(density_mat)
colnames(density_df) <- all_subjects

# export
saveRDS(density_df, file = "data/density_df.RData")
saveRDS(template, file = "data/template.RData")


# =============================================================================