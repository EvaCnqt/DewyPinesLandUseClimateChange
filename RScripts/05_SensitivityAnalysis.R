############################################################################
#
# This script processes the results of the projections of the dewy-pine
# IBMs used to perform a sensitivity analysis. We projected each natural 
# and anthropogenic population while perturbing (i.e. using projected values
# for climate variables) one vital rate at a time and leaving all other 
# vital rates unperturbed (i.e. using current values for climate variables).
#
# Author: Eva Conquet
#
###########################################################################

###########################################################################
#
# 1. House keeping and loading libraries and data ----
#
###########################################################################

## 1.1. House keeping ----
# -------------------

rm(list = ls())


## 1.2. Loading libraries ----
# -----------------------

library(ggplot2)
library(patchwork)


## 1.3. Loading data ----
# ------------------

ibm_results = read.csv("Output/Projections/IBM_Results.csv")
ibm_results = ibm_results[which(ibm_results$climate == "Control" &
                                ibm_results$scenario == "RCP 8.5"), ]



###########################################################################
#
# 2. Processing results ----
#
###########################################################################

n_sim = 100
n_years = 30

climate_change_models = c("CanESM5", 
                          "EC_Earth3", "FGOALS_G3", "GFDL_ESM4",
                          "GISS_E2_1_G", "INM_CM4_8", 
                          "IPSL_CM6A_LR", "MIROC6", 
                          "MPI_ESM1_2_LR", "MRI_ESM2_0", 
                          "NorESM2_MM")

populations_anthropogenic = c("Retin", "Prisoneros", "Bujeo", 
                              "MonteraTorero",
                              "SCarbDist")
populations_natural = c("SierraCarboneraY5", "SierraRetinY5", "Vertedero")


## 2.1. Merging results ----
# ---------------------

# Results data frame

sensitivityanalysis_results = expand.grid(year = seq(1, n_years),
                          simulation = seq(1, n_sim), 
                          climate = c(rep("Climate change", 
                                          times = length(climate_change_models))),
                          vital_rate = c("survival", "growth", 
                                         "flowering", "nbFlowers",
                                         "seedlingSize"),
                          habitat = c(rep("Natural",
                                          times = length(populations_natural)), 
                                      rep("Anthropogenic", 
                                          times = length(populations_anthropogenic))))


# Add populations

sensitivityanalysis_results$population = c(rep(populations_natural, each = n_years * n_sim * (length(climate_change_models)) * 5),
                           rep(populations_anthropogenic, each = n_years * n_sim * (length(climate_change_models)) * 5))


# Add climate models

sensitivityanalysis_results$climate_model = NA
sensitivityanalysis_results$climate_model[which(sensitivityanalysis_results$climate == "Climate change")] = rep(rep(climate_change_models, each = n_years * n_sim), 5)

sensitivityanalysis_results$treatment = paste(sensitivityanalysis_results$population,
                              sensitivityanalysis_results$climate_model,
                              sensitivityanalysis_results$vital_rate,
                              # sensitivityanalysis_results$simulation,
                              sep = "_")


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))][grep("survival", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))][grep("survival", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda = NA

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Retin" &
                             sensitivityanalysis_results$vital_rate == "survival")] = 
  c(t(rbind(ibm_retin_canesm5_sensitivity_survival[[1]]$log_lambda,
            ibm_retin_canesm5_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_retin_canesm5_sensitivity_survival[[3]]$log_lambda,
            # ibm_retin_canesm5_sensitivity_survival[[4]]$log_lambda,
            # ibm_retin_canesm5_sensitivity_survival[[5]]$log_lambda
            )),
    t(rbind(ibm_retin_ec_earth3_sensitivity_survival[[1]]$log_lambda,
            ibm_retin_ec_earth3_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_retin_ec_earth3_sensitivity_survival[[3]]$log_lambda,
            # ibm_retin_ec_earth3_sensitivity_survival[[4]]$log_lambda,
            # ibm_retin_ec_earth3_sensitivity_survival[[5]]$log_lambda
            )),
    t(rbind(ibm_retin_fgoals_g3_sensitivity_survival[[1]]$log_lambda,
            ibm_retin_fgoals_g3_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_retin_fgoals_g3_sensitivity_survival[[3]]$log_lambda,
            # ibm_retin_fgoals_g3_sensitivity_survival[[4]]$log_lambda,
            # ibm_retin_fgoals_g3_sensitivity_survival[[5]]$log_lambda
            )),
    t(rbind(ibm_retin_gfdl_esm4_sensitivity_survival[[1]]$log_lambda,
            ibm_retin_gfdl_esm4_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_retin_gfdl_esm4_sensitivity_survival[[3]]$log_lambda,
            # ibm_retin_gfdl_esm4_sensitivity_survival[[4]]$log_lambda,
            # ibm_retin_gfdl_esm4_sensitivity_survival[[5]]$log_lambda
            )),
    t(rbind(ibm_retin_giss_e2_1_g_sensitivity_survival[[1]]$log_lambda,
            ibm_retin_giss_e2_1_g_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_retin_giss_e2_1_g_sensitivity_survival[[3]]$log_lambda,
            # ibm_retin_giss_e2_1_g_sensitivity_survival[[4]]$log_lambda,
            # ibm_retin_giss_e2_1_g_sensitivity_survival[[5]]$log_lambda
            )),
    t(rbind(ibm_retin_inm_cm4_8_sensitivity_survival[[1]]$log_lambda,
            ibm_retin_inm_cm4_8_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_retin_inm_cm4_8_sensitivity_survival[[3]]$log_lambda,
            # ibm_retin_inm_cm4_8_sensitivity_survival[[4]]$log_lambda,
            # ibm_retin_inm_cm4_8_sensitivity_survival[[5]]$log_lambda
            )),
    t(rbind(ibm_retin_ipsl_cm6a_lr_sensitivity_survival[[1]]$log_lambda,
            ibm_retin_ipsl_cm6a_lr_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_retin_ipsl_cm6a_lr_sensitivity_survival[[3]]$log_lambda,
            # ibm_retin_ipsl_cm6a_lr_sensitivity_survival[[4]]$log_lambda,
            # ibm_retin_ipsl_cm6a_lr_sensitivity_survival[[5]]$log_lambda
            )),
    t(rbind(ibm_retin_miroc6_sensitivity_survival[[1]]$log_lambda,
            ibm_retin_miroc6_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_retin_miroc6_sensitivity_survival[[3]]$log_lambda,
            # ibm_retin_miroc6_sensitivity_survival[[4]]$log_lambda,
            # ibm_retin_miroc6_sensitivity_survival[[5]]$log_lambda
            )),
    t(rbind(ibm_retin_mpi_esm1_2_lr_sensitivity_survival[[1]]$log_lambda,
            ibm_retin_mpi_esm1_2_lr_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_retin_mpi_esm1_2_lr_sensitivity_survival[[3]]$log_lambda,
            # ibm_retin_mpi_esm1_2_lr_sensitivity_survival[[4]]$log_lambda,
            # ibm_retin_mpi_esm1_2_lr_sensitivity_survival[[5]]$log_lambda
            )),
    t(rbind(ibm_retin_mri_esm2_0_sensitivity_survival[[1]]$log_lambda,
            ibm_retin_mri_esm2_0_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_retin_mri_esm2_0_sensitivity_survival[[3]]$log_lambda,
            # ibm_retin_mri_esm2_0_sensitivity_survival[[4]]$log_lambda,
            # ibm_retin_mri_esm2_0_sensitivity_survival[[5]]$log_lambda
            )),
    t(rbind(ibm_retin_noresm2_mm_sensitivity_survival[[1]]$log_lambda,
            ibm_retin_noresm2_mm_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_retin_noresm2_mm_sensitivity_survival[[3]]$log_lambda,
            # ibm_retin_noresm2_mm_sensitivity_survival[[4]]$log_lambda,
            # ibm_retin_noresm2_mm_sensitivity_survival[[5]]$log_lambda
            )))

rm(list = grep("survival", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))][grep("growth", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))][grep("growth", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Retin" &
                               sensitivityanalysis_results$vital_rate == "growth")] = 
  c(t(rbind(ibm_retin_canesm5_sensitivity_growth[[1]]$log_lambda,
            ibm_retin_canesm5_sensitivity_growth[[2]]$log_lambda
            # ,
            # ibm_retin_canesm5_sensitivity_growth[[3]]$log_lambda,
            # ibm_retin_canesm5_sensitivity_growth[[4]]$log_lambda,
            # ibm_retin_canesm5_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_ec_earth3_sensitivity_growth[[1]]$log_lambda,
          ibm_retin_ec_earth3_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_retin_ec_earth3_sensitivity_growth[[3]]$log_lambda,
          # ibm_retin_ec_earth3_sensitivity_growth[[4]]$log_lambda,
          # ibm_retin_ec_earth3_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_fgoals_g3_sensitivity_growth[[1]]$log_lambda,
          ibm_retin_fgoals_g3_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_retin_fgoals_g3_sensitivity_growth[[3]]$log_lambda,
          # ibm_retin_fgoals_g3_sensitivity_growth[[4]]$log_lambda,
          # ibm_retin_fgoals_g3_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_gfdl_esm4_sensitivity_growth[[1]]$log_lambda,
          ibm_retin_gfdl_esm4_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_retin_gfdl_esm4_sensitivity_growth[[3]]$log_lambda,
          # ibm_retin_gfdl_esm4_sensitivity_growth[[4]]$log_lambda,
          # ibm_retin_gfdl_esm4_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_giss_e2_1_g_sensitivity_growth[[1]]$log_lambda,
          ibm_retin_giss_e2_1_g_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_retin_giss_e2_1_g_sensitivity_growth[[3]]$log_lambda,
          # ibm_retin_giss_e2_1_g_sensitivity_growth[[4]]$log_lambda,
          # ibm_retin_giss_e2_1_g_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_inm_cm4_8_sensitivity_growth[[1]]$log_lambda,
          ibm_retin_inm_cm4_8_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_retin_inm_cm4_8_sensitivity_growth[[3]]$log_lambda,
          # ibm_retin_inm_cm4_8_sensitivity_growth[[4]]$log_lambda,
          # ibm_retin_inm_cm4_8_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_ipsl_cm6a_lr_sensitivity_growth[[1]]$log_lambda,
          ibm_retin_ipsl_cm6a_lr_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_retin_ipsl_cm6a_lr_sensitivity_growth[[3]]$log_lambda,
          # ibm_retin_ipsl_cm6a_lr_sensitivity_growth[[4]]$log_lambda,
          # ibm_retin_ipsl_cm6a_lr_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_miroc6_sensitivity_growth[[1]]$log_lambda,
          ibm_retin_miroc6_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_retin_miroc6_sensitivity_growth[[3]]$log_lambda,
          # ibm_retin_miroc6_sensitivity_growth[[4]]$log_lambda,
          # ibm_retin_miroc6_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_mpi_esm1_2_lr_sensitivity_growth[[1]]$log_lambda,
          ibm_retin_mpi_esm1_2_lr_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_retin_mpi_esm1_2_lr_sensitivity_growth[[3]]$log_lambda,
          # ibm_retin_mpi_esm1_2_lr_sensitivity_growth[[4]]$log_lambda,
          # ibm_retin_mpi_esm1_2_lr_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_mri_esm2_0_sensitivity_growth[[1]]$log_lambda,
          ibm_retin_mri_esm2_0_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_retin_mri_esm2_0_sensitivity_growth[[3]]$log_lambda,
          # ibm_retin_mri_esm2_0_sensitivity_growth[[4]]$log_lambda,
          # ibm_retin_mri_esm2_0_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_noresm2_mm_sensitivity_growth[[1]]$log_lambda,
          ibm_retin_noresm2_mm_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_retin_noresm2_mm_sensitivity_growth[[3]]$log_lambda,
          # ibm_retin_noresm2_mm_sensitivity_growth[[4]]$log_lambda,
          # ibm_retin_noresm2_mm_sensitivity_growth[[5]]$log_lambda
  )))

rm(list = grep("growth", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))][grep("flowering", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))][grep("flowering", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Retin" &
                               sensitivityanalysis_results$vital_rate == "flowering")] = 
  c(t(rbind(ibm_retin_canesm5_sensitivity_flowering[[1]]$log_lambda,
            ibm_retin_canesm5_sensitivity_flowering[[2]]$log_lambda
            # ,
            # ibm_retin_canesm5_sensitivity_flowering[[3]]$log_lambda,
            # ibm_retin_canesm5_sensitivity_flowering[[4]]$log_lambda,
            # ibm_retin_canesm5_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_ec_earth3_sensitivity_flowering[[1]]$log_lambda,
          ibm_retin_ec_earth3_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_retin_ec_earth3_sensitivity_flowering[[3]]$log_lambda,
          # ibm_retin_ec_earth3_sensitivity_flowering[[4]]$log_lambda,
          # ibm_retin_ec_earth3_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_fgoals_g3_sensitivity_flowering[[1]]$log_lambda,
          ibm_retin_fgoals_g3_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_retin_fgoals_g3_sensitivity_flowering[[3]]$log_lambda,
          # ibm_retin_fgoals_g3_sensitivity_flowering[[4]]$log_lambda,
          # ibm_retin_fgoals_g3_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_gfdl_esm4_sensitivity_flowering[[1]]$log_lambda,
          ibm_retin_gfdl_esm4_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_retin_gfdl_esm4_sensitivity_flowering[[3]]$log_lambda,
          # ibm_retin_gfdl_esm4_sensitivity_flowering[[4]]$log_lambda,
          # ibm_retin_gfdl_esm4_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_giss_e2_1_g_sensitivity_flowering[[1]]$log_lambda,
          ibm_retin_giss_e2_1_g_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_retin_giss_e2_1_g_sensitivity_flowering[[3]]$log_lambda,
          # ibm_retin_giss_e2_1_g_sensitivity_flowering[[4]]$log_lambda,
          # ibm_retin_giss_e2_1_g_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_inm_cm4_8_sensitivity_flowering[[1]]$log_lambda,
          ibm_retin_inm_cm4_8_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_retin_inm_cm4_8_sensitivity_flowering[[3]]$log_lambda,
          # ibm_retin_inm_cm4_8_sensitivity_flowering[[4]]$log_lambda,
          # ibm_retin_inm_cm4_8_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_ipsl_cm6a_lr_sensitivity_flowering[[1]]$log_lambda,
          ibm_retin_ipsl_cm6a_lr_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_retin_ipsl_cm6a_lr_sensitivity_flowering[[3]]$log_lambda,
          # ibm_retin_ipsl_cm6a_lr_sensitivity_flowering[[4]]$log_lambda,
          # ibm_retin_ipsl_cm6a_lr_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_miroc6_sensitivity_flowering[[1]]$log_lambda,
          ibm_retin_miroc6_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_retin_miroc6_sensitivity_flowering[[3]]$log_lambda,
          # ibm_retin_miroc6_sensitivity_flowering[[4]]$log_lambda,
          # ibm_retin_miroc6_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_mpi_esm1_2_lr_sensitivity_flowering[[1]]$log_lambda,
          ibm_retin_mpi_esm1_2_lr_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_retin_mpi_esm1_2_lr_sensitivity_flowering[[3]]$log_lambda,
          # ibm_retin_mpi_esm1_2_lr_sensitivity_flowering[[4]]$log_lambda,
          # ibm_retin_mpi_esm1_2_lr_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_mri_esm2_0_sensitivity_flowering[[1]]$log_lambda,
          ibm_retin_mri_esm2_0_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_retin_mri_esm2_0_sensitivity_flowering[[3]]$log_lambda,
          # ibm_retin_mri_esm2_0_sensitivity_flowering[[4]]$log_lambda,
          # ibm_retin_mri_esm2_0_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_noresm2_mm_sensitivity_flowering[[1]]$log_lambda,
          ibm_retin_noresm2_mm_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_retin_noresm2_mm_sensitivity_flowering[[3]]$log_lambda,
          # ibm_retin_noresm2_mm_sensitivity_flowering[[4]]$log_lambda,
          # ibm_retin_noresm2_mm_sensitivity_flowering[[5]]$log_lambda
  )))

rm(list = grep("flowering", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))][grep("nbFlowers", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))][grep("nbFlowers", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Retin" &
                               sensitivityanalysis_results$vital_rate == "nbFlowers")] = 
  c(t(rbind(ibm_retin_canesm5_sensitivity_nbFlowers[[1]]$log_lambda,
            ibm_retin_canesm5_sensitivity_nbFlowers[[2]]$log_lambda
            # ,
            # ibm_retin_canesm5_sensitivity_nbFlowers[[3]]$log_lambda,
            # ibm_retin_canesm5_sensitivity_nbFlowers[[4]]$log_lambda,
            # ibm_retin_canesm5_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_ec_earth3_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_retin_ec_earth3_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_retin_ec_earth3_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_retin_ec_earth3_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_retin_ec_earth3_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_fgoals_g3_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_retin_fgoals_g3_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_retin_fgoals_g3_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_retin_fgoals_g3_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_retin_fgoals_g3_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_gfdl_esm4_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_retin_gfdl_esm4_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_retin_gfdl_esm4_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_retin_gfdl_esm4_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_retin_gfdl_esm4_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_giss_e2_1_g_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_retin_giss_e2_1_g_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_retin_giss_e2_1_g_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_retin_giss_e2_1_g_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_retin_giss_e2_1_g_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_inm_cm4_8_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_retin_inm_cm4_8_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_retin_inm_cm4_8_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_retin_inm_cm4_8_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_retin_inm_cm4_8_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_ipsl_cm6a_lr_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_retin_ipsl_cm6a_lr_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_retin_ipsl_cm6a_lr_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_retin_ipsl_cm6a_lr_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_retin_ipsl_cm6a_lr_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_miroc6_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_retin_miroc6_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_retin_miroc6_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_retin_miroc6_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_retin_miroc6_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_mpi_esm1_2_lr_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_retin_mpi_esm1_2_lr_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_retin_mpi_esm1_2_lr_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_retin_mpi_esm1_2_lr_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_retin_mpi_esm1_2_lr_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_mri_esm2_0_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_retin_mri_esm2_0_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_retin_mri_esm2_0_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_retin_mri_esm2_0_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_retin_mri_esm2_0_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_noresm2_mm_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_retin_noresm2_mm_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_retin_noresm2_mm_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_retin_noresm2_mm_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_retin_noresm2_mm_sensitivity_nbFlowers[[5]]$log_lambda
  )))

rm(list = grep("nbFlowers", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))][grep("seedlingSize", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))][grep("seedlingSize", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Retin", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Retin" &
                               sensitivityanalysis_results$vital_rate == "seedlingSize")] = 
  c(t(rbind(ibm_retin_canesm5_sensitivity_seedlingSize[[1]]$log_lambda,
            ibm_retin_canesm5_sensitivity_seedlingSize[[2]]$log_lambda
            # ,
            # ibm_retin_canesm5_sensitivity_seedlingSize[[3]]$log_lambda,
            # ibm_retin_canesm5_sensitivity_seedlingSize[[4]]$log_lambda,
            # ibm_retin_canesm5_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_ec_earth3_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_retin_ec_earth3_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_retin_ec_earth3_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_retin_ec_earth3_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_retin_ec_earth3_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_fgoals_g3_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_retin_fgoals_g3_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_retin_fgoals_g3_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_retin_fgoals_g3_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_retin_fgoals_g3_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_gfdl_esm4_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_retin_gfdl_esm4_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_retin_gfdl_esm4_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_retin_gfdl_esm4_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_retin_gfdl_esm4_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_giss_e2_1_g_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_retin_giss_e2_1_g_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_retin_giss_e2_1_g_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_retin_giss_e2_1_g_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_retin_giss_e2_1_g_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_inm_cm4_8_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_retin_inm_cm4_8_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_retin_inm_cm4_8_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_retin_inm_cm4_8_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_retin_inm_cm4_8_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_ipsl_cm6a_lr_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_retin_ipsl_cm6a_lr_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_retin_ipsl_cm6a_lr_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_retin_ipsl_cm6a_lr_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_retin_ipsl_cm6a_lr_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_miroc6_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_retin_miroc6_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_retin_miroc6_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_retin_miroc6_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_retin_miroc6_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_mpi_esm1_2_lr_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_retin_mpi_esm1_2_lr_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_retin_mpi_esm1_2_lr_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_retin_mpi_esm1_2_lr_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_retin_mpi_esm1_2_lr_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_mri_esm2_0_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_retin_mri_esm2_0_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_retin_mri_esm2_0_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_retin_mri_esm2_0_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_retin_mri_esm2_0_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_retin_noresm2_mm_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_retin_noresm2_mm_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_retin_noresm2_mm_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_retin_noresm2_mm_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_retin_noresm2_mm_sensitivity_seedlingSize[[5]]$log_lambda
  )))

rm(list = grep("seedlingSize", ls() , value = TRUE, invert = FALSE))




for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))][grep("survival", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))][grep("survival", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Prisoneros" &
                               sensitivityanalysis_results$vital_rate == "survival")] = 
  c(t(rbind(ibm_prisoneros_canesm5_sensitivity_survival[[1]]$log_lambda,
            ibm_prisoneros_canesm5_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_prisoneros_canesm5_sensitivity_survival[[3]]$log_lambda,
            # ibm_prisoneros_canesm5_sensitivity_survival[[4]]$log_lambda,
            # ibm_prisoneros_canesm5_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_ec_earth3_sensitivity_survival[[1]]$log_lambda,
          ibm_prisoneros_ec_earth3_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_prisoneros_ec_earth3_sensitivity_survival[[3]]$log_lambda,
          # ibm_prisoneros_ec_earth3_sensitivity_survival[[4]]$log_lambda,
          # ibm_prisoneros_ec_earth3_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_fgoals_g3_sensitivity_survival[[1]]$log_lambda,
          ibm_prisoneros_fgoals_g3_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_prisoneros_fgoals_g3_sensitivity_survival[[3]]$log_lambda,
          # ibm_prisoneros_fgoals_g3_sensitivity_survival[[4]]$log_lambda,
          # ibm_prisoneros_fgoals_g3_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_gfdl_esm4_sensitivity_survival[[1]]$log_lambda,
          ibm_prisoneros_gfdl_esm4_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_prisoneros_gfdl_esm4_sensitivity_survival[[3]]$log_lambda,
          # ibm_prisoneros_gfdl_esm4_sensitivity_survival[[4]]$log_lambda,
          # ibm_prisoneros_gfdl_esm4_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_giss_e2_1_g_sensitivity_survival[[1]]$log_lambda,
          ibm_prisoneros_giss_e2_1_g_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_prisoneros_giss_e2_1_g_sensitivity_survival[[3]]$log_lambda,
          # ibm_prisoneros_giss_e2_1_g_sensitivity_survival[[4]]$log_lambda,
          # ibm_prisoneros_giss_e2_1_g_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_inm_cm4_8_sensitivity_survival[[1]]$log_lambda,
          ibm_prisoneros_inm_cm4_8_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_prisoneros_inm_cm4_8_sensitivity_survival[[3]]$log_lambda,
          # ibm_prisoneros_inm_cm4_8_sensitivity_survival[[4]]$log_lambda,
          # ibm_prisoneros_inm_cm4_8_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_ipsl_cm6a_lr_sensitivity_survival[[1]]$log_lambda,
          ibm_prisoneros_ipsl_cm6a_lr_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_prisoneros_ipsl_cm6a_lr_sensitivity_survival[[3]]$log_lambda,
          # ibm_prisoneros_ipsl_cm6a_lr_sensitivity_survival[[4]]$log_lambda,
          # ibm_prisoneros_ipsl_cm6a_lr_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_miroc6_sensitivity_survival[[1]]$log_lambda,
          ibm_prisoneros_miroc6_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_prisoneros_miroc6_sensitivity_survival[[3]]$log_lambda,
          # ibm_prisoneros_miroc6_sensitivity_survival[[4]]$log_lambda,
          # ibm_prisoneros_miroc6_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_mpi_esm1_2_lr_sensitivity_survival[[1]]$log_lambda,
          ibm_prisoneros_mpi_esm1_2_lr_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_prisoneros_mpi_esm1_2_lr_sensitivity_survival[[3]]$log_lambda,
          # ibm_prisoneros_mpi_esm1_2_lr_sensitivity_survival[[4]]$log_lambda,
          # ibm_prisoneros_mpi_esm1_2_lr_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_mri_esm2_0_sensitivity_survival[[1]]$log_lambda,
          ibm_prisoneros_mri_esm2_0_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_prisoneros_mri_esm2_0_sensitivity_survival[[3]]$log_lambda,
          # ibm_prisoneros_mri_esm2_0_sensitivity_survival[[4]]$log_lambda,
          # ibm_prisoneros_mri_esm2_0_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_noresm2_mm_sensitivity_survival[[1]]$log_lambda,
          ibm_prisoneros_noresm2_mm_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_prisoneros_noresm2_mm_sensitivity_survival[[3]]$log_lambda,
          # ibm_prisoneros_noresm2_mm_sensitivity_survival[[4]]$log_lambda,
          # ibm_prisoneros_noresm2_mm_sensitivity_survival[[5]]$log_lambda
  )))

rm(list = grep("survival", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))][grep("growth", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))][grep("growth", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Prisoneros" &
                               sensitivityanalysis_results$vital_rate == "growth")] = 
  c(t(rbind(ibm_prisoneros_canesm5_sensitivity_growth[[1]]$log_lambda,
            ibm_prisoneros_canesm5_sensitivity_growth[[2]]$log_lambda
            # ,
            # ibm_prisoneros_canesm5_sensitivity_growth[[3]]$log_lambda,
            # ibm_prisoneros_canesm5_sensitivity_growth[[4]]$log_lambda,
            # ibm_prisoneros_canesm5_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_ec_earth3_sensitivity_growth[[1]]$log_lambda,
          ibm_prisoneros_ec_earth3_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_prisoneros_ec_earth3_sensitivity_growth[[3]]$log_lambda,
          # ibm_prisoneros_ec_earth3_sensitivity_growth[[4]]$log_lambda,
          # ibm_prisoneros_ec_earth3_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_fgoals_g3_sensitivity_growth[[1]]$log_lambda,
          ibm_prisoneros_fgoals_g3_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_prisoneros_fgoals_g3_sensitivity_growth[[3]]$log_lambda,
          # ibm_prisoneros_fgoals_g3_sensitivity_growth[[4]]$log_lambda,
          # ibm_prisoneros_fgoals_g3_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_gfdl_esm4_sensitivity_growth[[1]]$log_lambda,
          ibm_prisoneros_gfdl_esm4_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_prisoneros_gfdl_esm4_sensitivity_growth[[3]]$log_lambda,
          # ibm_prisoneros_gfdl_esm4_sensitivity_growth[[4]]$log_lambda,
          # ibm_prisoneros_gfdl_esm4_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_giss_e2_1_g_sensitivity_growth[[1]]$log_lambda,
          ibm_prisoneros_giss_e2_1_g_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_prisoneros_giss_e2_1_g_sensitivity_growth[[3]]$log_lambda,
          # ibm_prisoneros_giss_e2_1_g_sensitivity_growth[[4]]$log_lambda,
          # ibm_prisoneros_giss_e2_1_g_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_inm_cm4_8_sensitivity_growth[[1]]$log_lambda,
          ibm_prisoneros_inm_cm4_8_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_prisoneros_inm_cm4_8_sensitivity_growth[[3]]$log_lambda,
          # ibm_prisoneros_inm_cm4_8_sensitivity_growth[[4]]$log_lambda,
          # ibm_prisoneros_inm_cm4_8_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_ipsl_cm6a_lr_sensitivity_growth[[1]]$log_lambda,
          ibm_prisoneros_ipsl_cm6a_lr_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_prisoneros_ipsl_cm6a_lr_sensitivity_growth[[3]]$log_lambda,
          # ibm_prisoneros_ipsl_cm6a_lr_sensitivity_growth[[4]]$log_lambda,
          # ibm_prisoneros_ipsl_cm6a_lr_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_miroc6_sensitivity_growth[[1]]$log_lambda,
          ibm_prisoneros_miroc6_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_prisoneros_miroc6_sensitivity_growth[[3]]$log_lambda,
          # ibm_prisoneros_miroc6_sensitivity_growth[[4]]$log_lambda,
          # ibm_prisoneros_miroc6_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_mpi_esm1_2_lr_sensitivity_growth[[1]]$log_lambda,
          ibm_prisoneros_mpi_esm1_2_lr_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_prisoneros_mpi_esm1_2_lr_sensitivity_growth[[3]]$log_lambda,
          # ibm_prisoneros_mpi_esm1_2_lr_sensitivity_growth[[4]]$log_lambda,
          # ibm_prisoneros_mpi_esm1_2_lr_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_mri_esm2_0_sensitivity_growth[[1]]$log_lambda,
          ibm_prisoneros_mri_esm2_0_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_prisoneros_mri_esm2_0_sensitivity_growth[[3]]$log_lambda,
          # ibm_prisoneros_mri_esm2_0_sensitivity_growth[[4]]$log_lambda,
          # ibm_prisoneros_mri_esm2_0_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_noresm2_mm_sensitivity_growth[[1]]$log_lambda,
          ibm_prisoneros_noresm2_mm_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_prisoneros_noresm2_mm_sensitivity_growth[[3]]$log_lambda,
          # ibm_prisoneros_noresm2_mm_sensitivity_growth[[4]]$log_lambda,
          # ibm_prisoneros_noresm2_mm_sensitivity_growth[[5]]$log_lambda
  )))

rm(list = grep("growth", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))][grep("flowering", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))][grep("flowering", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Prisoneros" &
                               sensitivityanalysis_results$vital_rate == "flowering")] = 
  c(t(rbind(ibm_prisoneros_canesm5_sensitivity_flowering[[1]]$log_lambda,
            ibm_prisoneros_canesm5_sensitivity_flowering[[2]]$log_lambda
            # ,
            # ibm_prisoneros_canesm5_sensitivity_flowering[[3]]$log_lambda,
            # ibm_prisoneros_canesm5_sensitivity_flowering[[4]]$log_lambda,
            # ibm_prisoneros_canesm5_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_ec_earth3_sensitivity_flowering[[1]]$log_lambda,
          ibm_prisoneros_ec_earth3_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_prisoneros_ec_earth3_sensitivity_flowering[[3]]$log_lambda,
          # ibm_prisoneros_ec_earth3_sensitivity_flowering[[4]]$log_lambda,
          # ibm_prisoneros_ec_earth3_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_fgoals_g3_sensitivity_flowering[[1]]$log_lambda,
          ibm_prisoneros_fgoals_g3_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_prisoneros_fgoals_g3_sensitivity_flowering[[3]]$log_lambda,
          # ibm_prisoneros_fgoals_g3_sensitivity_flowering[[4]]$log_lambda,
          # ibm_prisoneros_fgoals_g3_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_gfdl_esm4_sensitivity_flowering[[1]]$log_lambda,
          ibm_prisoneros_gfdl_esm4_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_prisoneros_gfdl_esm4_sensitivity_flowering[[3]]$log_lambda,
          # ibm_prisoneros_gfdl_esm4_sensitivity_flowering[[4]]$log_lambda,
          # ibm_prisoneros_gfdl_esm4_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_giss_e2_1_g_sensitivity_flowering[[1]]$log_lambda,
          ibm_prisoneros_giss_e2_1_g_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_prisoneros_giss_e2_1_g_sensitivity_flowering[[3]]$log_lambda,
          # ibm_prisoneros_giss_e2_1_g_sensitivity_flowering[[4]]$log_lambda,
          # ibm_prisoneros_giss_e2_1_g_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_inm_cm4_8_sensitivity_flowering[[1]]$log_lambda,
          ibm_prisoneros_inm_cm4_8_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_prisoneros_inm_cm4_8_sensitivity_flowering[[3]]$log_lambda,
          # ibm_prisoneros_inm_cm4_8_sensitivity_flowering[[4]]$log_lambda,
          # ibm_prisoneros_inm_cm4_8_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_ipsl_cm6a_lr_sensitivity_flowering[[1]]$log_lambda,
          ibm_prisoneros_ipsl_cm6a_lr_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_prisoneros_ipsl_cm6a_lr_sensitivity_flowering[[3]]$log_lambda,
          # ibm_prisoneros_ipsl_cm6a_lr_sensitivity_flowering[[4]]$log_lambda,
          # ibm_prisoneros_ipsl_cm6a_lr_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_miroc6_sensitivity_flowering[[1]]$log_lambda,
          ibm_prisoneros_miroc6_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_prisoneros_miroc6_sensitivity_flowering[[3]]$log_lambda,
          # ibm_prisoneros_miroc6_sensitivity_flowering[[4]]$log_lambda,
          # ibm_prisoneros_miroc6_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_mpi_esm1_2_lr_sensitivity_flowering[[1]]$log_lambda,
          ibm_prisoneros_mpi_esm1_2_lr_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_prisoneros_mpi_esm1_2_lr_sensitivity_flowering[[3]]$log_lambda,
          # ibm_prisoneros_mpi_esm1_2_lr_sensitivity_flowering[[4]]$log_lambda,
          # ibm_prisoneros_mpi_esm1_2_lr_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_mri_esm2_0_sensitivity_flowering[[1]]$log_lambda,
          ibm_prisoneros_mri_esm2_0_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_prisoneros_mri_esm2_0_sensitivity_flowering[[3]]$log_lambda,
          # ibm_prisoneros_mri_esm2_0_sensitivity_flowering[[4]]$log_lambda,
          # ibm_prisoneros_mri_esm2_0_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_noresm2_mm_sensitivity_flowering[[1]]$log_lambda,
          ibm_prisoneros_noresm2_mm_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_prisoneros_noresm2_mm_sensitivity_flowering[[3]]$log_lambda,
          # ibm_prisoneros_noresm2_mm_sensitivity_flowering[[4]]$log_lambda,
          # ibm_prisoneros_noresm2_mm_sensitivity_flowering[[5]]$log_lambda
  )))

rm(list = grep("flowering", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))][grep("nbFlowers", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))][grep("nbFlowers", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Prisoneros" &
                               sensitivityanalysis_results$vital_rate == "nbFlowers")] = 
  c(t(rbind(ibm_prisoneros_canesm5_sensitivity_nbFlowers[[1]]$log_lambda,
            ibm_prisoneros_canesm5_sensitivity_nbFlowers[[2]]$log_lambda
            # ,
            # ibm_prisoneros_canesm5_sensitivity_nbFlowers[[3]]$log_lambda,
            # ibm_prisoneros_canesm5_sensitivity_nbFlowers[[4]]$log_lambda,
            # ibm_prisoneros_canesm5_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_ec_earth3_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_prisoneros_ec_earth3_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_prisoneros_ec_earth3_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_prisoneros_ec_earth3_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_prisoneros_ec_earth3_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_fgoals_g3_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_prisoneros_fgoals_g3_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_prisoneros_fgoals_g3_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_prisoneros_fgoals_g3_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_prisoneros_fgoals_g3_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_gfdl_esm4_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_prisoneros_gfdl_esm4_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_prisoneros_gfdl_esm4_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_prisoneros_gfdl_esm4_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_prisoneros_gfdl_esm4_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_giss_e2_1_g_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_prisoneros_giss_e2_1_g_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_prisoneros_giss_e2_1_g_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_prisoneros_giss_e2_1_g_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_prisoneros_giss_e2_1_g_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_inm_cm4_8_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_prisoneros_inm_cm4_8_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_prisoneros_inm_cm4_8_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_prisoneros_inm_cm4_8_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_prisoneros_inm_cm4_8_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_ipsl_cm6a_lr_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_prisoneros_ipsl_cm6a_lr_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_prisoneros_ipsl_cm6a_lr_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_prisoneros_ipsl_cm6a_lr_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_prisoneros_ipsl_cm6a_lr_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_miroc6_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_prisoneros_miroc6_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_prisoneros_miroc6_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_prisoneros_miroc6_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_prisoneros_miroc6_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_mpi_esm1_2_lr_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_prisoneros_mpi_esm1_2_lr_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_prisoneros_mpi_esm1_2_lr_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_prisoneros_mpi_esm1_2_lr_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_prisoneros_mpi_esm1_2_lr_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_mri_esm2_0_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_prisoneros_mri_esm2_0_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_prisoneros_mri_esm2_0_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_prisoneros_mri_esm2_0_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_prisoneros_mri_esm2_0_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_noresm2_mm_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_prisoneros_noresm2_mm_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_prisoneros_noresm2_mm_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_prisoneros_noresm2_mm_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_prisoneros_noresm2_mm_sensitivity_nbFlowers[[5]]$log_lambda
  )))

rm(list = grep("nbFlowers", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))][grep("seedlingSize", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))][grep("seedlingSize", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Prisoneros", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Prisoneros" &
                               sensitivityanalysis_results$vital_rate == "seedlingSize")] = 
  c(t(rbind(ibm_prisoneros_canesm5_sensitivity_seedlingSize[[1]]$log_lambda,
            ibm_prisoneros_canesm5_sensitivity_seedlingSize[[2]]$log_lambda
            # ,
            # ibm_prisoneros_canesm5_sensitivity_seedlingSize[[3]]$log_lambda,
            # ibm_prisoneros_canesm5_sensitivity_seedlingSize[[4]]$log_lambda,
            # ibm_prisoneros_canesm5_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_ec_earth3_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_prisoneros_ec_earth3_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_prisoneros_ec_earth3_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_prisoneros_ec_earth3_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_prisoneros_ec_earth3_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_fgoals_g3_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_prisoneros_fgoals_g3_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_prisoneros_fgoals_g3_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_prisoneros_fgoals_g3_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_prisoneros_fgoals_g3_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_gfdl_esm4_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_prisoneros_gfdl_esm4_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_prisoneros_gfdl_esm4_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_prisoneros_gfdl_esm4_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_prisoneros_gfdl_esm4_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_giss_e2_1_g_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_prisoneros_giss_e2_1_g_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_prisoneros_giss_e2_1_g_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_prisoneros_giss_e2_1_g_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_prisoneros_giss_e2_1_g_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_inm_cm4_8_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_prisoneros_inm_cm4_8_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_prisoneros_inm_cm4_8_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_prisoneros_inm_cm4_8_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_prisoneros_inm_cm4_8_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_ipsl_cm6a_lr_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_prisoneros_ipsl_cm6a_lr_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_prisoneros_ipsl_cm6a_lr_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_prisoneros_ipsl_cm6a_lr_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_prisoneros_ipsl_cm6a_lr_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_miroc6_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_prisoneros_miroc6_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_prisoneros_miroc6_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_prisoneros_miroc6_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_prisoneros_miroc6_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_mpi_esm1_2_lr_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_prisoneros_mpi_esm1_2_lr_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_prisoneros_mpi_esm1_2_lr_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_prisoneros_mpi_esm1_2_lr_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_prisoneros_mpi_esm1_2_lr_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_mri_esm2_0_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_prisoneros_mri_esm2_0_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_prisoneros_mri_esm2_0_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_prisoneros_mri_esm2_0_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_prisoneros_mri_esm2_0_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_prisoneros_noresm2_mm_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_prisoneros_noresm2_mm_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_prisoneros_noresm2_mm_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_prisoneros_noresm2_mm_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_prisoneros_noresm2_mm_sensitivity_seedlingSize[[5]]$log_lambda
  )))

rm(list = grep("seedlingSize", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))][grep("survival", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))][grep("survival", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Bujeo" &
                               sensitivityanalysis_results$vital_rate == "survival")] = 
  c(t(rbind(ibm_bujeo_canesm5_sensitivity_survival[[1]]$log_lambda,
            ibm_bujeo_canesm5_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_bujeo_canesm5_sensitivity_survival[[3]]$log_lambda,
            # ibm_bujeo_canesm5_sensitivity_survival[[4]]$log_lambda,
            # ibm_bujeo_canesm5_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_ec_earth3_sensitivity_survival[[1]]$log_lambda,
          ibm_bujeo_ec_earth3_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_bujeo_ec_earth3_sensitivity_survival[[3]]$log_lambda,
          # ibm_bujeo_ec_earth3_sensitivity_survival[[4]]$log_lambda,
          # ibm_bujeo_ec_earth3_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_fgoals_g3_sensitivity_survival[[1]]$log_lambda,
          ibm_bujeo_fgoals_g3_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_bujeo_fgoals_g3_sensitivity_survival[[3]]$log_lambda,
          # ibm_bujeo_fgoals_g3_sensitivity_survival[[4]]$log_lambda,
          # ibm_bujeo_fgoals_g3_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_gfdl_esm4_sensitivity_survival[[1]]$log_lambda,
          ibm_bujeo_gfdl_esm4_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_bujeo_gfdl_esm4_sensitivity_survival[[3]]$log_lambda,
          # ibm_bujeo_gfdl_esm4_sensitivity_survival[[4]]$log_lambda,
          # ibm_bujeo_gfdl_esm4_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_giss_e2_1_g_sensitivity_survival[[1]]$log_lambda,
          ibm_bujeo_giss_e2_1_g_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_bujeo_giss_e2_1_g_sensitivity_survival[[3]]$log_lambda,
          # ibm_bujeo_giss_e2_1_g_sensitivity_survival[[4]]$log_lambda,
          # ibm_bujeo_giss_e2_1_g_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_inm_cm4_8_sensitivity_survival[[1]]$log_lambda,
          ibm_bujeo_inm_cm4_8_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_bujeo_inm_cm4_8_sensitivity_survival[[3]]$log_lambda,
          # ibm_bujeo_inm_cm4_8_sensitivity_survival[[4]]$log_lambda,
          # ibm_bujeo_inm_cm4_8_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_ipsl_cm6a_lr_sensitivity_survival[[1]]$log_lambda,
          ibm_bujeo_ipsl_cm6a_lr_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_bujeo_ipsl_cm6a_lr_sensitivity_survival[[3]]$log_lambda,
          # ibm_bujeo_ipsl_cm6a_lr_sensitivity_survival[[4]]$log_lambda,
          # ibm_bujeo_ipsl_cm6a_lr_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_miroc6_sensitivity_survival[[1]]$log_lambda,
          ibm_bujeo_miroc6_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_bujeo_miroc6_sensitivity_survival[[3]]$log_lambda,
          # ibm_bujeo_miroc6_sensitivity_survival[[4]]$log_lambda,
          # ibm_bujeo_miroc6_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_mpi_esm1_2_lr_sensitivity_survival[[1]]$log_lambda,
          ibm_bujeo_mpi_esm1_2_lr_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_bujeo_mpi_esm1_2_lr_sensitivity_survival[[3]]$log_lambda,
          # ibm_bujeo_mpi_esm1_2_lr_sensitivity_survival[[4]]$log_lambda,
          # ibm_bujeo_mpi_esm1_2_lr_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_mri_esm2_0_sensitivity_survival[[1]]$log_lambda,
          ibm_bujeo_mri_esm2_0_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_bujeo_mri_esm2_0_sensitivity_survival[[3]]$log_lambda,
          # ibm_bujeo_mri_esm2_0_sensitivity_survival[[4]]$log_lambda,
          # ibm_bujeo_mri_esm2_0_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_noresm2_mm_sensitivity_survival[[1]]$log_lambda,
          ibm_bujeo_noresm2_mm_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_bujeo_noresm2_mm_sensitivity_survival[[3]]$log_lambda,
          # ibm_bujeo_noresm2_mm_sensitivity_survival[[4]]$log_lambda,
          # ibm_bujeo_noresm2_mm_sensitivity_survival[[5]]$log_lambda
  )))

rm(list = grep("survival", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))][grep("growth", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))][grep("growth", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Bujeo" &
                               sensitivityanalysis_results$vital_rate == "growth")] = 
  c(t(rbind(ibm_bujeo_canesm5_sensitivity_growth[[1]]$log_lambda,
            ibm_bujeo_canesm5_sensitivity_growth[[2]]$log_lambda
            # ,
            # ibm_bujeo_canesm5_sensitivity_growth[[3]]$log_lambda,
            # ibm_bujeo_canesm5_sensitivity_growth[[4]]$log_lambda,
            # ibm_bujeo_canesm5_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_ec_earth3_sensitivity_growth[[1]]$log_lambda,
          ibm_bujeo_ec_earth3_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_bujeo_ec_earth3_sensitivity_growth[[3]]$log_lambda,
          # ibm_bujeo_ec_earth3_sensitivity_growth[[4]]$log_lambda,
          # ibm_bujeo_ec_earth3_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_fgoals_g3_sensitivity_growth[[1]]$log_lambda,
          ibm_bujeo_fgoals_g3_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_bujeo_fgoals_g3_sensitivity_growth[[3]]$log_lambda,
          # ibm_bujeo_fgoals_g3_sensitivity_growth[[4]]$log_lambda,
          # ibm_bujeo_fgoals_g3_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_gfdl_esm4_sensitivity_growth[[1]]$log_lambda,
          ibm_bujeo_gfdl_esm4_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_bujeo_gfdl_esm4_sensitivity_growth[[3]]$log_lambda,
          # ibm_bujeo_gfdl_esm4_sensitivity_growth[[4]]$log_lambda,
          # ibm_bujeo_gfdl_esm4_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_giss_e2_1_g_sensitivity_growth[[1]]$log_lambda,
          ibm_bujeo_giss_e2_1_g_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_bujeo_giss_e2_1_g_sensitivity_growth[[3]]$log_lambda,
          # ibm_bujeo_giss_e2_1_g_sensitivity_growth[[4]]$log_lambda,
          # ibm_bujeo_giss_e2_1_g_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_inm_cm4_8_sensitivity_growth[[1]]$log_lambda,
          ibm_bujeo_inm_cm4_8_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_bujeo_inm_cm4_8_sensitivity_growth[[3]]$log_lambda,
          # ibm_bujeo_inm_cm4_8_sensitivity_growth[[4]]$log_lambda,
          # ibm_bujeo_inm_cm4_8_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_ipsl_cm6a_lr_sensitivity_growth[[1]]$log_lambda,
          ibm_bujeo_ipsl_cm6a_lr_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_bujeo_ipsl_cm6a_lr_sensitivity_growth[[3]]$log_lambda,
          # ibm_bujeo_ipsl_cm6a_lr_sensitivity_growth[[4]]$log_lambda,
          # ibm_bujeo_ipsl_cm6a_lr_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_miroc6_sensitivity_growth[[1]]$log_lambda,
          ibm_bujeo_miroc6_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_bujeo_miroc6_sensitivity_growth[[3]]$log_lambda,
          # ibm_bujeo_miroc6_sensitivity_growth[[4]]$log_lambda,
          # ibm_bujeo_miroc6_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_mpi_esm1_2_lr_sensitivity_growth[[1]]$log_lambda,
          ibm_bujeo_mpi_esm1_2_lr_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_bujeo_mpi_esm1_2_lr_sensitivity_growth[[3]]$log_lambda,
          # ibm_bujeo_mpi_esm1_2_lr_sensitivity_growth[[4]]$log_lambda,
          # ibm_bujeo_mpi_esm1_2_lr_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_mri_esm2_0_sensitivity_growth[[1]]$log_lambda,
          ibm_bujeo_mri_esm2_0_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_bujeo_mri_esm2_0_sensitivity_growth[[3]]$log_lambda,
          # ibm_bujeo_mri_esm2_0_sensitivity_growth[[4]]$log_lambda,
          # ibm_bujeo_mri_esm2_0_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_noresm2_mm_sensitivity_growth[[1]]$log_lambda,
          ibm_bujeo_noresm2_mm_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_bujeo_noresm2_mm_sensitivity_growth[[3]]$log_lambda,
          # ibm_bujeo_noresm2_mm_sensitivity_growth[[4]]$log_lambda,
          # ibm_bujeo_noresm2_mm_sensitivity_growth[[5]]$log_lambda
  )))

rm(list = grep("growth", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))][grep("flowering", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))][grep("flowering", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Bujeo" &
                               sensitivityanalysis_results$vital_rate == "flowering")] = 
  c(t(rbind(ibm_bujeo_canesm5_sensitivity_flowering[[1]]$log_lambda,
            ibm_bujeo_canesm5_sensitivity_flowering[[2]]$log_lambda
            # ,
            # ibm_bujeo_canesm5_sensitivity_flowering[[3]]$log_lambda,
            # ibm_bujeo_canesm5_sensitivity_flowering[[4]]$log_lambda,
            # ibm_bujeo_canesm5_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_ec_earth3_sensitivity_flowering[[1]]$log_lambda,
          ibm_bujeo_ec_earth3_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_bujeo_ec_earth3_sensitivity_flowering[[3]]$log_lambda,
          # ibm_bujeo_ec_earth3_sensitivity_flowering[[4]]$log_lambda,
          # ibm_bujeo_ec_earth3_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_fgoals_g3_sensitivity_flowering[[1]]$log_lambda,
          ibm_bujeo_fgoals_g3_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_bujeo_fgoals_g3_sensitivity_flowering[[3]]$log_lambda,
          # ibm_bujeo_fgoals_g3_sensitivity_flowering[[4]]$log_lambda,
          # ibm_bujeo_fgoals_g3_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_gfdl_esm4_sensitivity_flowering[[1]]$log_lambda,
          ibm_bujeo_gfdl_esm4_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_bujeo_gfdl_esm4_sensitivity_flowering[[3]]$log_lambda,
          # ibm_bujeo_gfdl_esm4_sensitivity_flowering[[4]]$log_lambda,
          # ibm_bujeo_gfdl_esm4_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_giss_e2_1_g_sensitivity_flowering[[1]]$log_lambda,
          ibm_bujeo_giss_e2_1_g_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_bujeo_giss_e2_1_g_sensitivity_flowering[[3]]$log_lambda,
          # ibm_bujeo_giss_e2_1_g_sensitivity_flowering[[4]]$log_lambda,
          # ibm_bujeo_giss_e2_1_g_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_inm_cm4_8_sensitivity_flowering[[1]]$log_lambda,
          ibm_bujeo_inm_cm4_8_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_bujeo_inm_cm4_8_sensitivity_flowering[[3]]$log_lambda,
          # ibm_bujeo_inm_cm4_8_sensitivity_flowering[[4]]$log_lambda,
          # ibm_bujeo_inm_cm4_8_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_ipsl_cm6a_lr_sensitivity_flowering[[1]]$log_lambda,
          ibm_bujeo_ipsl_cm6a_lr_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_bujeo_ipsl_cm6a_lr_sensitivity_flowering[[3]]$log_lambda,
          # ibm_bujeo_ipsl_cm6a_lr_sensitivity_flowering[[4]]$log_lambda,
          # ibm_bujeo_ipsl_cm6a_lr_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_miroc6_sensitivity_flowering[[1]]$log_lambda,
          ibm_bujeo_miroc6_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_bujeo_miroc6_sensitivity_flowering[[3]]$log_lambda,
          # ibm_bujeo_miroc6_sensitivity_flowering[[4]]$log_lambda,
          # ibm_bujeo_miroc6_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_mpi_esm1_2_lr_sensitivity_flowering[[1]]$log_lambda,
          ibm_bujeo_mpi_esm1_2_lr_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_bujeo_mpi_esm1_2_lr_sensitivity_flowering[[3]]$log_lambda,
          # ibm_bujeo_mpi_esm1_2_lr_sensitivity_flowering[[4]]$log_lambda,
          # ibm_bujeo_mpi_esm1_2_lr_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_mri_esm2_0_sensitivity_flowering[[1]]$log_lambda,
          ibm_bujeo_mri_esm2_0_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_bujeo_mri_esm2_0_sensitivity_flowering[[3]]$log_lambda,
          # ibm_bujeo_mri_esm2_0_sensitivity_flowering[[4]]$log_lambda,
          # ibm_bujeo_mri_esm2_0_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_noresm2_mm_sensitivity_flowering[[1]]$log_lambda,
          ibm_bujeo_noresm2_mm_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_bujeo_noresm2_mm_sensitivity_flowering[[3]]$log_lambda,
          # ibm_bujeo_noresm2_mm_sensitivity_flowering[[4]]$log_lambda,
          # ibm_bujeo_noresm2_mm_sensitivity_flowering[[5]]$log_lambda
  )))

rm(list = grep("flowering", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))][grep("nbFlowers", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))][grep("nbFlowers", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Bujeo" &
                               sensitivityanalysis_results$vital_rate == "nbFlowers")] = 
  c(t(rbind(ibm_bujeo_canesm5_sensitivity_nbFlowers[[1]]$log_lambda,
            ibm_bujeo_canesm5_sensitivity_nbFlowers[[2]]$log_lambda
            # ,
            # ibm_bujeo_canesm5_sensitivity_nbFlowers[[3]]$log_lambda,
            # ibm_bujeo_canesm5_sensitivity_nbFlowers[[4]]$log_lambda,
            # ibm_bujeo_canesm5_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_ec_earth3_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_bujeo_ec_earth3_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_bujeo_ec_earth3_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_bujeo_ec_earth3_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_bujeo_ec_earth3_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_fgoals_g3_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_bujeo_fgoals_g3_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_bujeo_fgoals_g3_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_bujeo_fgoals_g3_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_bujeo_fgoals_g3_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_gfdl_esm4_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_bujeo_gfdl_esm4_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_bujeo_gfdl_esm4_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_bujeo_gfdl_esm4_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_bujeo_gfdl_esm4_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_giss_e2_1_g_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_bujeo_giss_e2_1_g_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_bujeo_giss_e2_1_g_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_bujeo_giss_e2_1_g_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_bujeo_giss_e2_1_g_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_inm_cm4_8_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_bujeo_inm_cm4_8_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_bujeo_inm_cm4_8_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_bujeo_inm_cm4_8_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_bujeo_inm_cm4_8_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_ipsl_cm6a_lr_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_bujeo_ipsl_cm6a_lr_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_bujeo_ipsl_cm6a_lr_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_bujeo_ipsl_cm6a_lr_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_bujeo_ipsl_cm6a_lr_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_miroc6_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_bujeo_miroc6_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_bujeo_miroc6_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_bujeo_miroc6_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_bujeo_miroc6_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_mpi_esm1_2_lr_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_bujeo_mpi_esm1_2_lr_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_bujeo_mpi_esm1_2_lr_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_bujeo_mpi_esm1_2_lr_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_bujeo_mpi_esm1_2_lr_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_mri_esm2_0_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_bujeo_mri_esm2_0_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_bujeo_mri_esm2_0_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_bujeo_mri_esm2_0_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_bujeo_mri_esm2_0_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_noresm2_mm_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_bujeo_noresm2_mm_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_bujeo_noresm2_mm_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_bujeo_noresm2_mm_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_bujeo_noresm2_mm_sensitivity_nbFlowers[[5]]$log_lambda
  )))

rm(list = grep("nbFlowers", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))][grep("seedlingSize", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))][grep("seedlingSize", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_Bujeo", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Bujeo" &
                               sensitivityanalysis_results$vital_rate == "seedlingSize")] = 
  c(t(rbind(ibm_bujeo_canesm5_sensitivity_seedlingSize[[1]]$log_lambda,
            ibm_bujeo_canesm5_sensitivity_seedlingSize[[2]]$log_lambda
            # ,
            # ibm_bujeo_canesm5_sensitivity_seedlingSize[[3]]$log_lambda,
            # ibm_bujeo_canesm5_sensitivity_seedlingSize[[4]]$log_lambda,
            # ibm_bujeo_canesm5_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_ec_earth3_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_bujeo_ec_earth3_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_bujeo_ec_earth3_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_bujeo_ec_earth3_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_bujeo_ec_earth3_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_fgoals_g3_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_bujeo_fgoals_g3_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_bujeo_fgoals_g3_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_bujeo_fgoals_g3_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_bujeo_fgoals_g3_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_gfdl_esm4_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_bujeo_gfdl_esm4_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_bujeo_gfdl_esm4_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_bujeo_gfdl_esm4_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_bujeo_gfdl_esm4_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_giss_e2_1_g_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_bujeo_giss_e2_1_g_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_bujeo_giss_e2_1_g_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_bujeo_giss_e2_1_g_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_bujeo_giss_e2_1_g_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_inm_cm4_8_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_bujeo_inm_cm4_8_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_bujeo_inm_cm4_8_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_bujeo_inm_cm4_8_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_bujeo_inm_cm4_8_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_ipsl_cm6a_lr_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_bujeo_ipsl_cm6a_lr_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_bujeo_ipsl_cm6a_lr_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_bujeo_ipsl_cm6a_lr_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_bujeo_ipsl_cm6a_lr_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_miroc6_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_bujeo_miroc6_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_bujeo_miroc6_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_bujeo_miroc6_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_bujeo_miroc6_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_mpi_esm1_2_lr_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_bujeo_mpi_esm1_2_lr_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_bujeo_mpi_esm1_2_lr_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_bujeo_mpi_esm1_2_lr_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_bujeo_mpi_esm1_2_lr_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_mri_esm2_0_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_bujeo_mri_esm2_0_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_bujeo_mri_esm2_0_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_bujeo_mri_esm2_0_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_bujeo_mri_esm2_0_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_bujeo_noresm2_mm_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_bujeo_noresm2_mm_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_bujeo_noresm2_mm_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_bujeo_noresm2_mm_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_bujeo_noresm2_mm_sensitivity_seedlingSize[[5]]$log_lambda
  )))

rm(list = grep("seedlingSize", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))][grep("survival", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))][grep("survival", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "MonteraTorero" &
                               sensitivityanalysis_results$vital_rate == "survival")] = 
  c(t(rbind(ibm_monteratorero_canesm5_sensitivity_survival[[1]]$log_lambda,
            ibm_monteratorero_canesm5_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_monteratorero_canesm5_sensitivity_survival[[3]]$log_lambda,
            # ibm_monteratorero_canesm5_sensitivity_survival[[4]]$log_lambda,
            # ibm_monteratorero_canesm5_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_ec_earth3_sensitivity_survival[[1]]$log_lambda,
          ibm_monteratorero_ec_earth3_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_monteratorero_ec_earth3_sensitivity_survival[[3]]$log_lambda,
          # ibm_monteratorero_ec_earth3_sensitivity_survival[[4]]$log_lambda,
          # ibm_monteratorero_ec_earth3_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_fgoals_g3_sensitivity_survival[[1]]$log_lambda,
          ibm_monteratorero_fgoals_g3_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_monteratorero_fgoals_g3_sensitivity_survival[[3]]$log_lambda,
          # ibm_monteratorero_fgoals_g3_sensitivity_survival[[4]]$log_lambda,
          # ibm_monteratorero_fgoals_g3_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_gfdl_esm4_sensitivity_survival[[1]]$log_lambda,
          ibm_monteratorero_gfdl_esm4_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_monteratorero_gfdl_esm4_sensitivity_survival[[3]]$log_lambda,
          # ibm_monteratorero_gfdl_esm4_sensitivity_survival[[4]]$log_lambda,
          # ibm_monteratorero_gfdl_esm4_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_giss_e2_1_g_sensitivity_survival[[1]]$log_lambda,
          ibm_monteratorero_giss_e2_1_g_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_monteratorero_giss_e2_1_g_sensitivity_survival[[3]]$log_lambda,
          # ibm_monteratorero_giss_e2_1_g_sensitivity_survival[[4]]$log_lambda,
          # ibm_monteratorero_giss_e2_1_g_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_inm_cm4_8_sensitivity_survival[[1]]$log_lambda,
          ibm_monteratorero_inm_cm4_8_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_monteratorero_inm_cm4_8_sensitivity_survival[[3]]$log_lambda,
          # ibm_monteratorero_inm_cm4_8_sensitivity_survival[[4]]$log_lambda,
          # ibm_monteratorero_inm_cm4_8_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_ipsl_cm6a_lr_sensitivity_survival[[1]]$log_lambda,
          ibm_monteratorero_ipsl_cm6a_lr_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_monteratorero_ipsl_cm6a_lr_sensitivity_survival[[3]]$log_lambda,
          # ibm_monteratorero_ipsl_cm6a_lr_sensitivity_survival[[4]]$log_lambda,
          # ibm_monteratorero_ipsl_cm6a_lr_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_miroc6_sensitivity_survival[[1]]$log_lambda,
          ibm_monteratorero_miroc6_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_monteratorero_miroc6_sensitivity_survival[[3]]$log_lambda,
          # ibm_monteratorero_miroc6_sensitivity_survival[[4]]$log_lambda,
          # ibm_monteratorero_miroc6_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_mpi_esm1_2_lr_sensitivity_survival[[1]]$log_lambda,
          ibm_monteratorero_mpi_esm1_2_lr_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_monteratorero_mpi_esm1_2_lr_sensitivity_survival[[3]]$log_lambda,
          # ibm_monteratorero_mpi_esm1_2_lr_sensitivity_survival[[4]]$log_lambda,
          # ibm_monteratorero_mpi_esm1_2_lr_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_mri_esm2_0_sensitivity_survival[[1]]$log_lambda,
          ibm_monteratorero_mri_esm2_0_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_monteratorero_mri_esm2_0_sensitivity_survival[[3]]$log_lambda,
          # ibm_monteratorero_mri_esm2_0_sensitivity_survival[[4]]$log_lambda,
          # ibm_monteratorero_mri_esm2_0_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_noresm2_mm_sensitivity_survival[[1]]$log_lambda,
          ibm_monteratorero_noresm2_mm_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_monteratorero_noresm2_mm_sensitivity_survival[[3]]$log_lambda,
          # ibm_monteratorero_noresm2_mm_sensitivity_survival[[4]]$log_lambda,
          # ibm_monteratorero_noresm2_mm_sensitivity_survival[[5]]$log_lambda
  )))

rm(list = grep("survival", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))][grep("growth", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))][grep("growth", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "MonteraTorero" &
                               sensitivityanalysis_results$vital_rate == "growth")] = 
  c(t(rbind(ibm_monteratorero_canesm5_sensitivity_growth[[1]]$log_lambda,
            ibm_monteratorero_canesm5_sensitivity_growth[[2]]$log_lambda
            # ,
            # ibm_monteratorero_canesm5_sensitivity_growth[[3]]$log_lambda,
            # ibm_monteratorero_canesm5_sensitivity_growth[[4]]$log_lambda,
            # ibm_monteratorero_canesm5_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_ec_earth3_sensitivity_growth[[1]]$log_lambda,
          ibm_monteratorero_ec_earth3_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_monteratorero_ec_earth3_sensitivity_growth[[3]]$log_lambda,
          # ibm_monteratorero_ec_earth3_sensitivity_growth[[4]]$log_lambda,
          # ibm_monteratorero_ec_earth3_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_fgoals_g3_sensitivity_growth[[1]]$log_lambda,
          ibm_monteratorero_fgoals_g3_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_monteratorero_fgoals_g3_sensitivity_growth[[3]]$log_lambda,
          # ibm_monteratorero_fgoals_g3_sensitivity_growth[[4]]$log_lambda,
          # ibm_monteratorero_fgoals_g3_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_gfdl_esm4_sensitivity_growth[[1]]$log_lambda,
          ibm_monteratorero_gfdl_esm4_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_monteratorero_gfdl_esm4_sensitivity_growth[[3]]$log_lambda,
          # ibm_monteratorero_gfdl_esm4_sensitivity_growth[[4]]$log_lambda,
          # ibm_monteratorero_gfdl_esm4_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_giss_e2_1_g_sensitivity_growth[[1]]$log_lambda,
          ibm_monteratorero_giss_e2_1_g_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_monteratorero_giss_e2_1_g_sensitivity_growth[[3]]$log_lambda,
          # ibm_monteratorero_giss_e2_1_g_sensitivity_growth[[4]]$log_lambda,
          # ibm_monteratorero_giss_e2_1_g_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_inm_cm4_8_sensitivity_growth[[1]]$log_lambda,
          ibm_monteratorero_inm_cm4_8_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_monteratorero_inm_cm4_8_sensitivity_growth[[3]]$log_lambda,
          # ibm_monteratorero_inm_cm4_8_sensitivity_growth[[4]]$log_lambda,
          # ibm_monteratorero_inm_cm4_8_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_ipsl_cm6a_lr_sensitivity_growth[[1]]$log_lambda,
          ibm_monteratorero_ipsl_cm6a_lr_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_monteratorero_ipsl_cm6a_lr_sensitivity_growth[[3]]$log_lambda,
          # ibm_monteratorero_ipsl_cm6a_lr_sensitivity_growth[[4]]$log_lambda,
          # ibm_monteratorero_ipsl_cm6a_lr_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_miroc6_sensitivity_growth[[1]]$log_lambda,
          ibm_monteratorero_miroc6_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_monteratorero_miroc6_sensitivity_growth[[3]]$log_lambda,
          # ibm_monteratorero_miroc6_sensitivity_growth[[4]]$log_lambda,
          # ibm_monteratorero_miroc6_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_mpi_esm1_2_lr_sensitivity_growth[[1]]$log_lambda,
          ibm_monteratorero_mpi_esm1_2_lr_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_monteratorero_mpi_esm1_2_lr_sensitivity_growth[[3]]$log_lambda,
          # ibm_monteratorero_mpi_esm1_2_lr_sensitivity_growth[[4]]$log_lambda,
          # ibm_monteratorero_mpi_esm1_2_lr_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_mri_esm2_0_sensitivity_growth[[1]]$log_lambda,
          ibm_monteratorero_mri_esm2_0_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_monteratorero_mri_esm2_0_sensitivity_growth[[3]]$log_lambda,
          # ibm_monteratorero_mri_esm2_0_sensitivity_growth[[4]]$log_lambda,
          # ibm_monteratorero_mri_esm2_0_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_noresm2_mm_sensitivity_growth[[1]]$log_lambda,
          ibm_monteratorero_noresm2_mm_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_monteratorero_noresm2_mm_sensitivity_growth[[3]]$log_lambda,
          # ibm_monteratorero_noresm2_mm_sensitivity_growth[[4]]$log_lambda,
          # ibm_monteratorero_noresm2_mm_sensitivity_growth[[5]]$log_lambda
  )))

rm(list = grep("growth", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))][grep("flowering", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))][grep("flowering", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "MonteraTorero" &
                               sensitivityanalysis_results$vital_rate == "flowering")] = 
  c(t(rbind(ibm_monteratorero_canesm5_sensitivity_flowering[[1]]$log_lambda,
            ibm_monteratorero_canesm5_sensitivity_flowering[[2]]$log_lambda
            # ,
            # ibm_monteratorero_canesm5_sensitivity_flowering[[3]]$log_lambda,
            # ibm_monteratorero_canesm5_sensitivity_flowering[[4]]$log_lambda,
            # ibm_monteratorero_canesm5_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_ec_earth3_sensitivity_flowering[[1]]$log_lambda,
          ibm_monteratorero_ec_earth3_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_monteratorero_ec_earth3_sensitivity_flowering[[3]]$log_lambda,
          # ibm_monteratorero_ec_earth3_sensitivity_flowering[[4]]$log_lambda,
          # ibm_monteratorero_ec_earth3_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_fgoals_g3_sensitivity_flowering[[1]]$log_lambda,
          ibm_monteratorero_fgoals_g3_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_monteratorero_fgoals_g3_sensitivity_flowering[[3]]$log_lambda,
          # ibm_monteratorero_fgoals_g3_sensitivity_flowering[[4]]$log_lambda,
          # ibm_monteratorero_fgoals_g3_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_gfdl_esm4_sensitivity_flowering[[1]]$log_lambda,
          ibm_monteratorero_gfdl_esm4_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_monteratorero_gfdl_esm4_sensitivity_flowering[[3]]$log_lambda,
          # ibm_monteratorero_gfdl_esm4_sensitivity_flowering[[4]]$log_lambda,
          # ibm_monteratorero_gfdl_esm4_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_giss_e2_1_g_sensitivity_flowering[[1]]$log_lambda,
          ibm_monteratorero_giss_e2_1_g_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_monteratorero_giss_e2_1_g_sensitivity_flowering[[3]]$log_lambda,
          # ibm_monteratorero_giss_e2_1_g_sensitivity_flowering[[4]]$log_lambda,
          # ibm_monteratorero_giss_e2_1_g_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_inm_cm4_8_sensitivity_flowering[[1]]$log_lambda,
          ibm_monteratorero_inm_cm4_8_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_monteratorero_inm_cm4_8_sensitivity_flowering[[3]]$log_lambda,
          # ibm_monteratorero_inm_cm4_8_sensitivity_flowering[[4]]$log_lambda,
          # ibm_monteratorero_inm_cm4_8_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_ipsl_cm6a_lr_sensitivity_flowering[[1]]$log_lambda,
          ibm_monteratorero_ipsl_cm6a_lr_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_monteratorero_ipsl_cm6a_lr_sensitivity_flowering[[3]]$log_lambda,
          # ibm_monteratorero_ipsl_cm6a_lr_sensitivity_flowering[[4]]$log_lambda,
          # ibm_monteratorero_ipsl_cm6a_lr_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_miroc6_sensitivity_flowering[[1]]$log_lambda,
          ibm_monteratorero_miroc6_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_monteratorero_miroc6_sensitivity_flowering[[3]]$log_lambda,
          # ibm_monteratorero_miroc6_sensitivity_flowering[[4]]$log_lambda,
          # ibm_monteratorero_miroc6_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_mpi_esm1_2_lr_sensitivity_flowering[[1]]$log_lambda,
          ibm_monteratorero_mpi_esm1_2_lr_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_monteratorero_mpi_esm1_2_lr_sensitivity_flowering[[3]]$log_lambda,
          # ibm_monteratorero_mpi_esm1_2_lr_sensitivity_flowering[[4]]$log_lambda,
          # ibm_monteratorero_mpi_esm1_2_lr_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_mri_esm2_0_sensitivity_flowering[[1]]$log_lambda,
          ibm_monteratorero_mri_esm2_0_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_monteratorero_mri_esm2_0_sensitivity_flowering[[3]]$log_lambda,
          # ibm_monteratorero_mri_esm2_0_sensitivity_flowering[[4]]$log_lambda,
          # ibm_monteratorero_mri_esm2_0_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_noresm2_mm_sensitivity_flowering[[1]]$log_lambda,
          ibm_monteratorero_noresm2_mm_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_monteratorero_noresm2_mm_sensitivity_flowering[[3]]$log_lambda,
          # ibm_monteratorero_noresm2_mm_sensitivity_flowering[[4]]$log_lambda,
          # ibm_monteratorero_noresm2_mm_sensitivity_flowering[[5]]$log_lambda
  )))

rm(list = grep("flowering", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))][grep("nbFlowers", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))][grep("nbFlowers", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "MonteraTorero" &
                               sensitivityanalysis_results$vital_rate == "nbFlowers")] = 
  c(t(rbind(ibm_monteratorero_canesm5_sensitivity_nbFlowers[[1]]$log_lambda,
            ibm_monteratorero_canesm5_sensitivity_nbFlowers[[2]]$log_lambda
            # ,
            # ibm_monteratorero_canesm5_sensitivity_nbFlowers[[3]]$log_lambda,
            # ibm_monteratorero_canesm5_sensitivity_nbFlowers[[4]]$log_lambda,
            # ibm_monteratorero_canesm5_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_ec_earth3_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_monteratorero_ec_earth3_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_monteratorero_ec_earth3_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_monteratorero_ec_earth3_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_monteratorero_ec_earth3_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_fgoals_g3_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_monteratorero_fgoals_g3_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_monteratorero_fgoals_g3_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_monteratorero_fgoals_g3_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_monteratorero_fgoals_g3_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_gfdl_esm4_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_monteratorero_gfdl_esm4_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_monteratorero_gfdl_esm4_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_monteratorero_gfdl_esm4_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_monteratorero_gfdl_esm4_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_giss_e2_1_g_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_monteratorero_giss_e2_1_g_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_monteratorero_giss_e2_1_g_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_monteratorero_giss_e2_1_g_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_monteratorero_giss_e2_1_g_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_inm_cm4_8_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_monteratorero_inm_cm4_8_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_monteratorero_inm_cm4_8_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_monteratorero_inm_cm4_8_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_monteratorero_inm_cm4_8_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_ipsl_cm6a_lr_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_monteratorero_ipsl_cm6a_lr_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_monteratorero_ipsl_cm6a_lr_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_monteratorero_ipsl_cm6a_lr_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_monteratorero_ipsl_cm6a_lr_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_miroc6_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_monteratorero_miroc6_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_monteratorero_miroc6_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_monteratorero_miroc6_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_monteratorero_miroc6_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_mpi_esm1_2_lr_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_monteratorero_mpi_esm1_2_lr_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_monteratorero_mpi_esm1_2_lr_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_monteratorero_mpi_esm1_2_lr_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_monteratorero_mpi_esm1_2_lr_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_mri_esm2_0_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_monteratorero_mri_esm2_0_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_monteratorero_mri_esm2_0_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_monteratorero_mri_esm2_0_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_monteratorero_mri_esm2_0_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_noresm2_mm_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_monteratorero_noresm2_mm_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_monteratorero_noresm2_mm_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_monteratorero_noresm2_mm_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_monteratorero_noresm2_mm_sensitivity_nbFlowers[[5]]$log_lambda
  )))

rm(list = grep("nbFlowers", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))][grep("seedlingSize", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))][grep("seedlingSize", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_MonteraTorero", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "MonteraTorero" &
                               sensitivityanalysis_results$vital_rate == "seedlingSize")] = 
  c(t(rbind(ibm_monteratorero_canesm5_sensitivity_seedlingSize[[1]]$log_lambda,
            ibm_monteratorero_canesm5_sensitivity_seedlingSize[[2]]$log_lambda
            # ,
            # ibm_monteratorero_canesm5_sensitivity_seedlingSize[[3]]$log_lambda,
            # ibm_monteratorero_canesm5_sensitivity_seedlingSize[[4]]$log_lambda,
            # ibm_monteratorero_canesm5_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_ec_earth3_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_monteratorero_ec_earth3_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_monteratorero_ec_earth3_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_monteratorero_ec_earth3_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_monteratorero_ec_earth3_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_fgoals_g3_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_monteratorero_fgoals_g3_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_monteratorero_fgoals_g3_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_monteratorero_fgoals_g3_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_monteratorero_fgoals_g3_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_gfdl_esm4_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_monteratorero_gfdl_esm4_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_monteratorero_gfdl_esm4_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_monteratorero_gfdl_esm4_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_monteratorero_gfdl_esm4_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_giss_e2_1_g_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_monteratorero_giss_e2_1_g_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_monteratorero_giss_e2_1_g_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_monteratorero_giss_e2_1_g_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_monteratorero_giss_e2_1_g_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_inm_cm4_8_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_monteratorero_inm_cm4_8_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_monteratorero_inm_cm4_8_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_monteratorero_inm_cm4_8_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_monteratorero_inm_cm4_8_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_ipsl_cm6a_lr_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_monteratorero_ipsl_cm6a_lr_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_monteratorero_ipsl_cm6a_lr_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_monteratorero_ipsl_cm6a_lr_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_monteratorero_ipsl_cm6a_lr_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_miroc6_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_monteratorero_miroc6_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_monteratorero_miroc6_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_monteratorero_miroc6_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_monteratorero_miroc6_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_mpi_esm1_2_lr_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_monteratorero_mpi_esm1_2_lr_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_monteratorero_mpi_esm1_2_lr_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_monteratorero_mpi_esm1_2_lr_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_monteratorero_mpi_esm1_2_lr_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_mri_esm2_0_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_monteratorero_mri_esm2_0_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_monteratorero_mri_esm2_0_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_monteratorero_mri_esm2_0_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_monteratorero_mri_esm2_0_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_monteratorero_noresm2_mm_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_monteratorero_noresm2_mm_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_monteratorero_noresm2_mm_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_monteratorero_noresm2_mm_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_monteratorero_noresm2_mm_sensitivity_seedlingSize[[5]]$log_lambda
  )))

rm(list = grep("seedlingSize", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))][grep("survival", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))][grep("survival", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "SCarbDist" &
                               sensitivityanalysis_results$vital_rate == "survival")] = 
  c(t(rbind(ibm_scarbdist_canesm5_sensitivity_survival[[1]]$log_lambda,
            ibm_scarbdist_canesm5_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_scarbdist_canesm5_sensitivity_survival[[3]]$log_lambda,
            # ibm_scarbdist_canesm5_sensitivity_survival[[4]]$log_lambda,
            # ibm_scarbdist_canesm5_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_ec_earth3_sensitivity_survival[[1]]$log_lambda,
          ibm_scarbdist_ec_earth3_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_scarbdist_ec_earth3_sensitivity_survival[[3]]$log_lambda,
          # ibm_scarbdist_ec_earth3_sensitivity_survival[[4]]$log_lambda,
          # ibm_scarbdist_ec_earth3_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_fgoals_g3_sensitivity_survival[[1]]$log_lambda,
          ibm_scarbdist_fgoals_g3_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_scarbdist_fgoals_g3_sensitivity_survival[[3]]$log_lambda,
          # ibm_scarbdist_fgoals_g3_sensitivity_survival[[4]]$log_lambda,
          # ibm_scarbdist_fgoals_g3_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_gfdl_esm4_sensitivity_survival[[1]]$log_lambda,
          ibm_scarbdist_gfdl_esm4_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_scarbdist_gfdl_esm4_sensitivity_survival[[3]]$log_lambda,
          # ibm_scarbdist_gfdl_esm4_sensitivity_survival[[4]]$log_lambda,
          # ibm_scarbdist_gfdl_esm4_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_giss_e2_1_g_sensitivity_survival[[1]]$log_lambda,
          ibm_scarbdist_giss_e2_1_g_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_scarbdist_giss_e2_1_g_sensitivity_survival[[3]]$log_lambda,
          # ibm_scarbdist_giss_e2_1_g_sensitivity_survival[[4]]$log_lambda,
          # ibm_scarbdist_giss_e2_1_g_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_inm_cm4_8_sensitivity_survival[[1]]$log_lambda,
          ibm_scarbdist_inm_cm4_8_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_scarbdist_inm_cm4_8_sensitivity_survival[[3]]$log_lambda,
          # ibm_scarbdist_inm_cm4_8_sensitivity_survival[[4]]$log_lambda,
          # ibm_scarbdist_inm_cm4_8_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_ipsl_cm6a_lr_sensitivity_survival[[1]]$log_lambda,
          ibm_scarbdist_ipsl_cm6a_lr_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_scarbdist_ipsl_cm6a_lr_sensitivity_survival[[3]]$log_lambda,
          # ibm_scarbdist_ipsl_cm6a_lr_sensitivity_survival[[4]]$log_lambda,
          # ibm_scarbdist_ipsl_cm6a_lr_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_miroc6_sensitivity_survival[[1]]$log_lambda,
          ibm_scarbdist_miroc6_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_scarbdist_miroc6_sensitivity_survival[[3]]$log_lambda,
          # ibm_scarbdist_miroc6_sensitivity_survival[[4]]$log_lambda,
          # ibm_scarbdist_miroc6_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_mpi_esm1_2_lr_sensitivity_survival[[1]]$log_lambda,
          ibm_scarbdist_mpi_esm1_2_lr_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_scarbdist_mpi_esm1_2_lr_sensitivity_survival[[3]]$log_lambda,
          # ibm_scarbdist_mpi_esm1_2_lr_sensitivity_survival[[4]]$log_lambda,
          # ibm_scarbdist_mpi_esm1_2_lr_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_mri_esm2_0_sensitivity_survival[[1]]$log_lambda,
          ibm_scarbdist_mri_esm2_0_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_scarbdist_mri_esm2_0_sensitivity_survival[[3]]$log_lambda,
          # ibm_scarbdist_mri_esm2_0_sensitivity_survival[[4]]$log_lambda,
          # ibm_scarbdist_mri_esm2_0_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_noresm2_mm_sensitivity_survival[[1]]$log_lambda,
          ibm_scarbdist_noresm2_mm_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_scarbdist_noresm2_mm_sensitivity_survival[[3]]$log_lambda,
          # ibm_scarbdist_noresm2_mm_sensitivity_survival[[4]]$log_lambda,
          # ibm_scarbdist_noresm2_mm_sensitivity_survival[[5]]$log_lambda
  )))

rm(list = grep("survival", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))][grep("growth", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))][grep("growth", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "SCarbDist" &
                               sensitivityanalysis_results$vital_rate == "growth")] = 
  c(t(rbind(ibm_scarbdist_canesm5_sensitivity_growth[[1]]$log_lambda,
            ibm_scarbdist_canesm5_sensitivity_growth[[2]]$log_lambda
            # ,
            # ibm_scarbdist_canesm5_sensitivity_growth[[3]]$log_lambda,
            # ibm_scarbdist_canesm5_sensitivity_growth[[4]]$log_lambda,
            # ibm_scarbdist_canesm5_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_ec_earth3_sensitivity_growth[[1]]$log_lambda,
          ibm_scarbdist_ec_earth3_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_scarbdist_ec_earth3_sensitivity_growth[[3]]$log_lambda,
          # ibm_scarbdist_ec_earth3_sensitivity_growth[[4]]$log_lambda,
          # ibm_scarbdist_ec_earth3_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_fgoals_g3_sensitivity_growth[[1]]$log_lambda,
          ibm_scarbdist_fgoals_g3_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_scarbdist_fgoals_g3_sensitivity_growth[[3]]$log_lambda,
          # ibm_scarbdist_fgoals_g3_sensitivity_growth[[4]]$log_lambda,
          # ibm_scarbdist_fgoals_g3_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_gfdl_esm4_sensitivity_growth[[1]]$log_lambda,
          ibm_scarbdist_gfdl_esm4_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_scarbdist_gfdl_esm4_sensitivity_growth[[3]]$log_lambda,
          # ibm_scarbdist_gfdl_esm4_sensitivity_growth[[4]]$log_lambda,
          # ibm_scarbdist_gfdl_esm4_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_giss_e2_1_g_sensitivity_growth[[1]]$log_lambda,
          ibm_scarbdist_giss_e2_1_g_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_scarbdist_giss_e2_1_g_sensitivity_growth[[3]]$log_lambda,
          # ibm_scarbdist_giss_e2_1_g_sensitivity_growth[[4]]$log_lambda,
          # ibm_scarbdist_giss_e2_1_g_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_inm_cm4_8_sensitivity_growth[[1]]$log_lambda,
          ibm_scarbdist_inm_cm4_8_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_scarbdist_inm_cm4_8_sensitivity_growth[[3]]$log_lambda,
          # ibm_scarbdist_inm_cm4_8_sensitivity_growth[[4]]$log_lambda,
          # ibm_scarbdist_inm_cm4_8_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_ipsl_cm6a_lr_sensitivity_growth[[1]]$log_lambda,
          ibm_scarbdist_ipsl_cm6a_lr_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_scarbdist_ipsl_cm6a_lr_sensitivity_growth[[3]]$log_lambda,
          # ibm_scarbdist_ipsl_cm6a_lr_sensitivity_growth[[4]]$log_lambda,
          # ibm_scarbdist_ipsl_cm6a_lr_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_miroc6_sensitivity_growth[[1]]$log_lambda,
          ibm_scarbdist_miroc6_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_scarbdist_miroc6_sensitivity_growth[[3]]$log_lambda,
          # ibm_scarbdist_miroc6_sensitivity_growth[[4]]$log_lambda,
          # ibm_scarbdist_miroc6_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_mpi_esm1_2_lr_sensitivity_growth[[1]]$log_lambda,
          ibm_scarbdist_mpi_esm1_2_lr_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_scarbdist_mpi_esm1_2_lr_sensitivity_growth[[3]]$log_lambda,
          # ibm_scarbdist_mpi_esm1_2_lr_sensitivity_growth[[4]]$log_lambda,
          # ibm_scarbdist_mpi_esm1_2_lr_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_mri_esm2_0_sensitivity_growth[[1]]$log_lambda,
          ibm_scarbdist_mri_esm2_0_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_scarbdist_mri_esm2_0_sensitivity_growth[[3]]$log_lambda,
          # ibm_scarbdist_mri_esm2_0_sensitivity_growth[[4]]$log_lambda,
          # ibm_scarbdist_mri_esm2_0_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_noresm2_mm_sensitivity_growth[[1]]$log_lambda,
          ibm_scarbdist_noresm2_mm_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_scarbdist_noresm2_mm_sensitivity_growth[[3]]$log_lambda,
          # ibm_scarbdist_noresm2_mm_sensitivity_growth[[4]]$log_lambda,
          # ibm_scarbdist_noresm2_mm_sensitivity_growth[[5]]$log_lambda
  )))

rm(list = grep("growth", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))][grep("flowering", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))][grep("flowering", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "SCarbDist" &
                               sensitivityanalysis_results$vital_rate == "flowering")] = 
  c(t(rbind(ibm_scarbdist_canesm5_sensitivity_flowering[[1]]$log_lambda,
            ibm_scarbdist_canesm5_sensitivity_flowering[[2]]$log_lambda
            # ,
            # ibm_scarbdist_canesm5_sensitivity_flowering[[3]]$log_lambda,
            # ibm_scarbdist_canesm5_sensitivity_flowering[[4]]$log_lambda,
            # ibm_scarbdist_canesm5_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_ec_earth3_sensitivity_flowering[[1]]$log_lambda,
          ibm_scarbdist_ec_earth3_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_scarbdist_ec_earth3_sensitivity_flowering[[3]]$log_lambda,
          # ibm_scarbdist_ec_earth3_sensitivity_flowering[[4]]$log_lambda,
          # ibm_scarbdist_ec_earth3_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_fgoals_g3_sensitivity_flowering[[1]]$log_lambda,
          ibm_scarbdist_fgoals_g3_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_scarbdist_fgoals_g3_sensitivity_flowering[[3]]$log_lambda,
          # ibm_scarbdist_fgoals_g3_sensitivity_flowering[[4]]$log_lambda,
          # ibm_scarbdist_fgoals_g3_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_gfdl_esm4_sensitivity_flowering[[1]]$log_lambda,
          ibm_scarbdist_gfdl_esm4_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_scarbdist_gfdl_esm4_sensitivity_flowering[[3]]$log_lambda,
          # ibm_scarbdist_gfdl_esm4_sensitivity_flowering[[4]]$log_lambda,
          # ibm_scarbdist_gfdl_esm4_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_giss_e2_1_g_sensitivity_flowering[[1]]$log_lambda,
          ibm_scarbdist_giss_e2_1_g_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_scarbdist_giss_e2_1_g_sensitivity_flowering[[3]]$log_lambda,
          # ibm_scarbdist_giss_e2_1_g_sensitivity_flowering[[4]]$log_lambda,
          # ibm_scarbdist_giss_e2_1_g_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_inm_cm4_8_sensitivity_flowering[[1]]$log_lambda,
          ibm_scarbdist_inm_cm4_8_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_scarbdist_inm_cm4_8_sensitivity_flowering[[3]]$log_lambda,
          # ibm_scarbdist_inm_cm4_8_sensitivity_flowering[[4]]$log_lambda,
          # ibm_scarbdist_inm_cm4_8_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_ipsl_cm6a_lr_sensitivity_flowering[[1]]$log_lambda,
          ibm_scarbdist_ipsl_cm6a_lr_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_scarbdist_ipsl_cm6a_lr_sensitivity_flowering[[3]]$log_lambda,
          # ibm_scarbdist_ipsl_cm6a_lr_sensitivity_flowering[[4]]$log_lambda,
          # ibm_scarbdist_ipsl_cm6a_lr_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_miroc6_sensitivity_flowering[[1]]$log_lambda,
          ibm_scarbdist_miroc6_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_scarbdist_miroc6_sensitivity_flowering[[3]]$log_lambda,
          # ibm_scarbdist_miroc6_sensitivity_flowering[[4]]$log_lambda,
          # ibm_scarbdist_miroc6_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_mpi_esm1_2_lr_sensitivity_flowering[[1]]$log_lambda,
          ibm_scarbdist_mpi_esm1_2_lr_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_scarbdist_mpi_esm1_2_lr_sensitivity_flowering[[3]]$log_lambda,
          # ibm_scarbdist_mpi_esm1_2_lr_sensitivity_flowering[[4]]$log_lambda,
          # ibm_scarbdist_mpi_esm1_2_lr_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_mri_esm2_0_sensitivity_flowering[[1]]$log_lambda,
          ibm_scarbdist_mri_esm2_0_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_scarbdist_mri_esm2_0_sensitivity_flowering[[3]]$log_lambda,
          # ibm_scarbdist_mri_esm2_0_sensitivity_flowering[[4]]$log_lambda,
          # ibm_scarbdist_mri_esm2_0_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_noresm2_mm_sensitivity_flowering[[1]]$log_lambda,
          ibm_scarbdist_noresm2_mm_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_scarbdist_noresm2_mm_sensitivity_flowering[[3]]$log_lambda,
          # ibm_scarbdist_noresm2_mm_sensitivity_flowering[[4]]$log_lambda,
          # ibm_scarbdist_noresm2_mm_sensitivity_flowering[[5]]$log_lambda
  )))

rm(list = grep("flowering", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))][grep("nbFlowers", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))][grep("nbFlowers", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "SCarbDist" &
                               sensitivityanalysis_results$vital_rate == "nbFlowers")] = 
  c(t(rbind(ibm_scarbdist_canesm5_sensitivity_nbFlowers[[1]]$log_lambda,
            ibm_scarbdist_canesm5_sensitivity_nbFlowers[[2]]$log_lambda
            # ,
            # ibm_scarbdist_canesm5_sensitivity_nbFlowers[[3]]$log_lambda,
            # ibm_scarbdist_canesm5_sensitivity_nbFlowers[[4]]$log_lambda,
            # ibm_scarbdist_canesm5_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_ec_earth3_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_scarbdist_ec_earth3_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_scarbdist_ec_earth3_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_scarbdist_ec_earth3_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_scarbdist_ec_earth3_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_fgoals_g3_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_scarbdist_fgoals_g3_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_scarbdist_fgoals_g3_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_scarbdist_fgoals_g3_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_scarbdist_fgoals_g3_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_gfdl_esm4_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_scarbdist_gfdl_esm4_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_scarbdist_gfdl_esm4_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_scarbdist_gfdl_esm4_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_scarbdist_gfdl_esm4_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_giss_e2_1_g_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_scarbdist_giss_e2_1_g_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_scarbdist_giss_e2_1_g_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_scarbdist_giss_e2_1_g_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_scarbdist_giss_e2_1_g_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_inm_cm4_8_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_scarbdist_inm_cm4_8_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_scarbdist_inm_cm4_8_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_scarbdist_inm_cm4_8_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_scarbdist_inm_cm4_8_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_ipsl_cm6a_lr_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_scarbdist_ipsl_cm6a_lr_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_scarbdist_ipsl_cm6a_lr_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_scarbdist_ipsl_cm6a_lr_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_scarbdist_ipsl_cm6a_lr_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_miroc6_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_scarbdist_miroc6_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_scarbdist_miroc6_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_scarbdist_miroc6_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_scarbdist_miroc6_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_mpi_esm1_2_lr_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_scarbdist_mpi_esm1_2_lr_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_scarbdist_mpi_esm1_2_lr_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_scarbdist_mpi_esm1_2_lr_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_scarbdist_mpi_esm1_2_lr_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_mri_esm2_0_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_scarbdist_mri_esm2_0_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_scarbdist_mri_esm2_0_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_scarbdist_mri_esm2_0_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_scarbdist_mri_esm2_0_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_noresm2_mm_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_scarbdist_noresm2_mm_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_scarbdist_noresm2_mm_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_scarbdist_noresm2_mm_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_scarbdist_noresm2_mm_sensitivity_nbFlowers[[5]]$log_lambda
  )))

rm(list = grep("nbFlowers", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))][grep("seedlingSize", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))][grep("seedlingSize", list.files("Output/Projections/")[grep("Anthropogenic_SensitivityAnalysis_SCarbDist", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "SCarbDist" &
                               sensitivityanalysis_results$vital_rate == "seedlingSize")] = 
  c(t(rbind(ibm_scarbdist_canesm5_sensitivity_seedlingSize[[1]]$log_lambda,
            ibm_scarbdist_canesm5_sensitivity_seedlingSize[[2]]$log_lambda
            # ,
            # ibm_scarbdist_canesm5_sensitivity_seedlingSize[[3]]$log_lambda,
            # ibm_scarbdist_canesm5_sensitivity_seedlingSize[[4]]$log_lambda,
            # ibm_scarbdist_canesm5_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_ec_earth3_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_scarbdist_ec_earth3_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_scarbdist_ec_earth3_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_scarbdist_ec_earth3_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_scarbdist_ec_earth3_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_fgoals_g3_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_scarbdist_fgoals_g3_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_scarbdist_fgoals_g3_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_scarbdist_fgoals_g3_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_scarbdist_fgoals_g3_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_gfdl_esm4_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_scarbdist_gfdl_esm4_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_scarbdist_gfdl_esm4_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_scarbdist_gfdl_esm4_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_scarbdist_gfdl_esm4_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_giss_e2_1_g_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_scarbdist_giss_e2_1_g_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_scarbdist_giss_e2_1_g_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_scarbdist_giss_e2_1_g_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_scarbdist_giss_e2_1_g_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_inm_cm4_8_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_scarbdist_inm_cm4_8_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_scarbdist_inm_cm4_8_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_scarbdist_inm_cm4_8_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_scarbdist_inm_cm4_8_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_ipsl_cm6a_lr_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_scarbdist_ipsl_cm6a_lr_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_scarbdist_ipsl_cm6a_lr_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_scarbdist_ipsl_cm6a_lr_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_scarbdist_ipsl_cm6a_lr_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_miroc6_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_scarbdist_miroc6_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_scarbdist_miroc6_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_scarbdist_miroc6_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_scarbdist_miroc6_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_mpi_esm1_2_lr_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_scarbdist_mpi_esm1_2_lr_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_scarbdist_mpi_esm1_2_lr_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_scarbdist_mpi_esm1_2_lr_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_scarbdist_mpi_esm1_2_lr_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_mri_esm2_0_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_scarbdist_mri_esm2_0_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_scarbdist_mri_esm2_0_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_scarbdist_mri_esm2_0_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_scarbdist_mri_esm2_0_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_scarbdist_noresm2_mm_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_scarbdist_noresm2_mm_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_scarbdist_noresm2_mm_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_scarbdist_noresm2_mm_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_scarbdist_noresm2_mm_sensitivity_seedlingSize[[5]]$log_lambda
  )))

rm(list = grep("seedlingSize", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))][grep("survival", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))][grep("survival", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "SierraCarboneraY5" &
                               sensitivityanalysis_results$vital_rate == "survival")] = 
  c(t(rbind(ibm_sierracarboneray5_canesm5_sensitivity_survival[[1]]$log_lambda,
            ibm_sierracarboneray5_canesm5_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_sierracarboneray5_canesm5_sensitivity_survival[[3]]$log_lambda,
            # ibm_sierracarboneray5_canesm5_sensitivity_survival[[4]]$log_lambda,
            # ibm_sierracarboneray5_canesm5_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_ec_earth3_sensitivity_survival[[1]]$log_lambda,
          ibm_sierracarboneray5_ec_earth3_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_ec_earth3_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierracarboneray5_ec_earth3_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierracarboneray5_ec_earth3_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_fgoals_g3_sensitivity_survival[[1]]$log_lambda,
          ibm_sierracarboneray5_fgoals_g3_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_fgoals_g3_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierracarboneray5_fgoals_g3_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierracarboneray5_fgoals_g3_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_gfdl_esm4_sensitivity_survival[[1]]$log_lambda,
          ibm_sierracarboneray5_gfdl_esm4_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_gfdl_esm4_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierracarboneray5_gfdl_esm4_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierracarboneray5_gfdl_esm4_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_giss_e2_1_g_sensitivity_survival[[1]]$log_lambda,
          ibm_sierracarboneray5_giss_e2_1_g_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_giss_e2_1_g_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierracarboneray5_giss_e2_1_g_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierracarboneray5_giss_e2_1_g_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_inm_cm4_8_sensitivity_survival[[1]]$log_lambda,
          ibm_sierracarboneray5_inm_cm4_8_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_inm_cm4_8_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierracarboneray5_inm_cm4_8_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierracarboneray5_inm_cm4_8_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_survival[[1]]$log_lambda,
          ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_miroc6_sensitivity_survival[[1]]$log_lambda,
          ibm_sierracarboneray5_miroc6_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_miroc6_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierracarboneray5_miroc6_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierracarboneray5_miroc6_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_survival[[1]]$log_lambda,
          ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_mri_esm2_0_sensitivity_survival[[1]]$log_lambda,
          ibm_sierracarboneray5_mri_esm2_0_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_mri_esm2_0_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierracarboneray5_mri_esm2_0_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierracarboneray5_mri_esm2_0_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_noresm2_mm_sensitivity_survival[[1]]$log_lambda,
          ibm_sierracarboneray5_noresm2_mm_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_noresm2_mm_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierracarboneray5_noresm2_mm_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierracarboneray5_noresm2_mm_sensitivity_survival[[5]]$log_lambda
  )))

rm(list = grep("survival", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))][grep("growth", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))][grep("growth", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "SierraCarboneraY5" &
                               sensitivityanalysis_results$vital_rate == "growth")] = 
  c(t(rbind(ibm_sierracarboneray5_canesm5_sensitivity_growth[[1]]$log_lambda,
            ibm_sierracarboneray5_canesm5_sensitivity_growth[[2]]$log_lambda
            # ,
            # ibm_sierracarboneray5_canesm5_sensitivity_growth[[3]]$log_lambda,
            # ibm_sierracarboneray5_canesm5_sensitivity_growth[[4]]$log_lambda,
            # ibm_sierracarboneray5_canesm5_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_ec_earth3_sensitivity_growth[[1]]$log_lambda,
          ibm_sierracarboneray5_ec_earth3_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_ec_earth3_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierracarboneray5_ec_earth3_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierracarboneray5_ec_earth3_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_fgoals_g3_sensitivity_growth[[1]]$log_lambda,
          ibm_sierracarboneray5_fgoals_g3_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_fgoals_g3_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierracarboneray5_fgoals_g3_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierracarboneray5_fgoals_g3_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_gfdl_esm4_sensitivity_growth[[1]]$log_lambda,
          ibm_sierracarboneray5_gfdl_esm4_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_gfdl_esm4_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierracarboneray5_gfdl_esm4_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierracarboneray5_gfdl_esm4_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_giss_e2_1_g_sensitivity_growth[[1]]$log_lambda,
          ibm_sierracarboneray5_giss_e2_1_g_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_giss_e2_1_g_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierracarboneray5_giss_e2_1_g_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierracarboneray5_giss_e2_1_g_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_inm_cm4_8_sensitivity_growth[[1]]$log_lambda,
          ibm_sierracarboneray5_inm_cm4_8_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_inm_cm4_8_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierracarboneray5_inm_cm4_8_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierracarboneray5_inm_cm4_8_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_growth[[1]]$log_lambda,
          ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_miroc6_sensitivity_growth[[1]]$log_lambda,
          ibm_sierracarboneray5_miroc6_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_miroc6_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierracarboneray5_miroc6_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierracarboneray5_miroc6_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_growth[[1]]$log_lambda,
          ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_mri_esm2_0_sensitivity_growth[[1]]$log_lambda,
          ibm_sierracarboneray5_mri_esm2_0_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_mri_esm2_0_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierracarboneray5_mri_esm2_0_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierracarboneray5_mri_esm2_0_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_noresm2_mm_sensitivity_growth[[1]]$log_lambda,
          ibm_sierracarboneray5_noresm2_mm_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_noresm2_mm_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierracarboneray5_noresm2_mm_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierracarboneray5_noresm2_mm_sensitivity_growth[[5]]$log_lambda
  )))

rm(list = grep("growth", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))][grep("flowering", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))][grep("flowering", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "SierraCarboneraY5" &
                               sensitivityanalysis_results$vital_rate == "flowering")] = 
  c(t(rbind(ibm_sierracarboneray5_canesm5_sensitivity_flowering[[1]]$log_lambda,
            ibm_sierracarboneray5_canesm5_sensitivity_flowering[[2]]$log_lambda
            # ,
            # ibm_sierracarboneray5_canesm5_sensitivity_flowering[[3]]$log_lambda,
            # ibm_sierracarboneray5_canesm5_sensitivity_flowering[[4]]$log_lambda,
            # ibm_sierracarboneray5_canesm5_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_ec_earth3_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierracarboneray5_ec_earth3_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_ec_earth3_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierracarboneray5_ec_earth3_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierracarboneray5_ec_earth3_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_fgoals_g3_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierracarboneray5_fgoals_g3_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_fgoals_g3_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierracarboneray5_fgoals_g3_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierracarboneray5_fgoals_g3_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_gfdl_esm4_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierracarboneray5_gfdl_esm4_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_gfdl_esm4_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierracarboneray5_gfdl_esm4_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierracarboneray5_gfdl_esm4_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_giss_e2_1_g_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierracarboneray5_giss_e2_1_g_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_giss_e2_1_g_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierracarboneray5_giss_e2_1_g_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierracarboneray5_giss_e2_1_g_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_inm_cm4_8_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierracarboneray5_inm_cm4_8_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_inm_cm4_8_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierracarboneray5_inm_cm4_8_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierracarboneray5_inm_cm4_8_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_miroc6_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierracarboneray5_miroc6_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_miroc6_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierracarboneray5_miroc6_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierracarboneray5_miroc6_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_mri_esm2_0_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierracarboneray5_mri_esm2_0_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_mri_esm2_0_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierracarboneray5_mri_esm2_0_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierracarboneray5_mri_esm2_0_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_noresm2_mm_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierracarboneray5_noresm2_mm_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_noresm2_mm_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierracarboneray5_noresm2_mm_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierracarboneray5_noresm2_mm_sensitivity_flowering[[5]]$log_lambda
  )))

rm(list = grep("flowering", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))][grep("nbFlowers", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))][grep("nbFlowers", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "SierraCarboneraY5" &
                               sensitivityanalysis_results$vital_rate == "nbFlowers")] = 
  c(t(rbind(ibm_sierracarboneray5_canesm5_sensitivity_nbFlowers[[1]]$log_lambda,
            ibm_sierracarboneray5_canesm5_sensitivity_nbFlowers[[2]]$log_lambda
            # ,
            # ibm_sierracarboneray5_canesm5_sensitivity_nbFlowers[[3]]$log_lambda,
            # ibm_sierracarboneray5_canesm5_sensitivity_nbFlowers[[4]]$log_lambda,
            # ibm_sierracarboneray5_canesm5_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_ec_earth3_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierracarboneray5_ec_earth3_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_ec_earth3_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierracarboneray5_ec_earth3_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierracarboneray5_ec_earth3_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_fgoals_g3_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierracarboneray5_fgoals_g3_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_fgoals_g3_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierracarboneray5_fgoals_g3_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierracarboneray5_fgoals_g3_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_gfdl_esm4_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierracarboneray5_gfdl_esm4_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_gfdl_esm4_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierracarboneray5_gfdl_esm4_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierracarboneray5_gfdl_esm4_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_giss_e2_1_g_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierracarboneray5_giss_e2_1_g_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_giss_e2_1_g_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierracarboneray5_giss_e2_1_g_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierracarboneray5_giss_e2_1_g_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_inm_cm4_8_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierracarboneray5_inm_cm4_8_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_inm_cm4_8_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierracarboneray5_inm_cm4_8_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierracarboneray5_inm_cm4_8_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_miroc6_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierracarboneray5_miroc6_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_miroc6_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierracarboneray5_miroc6_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierracarboneray5_miroc6_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_mri_esm2_0_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierracarboneray5_mri_esm2_0_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_mri_esm2_0_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierracarboneray5_mri_esm2_0_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierracarboneray5_mri_esm2_0_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_noresm2_mm_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierracarboneray5_noresm2_mm_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_noresm2_mm_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierracarboneray5_noresm2_mm_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierracarboneray5_noresm2_mm_sensitivity_nbFlowers[[5]]$log_lambda
  )))

rm(list = grep("nbFlowers", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))][grep("seedlingSize", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))][grep("seedlingSize", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraCarboneraY5", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "SierraCarboneraY5" &
                               sensitivityanalysis_results$vital_rate == "seedlingSize")] = 
  c(t(rbind(ibm_sierracarboneray5_canesm5_sensitivity_seedlingSize[[1]]$log_lambda,
            ibm_sierracarboneray5_canesm5_sensitivity_seedlingSize[[2]]$log_lambda
            # ,
            # ibm_sierracarboneray5_canesm5_sensitivity_seedlingSize[[3]]$log_lambda,
            # ibm_sierracarboneray5_canesm5_sensitivity_seedlingSize[[4]]$log_lambda,
            # ibm_sierracarboneray5_canesm5_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_ec_earth3_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierracarboneray5_ec_earth3_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_ec_earth3_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierracarboneray5_ec_earth3_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierracarboneray5_ec_earth3_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_fgoals_g3_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierracarboneray5_fgoals_g3_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_fgoals_g3_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierracarboneray5_fgoals_g3_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierracarboneray5_fgoals_g3_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_gfdl_esm4_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierracarboneray5_gfdl_esm4_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_gfdl_esm4_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierracarboneray5_gfdl_esm4_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierracarboneray5_gfdl_esm4_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_giss_e2_1_g_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierracarboneray5_giss_e2_1_g_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_giss_e2_1_g_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierracarboneray5_giss_e2_1_g_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierracarboneray5_giss_e2_1_g_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_inm_cm4_8_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierracarboneray5_inm_cm4_8_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_inm_cm4_8_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierracarboneray5_inm_cm4_8_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierracarboneray5_inm_cm4_8_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierracarboneray5_ipsl_cm6a_lr_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_miroc6_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierracarboneray5_miroc6_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_miroc6_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierracarboneray5_miroc6_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierracarboneray5_miroc6_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierracarboneray5_mpi_esm1_2_lr_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_mri_esm2_0_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierracarboneray5_mri_esm2_0_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_mri_esm2_0_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierracarboneray5_mri_esm2_0_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierracarboneray5_mri_esm2_0_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierracarboneray5_noresm2_mm_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierracarboneray5_noresm2_mm_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierracarboneray5_noresm2_mm_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierracarboneray5_noresm2_mm_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierracarboneray5_noresm2_mm_sensitivity_seedlingSize[[5]]$log_lambda
  )))

rm(list = grep("seedlingSize", ls() , value = TRUE, invert = FALSE))



for(i in 1:length(list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))][grep("survival", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))][grep("survival", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "SierraRetinY5" &
                               sensitivityanalysis_results$vital_rate == "survival")] = 
  c(t(rbind(ibm_sierraretiny5_canesm5_sensitivity_survival[[1]]$log_lambda,
            ibm_sierraretiny5_canesm5_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_sierraretiny5_canesm5_sensitivity_survival[[3]]$log_lambda,
            # ibm_sierraretiny5_canesm5_sensitivity_survival[[4]]$log_lambda,
            # ibm_sierraretiny5_canesm5_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_ec_earth3_sensitivity_survival[[1]]$log_lambda,
          ibm_sierraretiny5_ec_earth3_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_ec_earth3_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierraretiny5_ec_earth3_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierraretiny5_ec_earth3_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_fgoals_g3_sensitivity_survival[[1]]$log_lambda,
          ibm_sierraretiny5_fgoals_g3_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_fgoals_g3_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierraretiny5_fgoals_g3_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierraretiny5_fgoals_g3_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_gfdl_esm4_sensitivity_survival[[1]]$log_lambda,
          ibm_sierraretiny5_gfdl_esm4_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_gfdl_esm4_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierraretiny5_gfdl_esm4_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierraretiny5_gfdl_esm4_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_giss_e2_1_g_sensitivity_survival[[1]]$log_lambda,
          ibm_sierraretiny5_giss_e2_1_g_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_giss_e2_1_g_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierraretiny5_giss_e2_1_g_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierraretiny5_giss_e2_1_g_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_inm_cm4_8_sensitivity_survival[[1]]$log_lambda,
          ibm_sierraretiny5_inm_cm4_8_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_inm_cm4_8_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierraretiny5_inm_cm4_8_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierraretiny5_inm_cm4_8_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_survival[[1]]$log_lambda,
          ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_miroc6_sensitivity_survival[[1]]$log_lambda,
          ibm_sierraretiny5_miroc6_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_miroc6_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierraretiny5_miroc6_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierraretiny5_miroc6_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_survival[[1]]$log_lambda,
          ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_mri_esm2_0_sensitivity_survival[[1]]$log_lambda,
          ibm_sierraretiny5_mri_esm2_0_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_mri_esm2_0_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierraretiny5_mri_esm2_0_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierraretiny5_mri_esm2_0_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_noresm2_mm_sensitivity_survival[[1]]$log_lambda,
          ibm_sierraretiny5_noresm2_mm_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_noresm2_mm_sensitivity_survival[[3]]$log_lambda,
          # ibm_sierraretiny5_noresm2_mm_sensitivity_survival[[4]]$log_lambda,
          # ibm_sierraretiny5_noresm2_mm_sensitivity_survival[[5]]$log_lambda
  )))

rm(list = grep("survival", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))][grep("growth", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))][grep("growth", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "SierraRetinY5" &
                               sensitivityanalysis_results$vital_rate == "growth")] = 
  c(t(rbind(ibm_sierraretiny5_canesm5_sensitivity_growth[[1]]$log_lambda,
            ibm_sierraretiny5_canesm5_sensitivity_growth[[2]]$log_lambda
            # ,
            # ibm_sierraretiny5_canesm5_sensitivity_growth[[3]]$log_lambda,
            # ibm_sierraretiny5_canesm5_sensitivity_growth[[4]]$log_lambda,
            # ibm_sierraretiny5_canesm5_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_ec_earth3_sensitivity_growth[[1]]$log_lambda,
          ibm_sierraretiny5_ec_earth3_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_ec_earth3_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierraretiny5_ec_earth3_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierraretiny5_ec_earth3_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_fgoals_g3_sensitivity_growth[[1]]$log_lambda,
          ibm_sierraretiny5_fgoals_g3_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_fgoals_g3_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierraretiny5_fgoals_g3_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierraretiny5_fgoals_g3_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_gfdl_esm4_sensitivity_growth[[1]]$log_lambda,
          ibm_sierraretiny5_gfdl_esm4_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_gfdl_esm4_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierraretiny5_gfdl_esm4_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierraretiny5_gfdl_esm4_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_giss_e2_1_g_sensitivity_growth[[1]]$log_lambda,
          ibm_sierraretiny5_giss_e2_1_g_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_giss_e2_1_g_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierraretiny5_giss_e2_1_g_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierraretiny5_giss_e2_1_g_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_inm_cm4_8_sensitivity_growth[[1]]$log_lambda,
          ibm_sierraretiny5_inm_cm4_8_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_inm_cm4_8_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierraretiny5_inm_cm4_8_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierraretiny5_inm_cm4_8_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_growth[[1]]$log_lambda,
          ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_miroc6_sensitivity_growth[[1]]$log_lambda,
          ibm_sierraretiny5_miroc6_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_miroc6_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierraretiny5_miroc6_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierraretiny5_miroc6_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_growth[[1]]$log_lambda,
          ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_mri_esm2_0_sensitivity_growth[[1]]$log_lambda,
          ibm_sierraretiny5_mri_esm2_0_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_mri_esm2_0_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierraretiny5_mri_esm2_0_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierraretiny5_mri_esm2_0_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_noresm2_mm_sensitivity_growth[[1]]$log_lambda,
          ibm_sierraretiny5_noresm2_mm_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_noresm2_mm_sensitivity_growth[[3]]$log_lambda,
          # ibm_sierraretiny5_noresm2_mm_sensitivity_growth[[4]]$log_lambda,
          # ibm_sierraretiny5_noresm2_mm_sensitivity_growth[[5]]$log_lambda
  )))

rm(list = grep("growth", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))][grep("flowering", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))][grep("flowering", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "SierraRetinY5" &
                               sensitivityanalysis_results$vital_rate == "flowering")] = 
  c(t(rbind(ibm_sierraretiny5_canesm5_sensitivity_flowering[[1]]$log_lambda,
            ibm_sierraretiny5_canesm5_sensitivity_flowering[[2]]$log_lambda
            # ,
            # ibm_sierraretiny5_canesm5_sensitivity_flowering[[3]]$log_lambda,
            # ibm_sierraretiny5_canesm5_sensitivity_flowering[[4]]$log_lambda,
            # ibm_sierraretiny5_canesm5_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_ec_earth3_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierraretiny5_ec_earth3_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_ec_earth3_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierraretiny5_ec_earth3_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierraretiny5_ec_earth3_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_fgoals_g3_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierraretiny5_fgoals_g3_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_fgoals_g3_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierraretiny5_fgoals_g3_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierraretiny5_fgoals_g3_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_gfdl_esm4_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierraretiny5_gfdl_esm4_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_gfdl_esm4_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierraretiny5_gfdl_esm4_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierraretiny5_gfdl_esm4_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_giss_e2_1_g_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierraretiny5_giss_e2_1_g_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_giss_e2_1_g_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierraretiny5_giss_e2_1_g_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierraretiny5_giss_e2_1_g_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_inm_cm4_8_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierraretiny5_inm_cm4_8_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_inm_cm4_8_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierraretiny5_inm_cm4_8_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierraretiny5_inm_cm4_8_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_miroc6_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierraretiny5_miroc6_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_miroc6_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierraretiny5_miroc6_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierraretiny5_miroc6_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_mri_esm2_0_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierraretiny5_mri_esm2_0_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_mri_esm2_0_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierraretiny5_mri_esm2_0_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierraretiny5_mri_esm2_0_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_noresm2_mm_sensitivity_flowering[[1]]$log_lambda,
          ibm_sierraretiny5_noresm2_mm_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_noresm2_mm_sensitivity_flowering[[3]]$log_lambda,
          # ibm_sierraretiny5_noresm2_mm_sensitivity_flowering[[4]]$log_lambda,
          # ibm_sierraretiny5_noresm2_mm_sensitivity_flowering[[5]]$log_lambda
  )))

rm(list = grep("flowering", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))][grep("nbFlowers", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))][grep("nbFlowers", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "SierraRetinY5" &
                               sensitivityanalysis_results$vital_rate == "nbFlowers")] = 
  c(t(rbind(ibm_sierraretiny5_canesm5_sensitivity_nbFlowers[[1]]$log_lambda,
            ibm_sierraretiny5_canesm5_sensitivity_nbFlowers[[2]]$log_lambda
            # ,
            # ibm_sierraretiny5_canesm5_sensitivity_nbFlowers[[3]]$log_lambda,
            # ibm_sierraretiny5_canesm5_sensitivity_nbFlowers[[4]]$log_lambda,
            # ibm_sierraretiny5_canesm5_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_ec_earth3_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierraretiny5_ec_earth3_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_ec_earth3_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierraretiny5_ec_earth3_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierraretiny5_ec_earth3_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_fgoals_g3_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierraretiny5_fgoals_g3_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_fgoals_g3_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierraretiny5_fgoals_g3_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierraretiny5_fgoals_g3_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_gfdl_esm4_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierraretiny5_gfdl_esm4_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_gfdl_esm4_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierraretiny5_gfdl_esm4_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierraretiny5_gfdl_esm4_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_giss_e2_1_g_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierraretiny5_giss_e2_1_g_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_giss_e2_1_g_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierraretiny5_giss_e2_1_g_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierraretiny5_giss_e2_1_g_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_inm_cm4_8_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierraretiny5_inm_cm4_8_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_inm_cm4_8_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierraretiny5_inm_cm4_8_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierraretiny5_inm_cm4_8_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_miroc6_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierraretiny5_miroc6_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_miroc6_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierraretiny5_miroc6_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierraretiny5_miroc6_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_mri_esm2_0_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierraretiny5_mri_esm2_0_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_mri_esm2_0_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierraretiny5_mri_esm2_0_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierraretiny5_mri_esm2_0_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_noresm2_mm_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_sierraretiny5_noresm2_mm_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_noresm2_mm_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_sierraretiny5_noresm2_mm_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_sierraretiny5_noresm2_mm_sensitivity_nbFlowers[[5]]$log_lambda
  )))

rm(list = grep("nbFlowers", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))][grep("seedlingSize", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))][grep("seedlingSize", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_SierraRetinY5", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "SierraRetinY5" &
                               sensitivityanalysis_results$vital_rate == "seedlingSize")] = 
  c(t(rbind(ibm_sierraretiny5_canesm5_sensitivity_seedlingSize[[1]]$log_lambda,
            ibm_sierraretiny5_canesm5_sensitivity_seedlingSize[[2]]$log_lambda
            # ,
            # ibm_sierraretiny5_canesm5_sensitivity_seedlingSize[[3]]$log_lambda,
            # ibm_sierraretiny5_canesm5_sensitivity_seedlingSize[[4]]$log_lambda,
            # ibm_sierraretiny5_canesm5_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_ec_earth3_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierraretiny5_ec_earth3_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_ec_earth3_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierraretiny5_ec_earth3_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierraretiny5_ec_earth3_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_fgoals_g3_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierraretiny5_fgoals_g3_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_fgoals_g3_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierraretiny5_fgoals_g3_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierraretiny5_fgoals_g3_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_gfdl_esm4_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierraretiny5_gfdl_esm4_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_gfdl_esm4_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierraretiny5_gfdl_esm4_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierraretiny5_gfdl_esm4_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_giss_e2_1_g_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierraretiny5_giss_e2_1_g_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_giss_e2_1_g_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierraretiny5_giss_e2_1_g_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierraretiny5_giss_e2_1_g_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_inm_cm4_8_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierraretiny5_inm_cm4_8_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_inm_cm4_8_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierraretiny5_inm_cm4_8_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierraretiny5_inm_cm4_8_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierraretiny5_ipsl_cm6a_lr_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_miroc6_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierraretiny5_miroc6_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_miroc6_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierraretiny5_miroc6_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierraretiny5_miroc6_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierraretiny5_mpi_esm1_2_lr_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_mri_esm2_0_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierraretiny5_mri_esm2_0_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_mri_esm2_0_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierraretiny5_mri_esm2_0_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierraretiny5_mri_esm2_0_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_sierraretiny5_noresm2_mm_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_sierraretiny5_noresm2_mm_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_sierraretiny5_noresm2_mm_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_sierraretiny5_noresm2_mm_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_sierraretiny5_noresm2_mm_sensitivity_seedlingSize[[5]]$log_lambda
  )))

rm(list = grep("seedlingSize", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))][grep("survival", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))][grep("survival", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Vertedero" &
                               sensitivityanalysis_results$vital_rate == "survival")] = 
  c(t(rbind(ibm_vertedero_canesm5_sensitivity_survival[[1]]$log_lambda,
            ibm_vertedero_canesm5_sensitivity_survival[[2]]$log_lambda
            # ,
            # ibm_vertedero_canesm5_sensitivity_survival[[3]]$log_lambda,
            # ibm_vertedero_canesm5_sensitivity_survival[[4]]$log_lambda,
            # ibm_vertedero_canesm5_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_ec_earth3_sensitivity_survival[[1]]$log_lambda,
          ibm_vertedero_ec_earth3_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_vertedero_ec_earth3_sensitivity_survival[[3]]$log_lambda,
          # ibm_vertedero_ec_earth3_sensitivity_survival[[4]]$log_lambda,
          # ibm_vertedero_ec_earth3_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_fgoals_g3_sensitivity_survival[[1]]$log_lambda,
          ibm_vertedero_fgoals_g3_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_vertedero_fgoals_g3_sensitivity_survival[[3]]$log_lambda,
          # ibm_vertedero_fgoals_g3_sensitivity_survival[[4]]$log_lambda,
          # ibm_vertedero_fgoals_g3_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_gfdl_esm4_sensitivity_survival[[1]]$log_lambda,
          ibm_vertedero_gfdl_esm4_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_vertedero_gfdl_esm4_sensitivity_survival[[3]]$log_lambda,
          # ibm_vertedero_gfdl_esm4_sensitivity_survival[[4]]$log_lambda,
          # ibm_vertedero_gfdl_esm4_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_giss_e2_1_g_sensitivity_survival[[1]]$log_lambda,
          ibm_vertedero_giss_e2_1_g_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_vertedero_giss_e2_1_g_sensitivity_survival[[3]]$log_lambda,
          # ibm_vertedero_giss_e2_1_g_sensitivity_survival[[4]]$log_lambda,
          # ibm_vertedero_giss_e2_1_g_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_inm_cm4_8_sensitivity_survival[[1]]$log_lambda,
          ibm_vertedero_inm_cm4_8_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_vertedero_inm_cm4_8_sensitivity_survival[[3]]$log_lambda,
          # ibm_vertedero_inm_cm4_8_sensitivity_survival[[4]]$log_lambda,
          # ibm_vertedero_inm_cm4_8_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_ipsl_cm6a_lr_sensitivity_survival[[1]]$log_lambda,
          ibm_vertedero_ipsl_cm6a_lr_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_vertedero_ipsl_cm6a_lr_sensitivity_survival[[3]]$log_lambda,
          # ibm_vertedero_ipsl_cm6a_lr_sensitivity_survival[[4]]$log_lambda,
          # ibm_vertedero_ipsl_cm6a_lr_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_miroc6_sensitivity_survival[[1]]$log_lambda,
          ibm_vertedero_miroc6_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_vertedero_miroc6_sensitivity_survival[[3]]$log_lambda,
          # ibm_vertedero_miroc6_sensitivity_survival[[4]]$log_lambda,
          # ibm_vertedero_miroc6_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_mpi_esm1_2_lr_sensitivity_survival[[1]]$log_lambda,
          ibm_vertedero_mpi_esm1_2_lr_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_vertedero_mpi_esm1_2_lr_sensitivity_survival[[3]]$log_lambda,
          # ibm_vertedero_mpi_esm1_2_lr_sensitivity_survival[[4]]$log_lambda,
          # ibm_vertedero_mpi_esm1_2_lr_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_mri_esm2_0_sensitivity_survival[[1]]$log_lambda,
          ibm_vertedero_mri_esm2_0_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_vertedero_mri_esm2_0_sensitivity_survival[[3]]$log_lambda,
          # ibm_vertedero_mri_esm2_0_sensitivity_survival[[4]]$log_lambda,
          # ibm_vertedero_mri_esm2_0_sensitivity_survival[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_noresm2_mm_sensitivity_survival[[1]]$log_lambda,
          ibm_vertedero_noresm2_mm_sensitivity_survival[[2]]$log_lambda
          # ,
          # ibm_vertedero_noresm2_mm_sensitivity_survival[[3]]$log_lambda,
          # ibm_vertedero_noresm2_mm_sensitivity_survival[[4]]$log_lambda,
          # ibm_vertedero_noresm2_mm_sensitivity_survival[[5]]$log_lambda
  )))

rm(list = grep("survival", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))][grep("growth", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))][grep("growth", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Vertedero" &
                               sensitivityanalysis_results$vital_rate == "growth")] = 
  c(t(rbind(ibm_vertedero_canesm5_sensitivity_growth[[1]]$log_lambda,
            ibm_vertedero_canesm5_sensitivity_growth[[2]]$log_lambda
            # ,
            # ibm_vertedero_canesm5_sensitivity_growth[[3]]$log_lambda,
            # ibm_vertedero_canesm5_sensitivity_growth[[4]]$log_lambda,
            # ibm_vertedero_canesm5_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_ec_earth3_sensitivity_growth[[1]]$log_lambda,
          ibm_vertedero_ec_earth3_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_vertedero_ec_earth3_sensitivity_growth[[3]]$log_lambda,
          # ibm_vertedero_ec_earth3_sensitivity_growth[[4]]$log_lambda,
          # ibm_vertedero_ec_earth3_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_fgoals_g3_sensitivity_growth[[1]]$log_lambda,
          ibm_vertedero_fgoals_g3_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_vertedero_fgoals_g3_sensitivity_growth[[3]]$log_lambda,
          # ibm_vertedero_fgoals_g3_sensitivity_growth[[4]]$log_lambda,
          # ibm_vertedero_fgoals_g3_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_gfdl_esm4_sensitivity_growth[[1]]$log_lambda,
          ibm_vertedero_gfdl_esm4_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_vertedero_gfdl_esm4_sensitivity_growth[[3]]$log_lambda,
          # ibm_vertedero_gfdl_esm4_sensitivity_growth[[4]]$log_lambda,
          # ibm_vertedero_gfdl_esm4_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_giss_e2_1_g_sensitivity_growth[[1]]$log_lambda,
          ibm_vertedero_giss_e2_1_g_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_vertedero_giss_e2_1_g_sensitivity_growth[[3]]$log_lambda,
          # ibm_vertedero_giss_e2_1_g_sensitivity_growth[[4]]$log_lambda,
          # ibm_vertedero_giss_e2_1_g_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_inm_cm4_8_sensitivity_growth[[1]]$log_lambda,
          ibm_vertedero_inm_cm4_8_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_vertedero_inm_cm4_8_sensitivity_growth[[3]]$log_lambda,
          # ibm_vertedero_inm_cm4_8_sensitivity_growth[[4]]$log_lambda,
          # ibm_vertedero_inm_cm4_8_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_ipsl_cm6a_lr_sensitivity_growth[[1]]$log_lambda,
          ibm_vertedero_ipsl_cm6a_lr_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_vertedero_ipsl_cm6a_lr_sensitivity_growth[[3]]$log_lambda,
          # ibm_vertedero_ipsl_cm6a_lr_sensitivity_growth[[4]]$log_lambda,
          # ibm_vertedero_ipsl_cm6a_lr_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_miroc6_sensitivity_growth[[1]]$log_lambda,
          ibm_vertedero_miroc6_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_vertedero_miroc6_sensitivity_growth[[3]]$log_lambda,
          # ibm_vertedero_miroc6_sensitivity_growth[[4]]$log_lambda,
          # ibm_vertedero_miroc6_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_mpi_esm1_2_lr_sensitivity_growth[[1]]$log_lambda,
          ibm_vertedero_mpi_esm1_2_lr_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_vertedero_mpi_esm1_2_lr_sensitivity_growth[[3]]$log_lambda,
          # ibm_vertedero_mpi_esm1_2_lr_sensitivity_growth[[4]]$log_lambda,
          # ibm_vertedero_mpi_esm1_2_lr_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_mri_esm2_0_sensitivity_growth[[1]]$log_lambda,
          ibm_vertedero_mri_esm2_0_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_vertedero_mri_esm2_0_sensitivity_growth[[3]]$log_lambda,
          # ibm_vertedero_mri_esm2_0_sensitivity_growth[[4]]$log_lambda,
          # ibm_vertedero_mri_esm2_0_sensitivity_growth[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_noresm2_mm_sensitivity_growth[[1]]$log_lambda,
          ibm_vertedero_noresm2_mm_sensitivity_growth[[2]]$log_lambda
          # ,
          # ibm_vertedero_noresm2_mm_sensitivity_growth[[3]]$log_lambda,
          # ibm_vertedero_noresm2_mm_sensitivity_growth[[4]]$log_lambda,
          # ibm_vertedero_noresm2_mm_sensitivity_growth[[5]]$log_lambda
  )))

rm(list = grep("growth", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))][grep("flowering", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))][grep("flowering", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Vertedero" &
                               sensitivityanalysis_results$vital_rate == "flowering")] = 
  c(t(rbind(ibm_vertedero_canesm5_sensitivity_flowering[[1]]$log_lambda,
            ibm_vertedero_canesm5_sensitivity_flowering[[2]]$log_lambda
            # ,
            # ibm_vertedero_canesm5_sensitivity_flowering[[3]]$log_lambda,
            # ibm_vertedero_canesm5_sensitivity_flowering[[4]]$log_lambda,
            # ibm_vertedero_canesm5_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_ec_earth3_sensitivity_flowering[[1]]$log_lambda,
          ibm_vertedero_ec_earth3_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_vertedero_ec_earth3_sensitivity_flowering[[3]]$log_lambda,
          # ibm_vertedero_ec_earth3_sensitivity_flowering[[4]]$log_lambda,
          # ibm_vertedero_ec_earth3_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_fgoals_g3_sensitivity_flowering[[1]]$log_lambda,
          ibm_vertedero_fgoals_g3_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_vertedero_fgoals_g3_sensitivity_flowering[[3]]$log_lambda,
          # ibm_vertedero_fgoals_g3_sensitivity_flowering[[4]]$log_lambda,
          # ibm_vertedero_fgoals_g3_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_gfdl_esm4_sensitivity_flowering[[1]]$log_lambda,
          ibm_vertedero_gfdl_esm4_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_vertedero_gfdl_esm4_sensitivity_flowering[[3]]$log_lambda,
          # ibm_vertedero_gfdl_esm4_sensitivity_flowering[[4]]$log_lambda,
          # ibm_vertedero_gfdl_esm4_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_giss_e2_1_g_sensitivity_flowering[[1]]$log_lambda,
          ibm_vertedero_giss_e2_1_g_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_vertedero_giss_e2_1_g_sensitivity_flowering[[3]]$log_lambda,
          # ibm_vertedero_giss_e2_1_g_sensitivity_flowering[[4]]$log_lambda,
          # ibm_vertedero_giss_e2_1_g_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_inm_cm4_8_sensitivity_flowering[[1]]$log_lambda,
          ibm_vertedero_inm_cm4_8_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_vertedero_inm_cm4_8_sensitivity_flowering[[3]]$log_lambda,
          # ibm_vertedero_inm_cm4_8_sensitivity_flowering[[4]]$log_lambda,
          # ibm_vertedero_inm_cm4_8_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_ipsl_cm6a_lr_sensitivity_flowering[[1]]$log_lambda,
          ibm_vertedero_ipsl_cm6a_lr_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_vertedero_ipsl_cm6a_lr_sensitivity_flowering[[3]]$log_lambda,
          # ibm_vertedero_ipsl_cm6a_lr_sensitivity_flowering[[4]]$log_lambda,
          # ibm_vertedero_ipsl_cm6a_lr_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_miroc6_sensitivity_flowering[[1]]$log_lambda,
          ibm_vertedero_miroc6_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_vertedero_miroc6_sensitivity_flowering[[3]]$log_lambda,
          # ibm_vertedero_miroc6_sensitivity_flowering[[4]]$log_lambda,
          # ibm_vertedero_miroc6_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_mpi_esm1_2_lr_sensitivity_flowering[[1]]$log_lambda,
          ibm_vertedero_mpi_esm1_2_lr_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_vertedero_mpi_esm1_2_lr_sensitivity_flowering[[3]]$log_lambda,
          # ibm_vertedero_mpi_esm1_2_lr_sensitivity_flowering[[4]]$log_lambda,
          # ibm_vertedero_mpi_esm1_2_lr_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_mri_esm2_0_sensitivity_flowering[[1]]$log_lambda,
          ibm_vertedero_mri_esm2_0_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_vertedero_mri_esm2_0_sensitivity_flowering[[3]]$log_lambda,
          # ibm_vertedero_mri_esm2_0_sensitivity_flowering[[4]]$log_lambda,
          # ibm_vertedero_mri_esm2_0_sensitivity_flowering[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_noresm2_mm_sensitivity_flowering[[1]]$log_lambda,
          ibm_vertedero_noresm2_mm_sensitivity_flowering[[2]]$log_lambda
          # ,
          # ibm_vertedero_noresm2_mm_sensitivity_flowering[[3]]$log_lambda,
          # ibm_vertedero_noresm2_mm_sensitivity_flowering[[4]]$log_lambda,
          # ibm_vertedero_noresm2_mm_sensitivity_flowering[[5]]$log_lambda
  )))

rm(list = grep("flowering", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))][grep("nbFlowers", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))][grep("nbFlowers", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Vertedero" &
                               sensitivityanalysis_results$vital_rate == "nbFlowers")] = 
  c(t(rbind(ibm_vertedero_canesm5_sensitivity_nbFlowers[[1]]$log_lambda,
            ibm_vertedero_canesm5_sensitivity_nbFlowers[[2]]$log_lambda
            # ,
            # ibm_vertedero_canesm5_sensitivity_nbFlowers[[3]]$log_lambda,
            # ibm_vertedero_canesm5_sensitivity_nbFlowers[[4]]$log_lambda,
            # ibm_vertedero_canesm5_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_ec_earth3_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_vertedero_ec_earth3_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_vertedero_ec_earth3_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_vertedero_ec_earth3_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_vertedero_ec_earth3_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_fgoals_g3_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_vertedero_fgoals_g3_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_vertedero_fgoals_g3_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_vertedero_fgoals_g3_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_vertedero_fgoals_g3_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_gfdl_esm4_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_vertedero_gfdl_esm4_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_vertedero_gfdl_esm4_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_vertedero_gfdl_esm4_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_vertedero_gfdl_esm4_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_giss_e2_1_g_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_vertedero_giss_e2_1_g_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_vertedero_giss_e2_1_g_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_vertedero_giss_e2_1_g_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_vertedero_giss_e2_1_g_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_inm_cm4_8_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_vertedero_inm_cm4_8_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_vertedero_inm_cm4_8_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_vertedero_inm_cm4_8_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_vertedero_inm_cm4_8_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_ipsl_cm6a_lr_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_vertedero_ipsl_cm6a_lr_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_vertedero_ipsl_cm6a_lr_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_vertedero_ipsl_cm6a_lr_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_vertedero_ipsl_cm6a_lr_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_miroc6_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_vertedero_miroc6_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_vertedero_miroc6_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_vertedero_miroc6_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_vertedero_miroc6_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_mpi_esm1_2_lr_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_vertedero_mpi_esm1_2_lr_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_vertedero_mpi_esm1_2_lr_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_vertedero_mpi_esm1_2_lr_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_vertedero_mpi_esm1_2_lr_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_mri_esm2_0_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_vertedero_mri_esm2_0_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_vertedero_mri_esm2_0_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_vertedero_mri_esm2_0_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_vertedero_mri_esm2_0_sensitivity_nbFlowers[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_noresm2_mm_sensitivity_nbFlowers[[1]]$log_lambda,
          ibm_vertedero_noresm2_mm_sensitivity_nbFlowers[[2]]$log_lambda
          # ,
          # ibm_vertedero_noresm2_mm_sensitivity_nbFlowers[[3]]$log_lambda,
          # ibm_vertedero_noresm2_mm_sensitivity_nbFlowers[[4]]$log_lambda,
          # ibm_vertedero_noresm2_mm_sensitivity_nbFlowers[[5]]$log_lambda
  )))

rm(list = grep("nbFlowers", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))][grep("seedlingSize", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))])])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))][grep("seedlingSize", list.files("Output/Projections/")[grep("Natural_SensitivityAnalysis_Vertedero", list.files("Output/Projections/"))])][i]))
}

# Add log lambda (per simulation = per row)

sensitivityanalysis_results$log_lambda[which(sensitivityanalysis_results$population == "Vertedero" &
                               sensitivityanalysis_results$vital_rate == "seedlingSize")] = 
  c(t(rbind(ibm_vertedero_canesm5_sensitivity_seedlingSize[[1]]$log_lambda,
            ibm_vertedero_canesm5_sensitivity_seedlingSize[[2]]$log_lambda
            # ,
            # ibm_vertedero_canesm5_sensitivity_seedlingSize[[3]]$log_lambda,
            # ibm_vertedero_canesm5_sensitivity_seedlingSize[[4]]$log_lambda,
            # ibm_vertedero_canesm5_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_ec_earth3_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_vertedero_ec_earth3_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_vertedero_ec_earth3_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_vertedero_ec_earth3_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_vertedero_ec_earth3_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_fgoals_g3_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_vertedero_fgoals_g3_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_vertedero_fgoals_g3_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_vertedero_fgoals_g3_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_vertedero_fgoals_g3_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_gfdl_esm4_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_vertedero_gfdl_esm4_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_vertedero_gfdl_esm4_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_vertedero_gfdl_esm4_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_vertedero_gfdl_esm4_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_giss_e2_1_g_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_vertedero_giss_e2_1_g_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_vertedero_giss_e2_1_g_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_vertedero_giss_e2_1_g_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_vertedero_giss_e2_1_g_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_inm_cm4_8_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_vertedero_inm_cm4_8_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_vertedero_inm_cm4_8_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_vertedero_inm_cm4_8_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_vertedero_inm_cm4_8_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_ipsl_cm6a_lr_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_vertedero_ipsl_cm6a_lr_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_vertedero_ipsl_cm6a_lr_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_vertedero_ipsl_cm6a_lr_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_vertedero_ipsl_cm6a_lr_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_miroc6_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_vertedero_miroc6_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_vertedero_miroc6_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_vertedero_miroc6_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_vertedero_miroc6_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_mpi_esm1_2_lr_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_vertedero_mpi_esm1_2_lr_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_vertedero_mpi_esm1_2_lr_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_vertedero_mpi_esm1_2_lr_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_vertedero_mpi_esm1_2_lr_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_mri_esm2_0_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_vertedero_mri_esm2_0_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_vertedero_mri_esm2_0_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_vertedero_mri_esm2_0_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_vertedero_mri_esm2_0_sensitivity_seedlingSize[[5]]$log_lambda
  )),
  t(rbind(ibm_vertedero_noresm2_mm_sensitivity_seedlingSize[[1]]$log_lambda,
          ibm_vertedero_noresm2_mm_sensitivity_seedlingSize[[2]]$log_lambda
          # ,
          # ibm_vertedero_noresm2_mm_sensitivity_seedlingSize[[3]]$log_lambda,
          # ibm_vertedero_noresm2_mm_sensitivity_seedlingSize[[4]]$log_lambda,
          # ibm_vertedero_noresm2_mm_sensitivity_seedlingSize[[5]]$log_lambda
  )))

rm(list = grep("seedlingSize", ls() , value = TRUE, invert = FALSE))

write.csv(sensitivityanalysis_results, file = "Output/Projections/SensitivityAnalysis_Results.csv", row.names = F)


## 2.2. Summarising results ----
# -------------------------

sensitivityanalysis_results = read.csv("Output/Projections/SensitivityAnalysis_Results.csv")

sensitivityanalysis_results$habitat = factor(sensitivityanalysis_results$habitat,
                                             levels = c("Natural", "Anthropogenic"))

mean_lambda_simulation_control = 
  aggregate(log_lambda ~ habitat + population + simulation, 
            data = ibm_results[which(!is.infinite(ibm_results$log_lambda)), ],
            FUN = mean)

mean_lambda_simulation_control = rbind(mean_lambda_simulation_control[which(mean_lambda_simulation_control$habitat == "Natural"), ],
                                       mean_lambda_simulation_control[which(mean_lambda_simulation_control$habitat == "Anthropogenic"), ])

mean_lambda_simulation_control$population_full = mean_lambda_simulation_control$population

mean_lambda_simulation_control$population_full[which(mean_lambda_simulation_control$population == "SierraCarboneraY5")] = "Sierra\nCarbonera\nYoung"
mean_lambda_simulation_control$population_full[which(mean_lambda_simulation_control$population == "SierraRetinY5")] = "Sierra del\nRetín\nYoung"
mean_lambda_simulation_control$population_full[which(mean_lambda_simulation_control$population == "MonteraTorero")] = "Montera\ndel\nTorero"
mean_lambda_simulation_control$population_full[which(mean_lambda_simulation_control$population == "Retin")] = "Sierra del\nRetín\nDisturbed"
mean_lambda_simulation_control$population_full[which(mean_lambda_simulation_control$population == "Prisoneros")] = "Prisioneros"
mean_lambda_simulation_control$population_full[which(mean_lambda_simulation_control$population == "SCarbDist")] = "Sierra\nCarbonera\nDisturbed"



mean_lambda_simulation = 
  aggregate(log_lambda ~ habitat + climate_model + vital_rate + population + simulation, 
            data = sensitivityanalysis_results[which(!is.infinite(sensitivityanalysis_results$log_lambda)), ],
            FUN = mean)

mean_lambda_simulation = rbind(mean_lambda_simulation[which(mean_lambda_simulation$habitat == "Natural"), ],
                               mean_lambda_simulation[which(mean_lambda_simulation$habitat == "Anthropogenic"), ])

mean_lambda_simulation$population_full = mean_lambda_simulation$population

mean_lambda_simulation$population_full[which(mean_lambda_simulation$population == "SierraCarboneraY5")] = "Sierra\nCarbonera\nYoung"
mean_lambda_simulation$population_full[which(mean_lambda_simulation$population == "SierraRetinY5")] = "Sierra del\nRetín\nYoung"
mean_lambda_simulation$population_full[which(mean_lambda_simulation$population == "MonteraTorero")] = "Montera\ndel\nTorero"
mean_lambda_simulation$population_full[which(mean_lambda_simulation$population == "Retin")] = "Sierra del\nRetín\nDisturbed"
mean_lambda_simulation$population_full[which(mean_lambda_simulation$population == "Prisoneros")] = "Prisioneros"
mean_lambda_simulation$population_full[which(mean_lambda_simulation$population == "SCarbDist")] = "Sierra\nCarbonera\nDisturbed"



mean_lambda_control = 
  aggregate(log_lambda ~ habitat + population, 
            data = ibm_results[which(!is.infinite(ibm_results$log_lambda)), ],
            FUN = mean)

mean_lambda_control = rbind(mean_lambda_control[which(mean_lambda_control$habitat == "Natural"), ],
                            mean_lambda_control[which(mean_lambda_control$habitat == "Anthropogenic"), ])

mean_lambda_control$population_full = mean_lambda_control$population

mean_lambda_control$population_full[which(mean_lambda_control$population == "SierraCarboneraY5")] = "Sierra\nCarbonera\nYoung"
mean_lambda_control$population_full[which(mean_lambda_control$population == "SierraRetinY5")] = "Sierra del\nRetín\nYoung"
mean_lambda_control$population_full[which(mean_lambda_control$population == "MonteraTorero")] = "Montera\ndel\nTorero"
mean_lambda_control$population_full[which(mean_lambda_control$population == "Retin")] = "Sierra del\nRetín\nDisturbed"
mean_lambda_control$population_full[which(mean_lambda_control$population == "Prisoneros")] = "Prisioneros"
mean_lambda_control$population_full[which(mean_lambda_control$population == "SCarbDist")] = "Sierra\nCarbonera\nDisturbed"



mean_lambda = 
  aggregate(log_lambda ~ habitat + climate_model + vital_rate + population, 
            data = sensitivityanalysis_results[which(!is.infinite(sensitivityanalysis_results$log_lambda)), ],
            FUN = mean)

mean_lambda = rbind(mean_lambda[which(mean_lambda$habitat == "Natural"), ],
                    mean_lambda[which(mean_lambda$habitat == "Anthropogenic"), ])

mean_lambda$population_full = mean_lambda$population

mean_lambda$population_full[which(mean_lambda$population == "SierraCarboneraY5")] = "Sierra\nCarbonera\nYoung"
mean_lambda$population_full[which(mean_lambda$population == "SierraRetinY5")] = "Sierra del\nRetín\nYoung"
mean_lambda$population_full[which(mean_lambda$population == "MonteraTorero")] = "Montera\ndel\nTorero"
mean_lambda$population_full[which(mean_lambda$population == "Retin")] = "Sierra del\nRetín\nDisturbed"
mean_lambda$population_full[which(mean_lambda$population == "Prisoneros")] = "Prisioneros"
mean_lambda$population_full[which(mean_lambda$population == "SCarbDist")] = "Sierra\nCarbonera\nDisturbed"




###########################################################################
#
# 3. Computing sensitivity values ----
#
###########################################################################

n_iter = 500

sensitivity_results = expand.grid(iteration = seq(1, n_iter), 
                                  climate_model = climate_change_models,
                                  vital_rate = c("survival", "growth", 
                                                 "flowering", "nbFlowers",
                                                 "seedlingSize"),
                                  habitat = c(rep("Natural",
                                                  times = length(populations_natural)), 
                                              rep("Anthropogenic", 
                                                  times = length(populations_anthropogenic))))

sensitivity_results$population = c(rep(populations_natural, each = n_iter * (length(climate_change_models)) * 5),
                                   rep(populations_anthropogenic, each = n_iter * (length(climate_change_models)) * 5))

sensitivity_results = rbind(cbind(sensitivity_results, method = "Across simulations"),
                            cbind(sensitivity_results, method = "Per simulation"))


sensitivity_results$treatment = paste(sensitivity_results$population,
                                      sensitivity_results$climate_model,
                                      sensitivity_results$vital_rate,
                                      sensitivity_results$method,
                                      # sensitivity_results$simulation,
                                      sep = "_")

sensitivity_results$sensitivity = NA


## 3.1. Differences between mean across simulations ----
# ------------------------------------------------

# For each population and climate model: 
# We do XX iterations where we sample 100 control projections among the original 500
# and calculate the percentage difference between the mean stochastic log lambda 
# of these 100 projections and the mean stochastic log lambda of the projection
# with each vital rate perturbed.

for(population in unique(sensitivity_results$population)){
  
  for(climate_model in unique(sensitivity_results$climate_model)){
    
    for(vital_rate in unique(sensitivity_results$vital_rate)){
      
      print(paste(population, climate_model, vital_rate, sep = " "))
      
      for(iter in 1:n_iter){
        
        sampled_sim = sample(unique(mean_lambda_simulation_control$simulation), 100, replace = F)
        
        stoch_lambda_control = mean_lambda_simulation_control$log_lambda[which(mean_lambda_simulation_control$population == population & 
                                                                               mean_lambda_simulation_control$simulation %in% sampled_sim)]
        
        stoch_lambda_perturbed = mean_lambda_simulation$log_lambda[which(mean_lambda_simulation$population == population & 
                                                                         mean_lambda_simulation$vital_rate == vital_rate &
                                                                         mean_lambda_simulation$climate_model == climate_model)]
        
        sensitivity_value = (mean(stoch_lambda_control) - mean(stoch_lambda_perturbed))/mean(stoch_lambda_control)
        
        sensitivity_results$sensitivity[which(sensitivity_results$method == "Across simulations" &
                                              sensitivity_results$population == population & 
                                              sensitivity_results$vital_rate == vital_rate &
                                              sensitivity_results$climate_model == climate_model &
                                              sensitivity_results$iteration == iter)] = sensitivity_value
      }
    }
  }
}

write.csv(sensitivity_results, file = "Output/Projections/ProcessedSensitivityResults.csv", row.names = F)


## 3.2. Differences between mean per simulations ----
# ----------------------------------------------

# For each population and climate model: 
# We sample 100 control projections among the original 500 and calculate the 
# pairwise percentage difference between the stochastic log lambdas 
# of each 100 projections and the stochastic log lambdas of the projections
# with each vital rate perturbed. We change the pairs at each of the XX iterations.

sampled_sim = sample(unique(mean_lambda_simulation_control$simulation), 100, replace = F)

stoch_lambda_control = mean_lambda_simulation_control$log_lambda[which(mean_lambda_simulation_control$population == population & 
                                                                       mean_lambda_simulation_control$simulation %in% sampled_sim)]

for(population in unique(sensitivity_results$population)){
  
  for(climate_model in unique(sensitivity_results$climate_model)){
    
    for(vital_rate in unique(sensitivity_results$vital_rate)){
      
      print(paste(population, climate_model, vital_rate, sep = " "))
      
      for(iter in 1:n_iter){
        
        sampled_pair = sample(unique(mean_lambda_simulation$simulation), 100)
        
        stoch_lambda_perturbed = mean_lambda_simulation$log_lambda[which(mean_lambda_simulation$population == population & 
                                                                         mean_lambda_simulation$vital_rate == vital_rate &
                                                                         mean_lambda_simulation$climate_model == climate_model)][sampled_pair]
        
        sensitivity_values = (stoch_lambda_control - stoch_lambda_perturbed)/stoch_lambda_control
        
        sensitivity_results$sensitivity[which(sensitivity_results$method == "Per simulation" &
                                              sensitivity_results$population == population & 
                                              sensitivity_results$vital_rate == vital_rate &
                                              sensitivity_results$climate_model == climate_model &
                                              sensitivity_results$iteration == iter)] = mean(sensitivity_values)
      }
    }
  }
}


write.csv(sensitivity_results, file = "Output/Projections/ProcessedSensitivityResults.csv", row.names = F)


sensitivity_results = read.csv("Output/Projections/ProcessedSensitivityResults.csv")

mean_sensitivity = aggregate(sensitivity ~ habitat + vital_rate + population + method + iteration, 
                             sensitivity_results,
                             mean)

mean_sensitivity_agg = aggregate(sensitivity ~ habitat + population + method + vital_rate, 
                             sensitivity_results,
                             FUN = function(x) c(mean(x), quantile(x, probs = c(0.025, 0.975), na.rm = T)))


mean_sensitivity_habitat = aggregate(sensitivity ~ habitat + population + vital_rate + method + iteration, 
                             sensitivity_results,
                             mean)

mean_sensitivity_habitat_agg = aggregate(sensitivity ~ habitat + method + vital_rate, 
                                     sensitivity_results,
                                     FUN = function(x) c(mean(x), quantile(x, probs = c(0.025, 0.975), na.rm = T)))


mean_sensitivity$population_full = mean_sensitivity$population

mean_sensitivity$population_full[which(mean_sensitivity$population == "SierraCarboneraY5")] = "Sierra\nCarbonera\nYoung"
mean_sensitivity$population_full[which(mean_sensitivity$population == "SierraRetinY5")] = "Sierra del\nRetín\nYoung"
mean_sensitivity$population_full[which(mean_sensitivity$population == "MonteraTorero")] = "Montera\ndel\nTorero"
mean_sensitivity$population_full[which(mean_sensitivity$population == "Retin")] = "Sierra del\nRetín\nDisturbed"
mean_sensitivity$population_full[which(mean_sensitivity$population == "Prisoneros")] = "Prisioneros"
mean_sensitivity$population_full[which(mean_sensitivity$population == "SCarbDist")] = "Sierra\nCarbonera\nDisturbed"

mean_sensitivity$population_full = factor(mean_sensitivity$population_full,
                                              levels = c("Sierra\nCarbonera\nYoung",
                                                         "Sierra del\nRetín\nYoung", "Vertedero",
                                                         "Bujeo", "Montera\ndel\nTorero",
                                                         "Prisioneros", "Sierra\nCarbonera\nDisturbed",
                                                         "Sierra del\nRetín\nDisturbed"))


mean_sensitivity$vital_rate_full = mean_sensitivity$vital_rate

mean_sensitivity$vital_rate_full[which(mean_sensitivity$vital_rate == "survival")] = "Survival"
mean_sensitivity$vital_rate_full[which(mean_sensitivity$vital_rate == "growth")] = "Growth"
mean_sensitivity$vital_rate_full[which(mean_sensitivity$vital_rate == "flowering")] = "Flowering"
mean_sensitivity$vital_rate_full[which(mean_sensitivity$vital_rate == "nbFlowers")] = "Number of flowers"
mean_sensitivity$vital_rate_full[which(mean_sensitivity$vital_rate == "seedlingSize")] = "Seedling size"

mean_sensitivity$vital_rate_full = factor(mean_sensitivity$vital_rate_full,
                                          levels = c("Survival", "Growth", "Flowering",
                                                     "Number of flowers", "Seedling size"))



mean_sensitivity_agg$population_full = mean_sensitivity_agg$population

mean_sensitivity_agg$population_full[which(mean_sensitivity_agg$population == "SierraCarboneraY5")] = "Sierra\nCarbonera\nYoung"
mean_sensitivity_agg$population_full[which(mean_sensitivity_agg$population == "SierraRetinY5")] = "Sierra del\nRetín\nYoung"
mean_sensitivity_agg$population_full[which(mean_sensitivity_agg$population == "MonteraTorero")] = "Montera\ndel\nTorero"
mean_sensitivity_agg$population_full[which(mean_sensitivity_agg$population == "Retin")] = "Sierra del\nRetín\nDisturbed"
mean_sensitivity_agg$population_full[which(mean_sensitivity_agg$population == "Prisoneros")] = "Prisioneros"
mean_sensitivity_agg$population_full[which(mean_sensitivity_agg$population == "SCarbDist")] = "Sierra\nCarbonera\nDisturbed"

mean_sensitivity_agg$population_full = factor(mean_sensitivity_agg$population_full,
                                              levels = c("Sierra\nCarbonera\nYoung",
                                                         "Sierra del\nRetín\nYoung", "Vertedero",
                                                         "Bujeo", "Montera\ndel\nTorero",
                                                         "Prisioneros", "Sierra\nCarbonera\nDisturbed",
                                                         "Sierra del\nRetín\nDisturbed"))
mean_sensitivity_agg = mean_sensitivity_agg[order(mean_sensitivity_agg$population_full), ]


mean_sensitivity_agg$vital_rate_full = mean_sensitivity_agg$vital_rate

mean_sensitivity_agg$vital_rate_full[which(mean_sensitivity_agg$vital_rate == "survival")] = "Survival"
mean_sensitivity_agg$vital_rate_full[which(mean_sensitivity_agg$vital_rate == "growth")] = "Growth"
mean_sensitivity_agg$vital_rate_full[which(mean_sensitivity_agg$vital_rate == "flowering")] = "Flowering"
mean_sensitivity_agg$vital_rate_full[which(mean_sensitivity_agg$vital_rate == "nbFlowers")] = "Number of flowers"
mean_sensitivity_agg$vital_rate_full[which(mean_sensitivity_agg$vital_rate == "seedlingSize")] = "Seedling size"

mean_sensitivity_agg$vital_rate_full = factor(mean_sensitivity_agg$vital_rate_full,
                                              levels = c("Survival", "Growth", "Flowering",
                                                         "Number of flowers", "Seedling size"))
mean_sensitivity_agg = mean_sensitivity_agg[order(mean_sensitivity_agg$vital_rate_full), ]




mean_sensitivity_habitat$vital_rate_full = mean_sensitivity_habitat$vital_rate

mean_sensitivity_habitat$vital_rate_full[which(mean_sensitivity_habitat$vital_rate == "survival")] = "Survival"
mean_sensitivity_habitat$vital_rate_full[which(mean_sensitivity_habitat$vital_rate == "growth")] = "Growth"
mean_sensitivity_habitat$vital_rate_full[which(mean_sensitivity_habitat$vital_rate == "flowering")] = "Flowering"
mean_sensitivity_habitat$vital_rate_full[which(mean_sensitivity_habitat$vital_rate == "nbFlowers")] = "Number of flowers"
mean_sensitivity_habitat$vital_rate_full[which(mean_sensitivity_habitat$vital_rate == "seedlingSize")] = "Seedling size"

mean_sensitivity_habitat$vital_rate_full = factor(mean_sensitivity_habitat$vital_rate_full,
                                                  levels = c("Survival", "Growth", "Flowering",
                                                             "Number of flowers", "Seedling size"))


mean_sensitivity_habitat_agg$vital_rate_full = mean_sensitivity_habitat_agg$vital_rate

mean_sensitivity_habitat_agg$vital_rate_full[which(mean_sensitivity_habitat_agg$vital_rate == "survival")] = "Survival"
mean_sensitivity_habitat_agg$vital_rate_full[which(mean_sensitivity_habitat_agg$vital_rate == "growth")] = "Growth"
mean_sensitivity_habitat_agg$vital_rate_full[which(mean_sensitivity_habitat_agg$vital_rate == "flowering")] = "Flowering"
mean_sensitivity_habitat_agg$vital_rate_full[which(mean_sensitivity_habitat_agg$vital_rate == "nbFlowers")] = "Number of flowers"
mean_sensitivity_habitat_agg$vital_rate_full[which(mean_sensitivity_habitat_agg$vital_rate == "seedlingSize")] = "Seedling size"

mean_sensitivity_habitat_agg$vital_rate_full = factor(mean_sensitivity_habitat_agg$vital_rate_full,
                                                      levels = c("Survival", "Growth", "Flowering",
                                                                 "Number of flowers", "Seedling size"))
mean_sensitivity_habitat_agg = mean_sensitivity_habitat_agg[order(mean_sensitivity_habitat_agg$vital_rate_full), ]




###########################################################################
#
# 4. Plotting results ----
#
###########################################################################

cbbPalette <- c("#111111", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## 4.1. Sensitivity per population ----
# --------------------------------

mean_sensitivity$habitat = factor(mean_sensitivity$habitat,
                                  levels = c("Natural", "Anthropogenic"))

png(filename = "Output/Plots/Sensitivity_AcrossSimulations_PerPopulation.png", 
    width = 12,
    height = 20,
    units = "cm",
    bg = "white",
    res = 600)

ggplot(mean_sensitivity[which(mean_sensitivity$method == "Across simulations"), ], 
       aes(x = population_full, y = sensitivity, 
           fill = habitat, colour = habitat)) +
  facet_wrap(~ vital_rate_full, scales = "free", ncol = 1) +
  geom_boxplot(outlier.size = 0.2, size = 0.2, alpha = 0.5) +
  stat_summary(fun = mean, colour = cbbPalette[5], aes(shape = "Mean", group = habitat), geom = "point", shape = 17, size = 1, position = position_dodge(width = 0.75)) +
  labs(x = "Population", 
       y = "Sensitivity") +
  scale_fill_manual(values = cbbPalette[1:2],
                    name = "Habitat type") +
  scale_colour_manual(values = cbbPalette[1:2],
                      name = "Habitat type") +
  scale_shape_manual("", values = c("Mean" = 17)) +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title.x = element_text(size = 8, colour = "black", 
                                    margin = margin(t = 0, r = 0, b = 10, l = 0)), 
        axis.title.y = element_text(size = 8, colour = "black", angle = 90, 
                                    margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 2, r = 0, b = 5, l = 0)), 
        axis.text.y = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        panel.grid = element_blank(),
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 8), 
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(1, "lines"),
        strip.background = element_blank())

dev.off()



png(filename = "Output/Plots/Sensitivity_PerSimulation_PerPopulation.png", 
    width = 12,
    height = 20,
    units = "cm",
    bg = "white",
    res = 600)

ggplot(mean_sensitivity[which(mean_sensitivity$method == "Per simulation"), ], 
       aes(x = population_full, y = sensitivity, 
           fill = habitat, colour = habitat)) +
  facet_wrap(~ vital_rate_full, scales = "free", ncol = 1) +
  geom_boxplot(outlier.size = 0.2, size = 0.2, alpha = 0.5) +
  stat_summary(fun = mean, colour = cbbPalette[5], aes(shape = "Mean", group = habitat), geom = "point", shape = 17, size = 1, position = position_dodge(width = 0.75)) +
  labs(x = "Population", 
       y = "Sensitivity") +
  scale_fill_manual(values = cbbPalette[1:2],
                    name = "Habitat type") +
  scale_colour_manual(values = cbbPalette[1:2],
                      name = "Habitat type") +
  scale_shape_manual("", values = c("Mean" = 17)) +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title.x = element_text(size = 8, colour = "black", 
                                    margin = margin(t = 0, r = 0, b = 10, l = 0)), 
        axis.title.y = element_text(size = 8, colour = "black", angle = 90, 
                                    margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 2, r = 0, b = 5, l = 0)), 
        axis.text.y = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        panel.grid = element_blank(),
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 8), 
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(1, "lines"),
        strip.background = element_blank())

dev.off()


## 4.2. Sensitivity per habitat ----
# -----------------------------

mean_sensitivity_habitat$habitat = factor(mean_sensitivity_habitat$habitat,
                                          levels = c("Natural", "Anthropogenic"))

png(filename = "Output/Plots/Sensitivity_AcrossSimulations_PerHabitat.png", 
    width = 12,
    height = 10,
    units = "cm",
    bg = "white",
    res = 600)

ggplot(mean_sensitivity_habitat[which(mean_sensitivity_habitat$method == "Across simulations"), ], 
       aes(x = vital_rate_full, y = sensitivity, 
           fill = habitat, colour = habitat)) +
  geom_boxplot(outlier.size = 0.2, size = 0.2, alpha = 0.5) +
  stat_summary(fun = mean, colour = cbbPalette[5], aes(shape = "Mean", group = habitat), geom = "point", shape = 17, size = 1, position = position_dodge(width = 0.75)) +
  labs(x = "Vital rate", 
       y = "Sensitivity") +
  scale_fill_manual(values = cbbPalette[1:2],
                    name = "Habitat type") +
  scale_colour_manual(values = cbbPalette[1:2],
                      name = "Habitat type") +
  scale_shape_manual("", values = c("Mean" = 17)) +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title.x = element_text(size = 8, colour = "black", 
                                    margin = margin(t = 0, r = 0, b = 10, l = 0)), 
        axis.title.y = element_text(size = 8, colour = "black", angle = 90, 
                                    margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 2, r = 0, b = 5, l = 0)), 
        axis.text.y = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        panel.grid = element_blank(),
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 8), 
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(1, "lines"),
        strip.background = element_blank())

dev.off()



png(filename = "Output/Plots/Sensitivity_PerSimulation_PerHabitat.png", 
    width = 12,
    height = 10,
    units = "cm",
    bg = "white",
    res = 600)

ggplot(mean_sensitivity_habitat[which(mean_sensitivity_habitat$method == "Per simulation"), ], 
       aes(x = vital_rate_full, y = sensitivity, 
           fill = habitat, colour = habitat)) +
  geom_boxplot(outlier.size = 0.2, size = 0.2, alpha = 0.5) +
  stat_summary(fun = mean, colour = cbbPalette[5], aes(shape = "Mean", group = habitat), geom = "point", shape = 17, size = 1, position = position_dodge(width = 0.75)) +
  labs(x = "Population", 
       y = "Sensitivity") +
  scale_fill_manual(values = cbbPalette[1:2],
                    name = "Habitat type") +
  scale_colour_manual(values = cbbPalette[1:2],
                      name = "Habitat type") +
  scale_shape_manual("", values = c("Mean" = 17)) +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title.x = element_text(size = 8, colour = "black", 
                                    margin = margin(t = 0, r = 0, b = 10, l = 0)), 
        axis.title.y = element_text(size = 8, colour = "black", angle = 90, 
                                    margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 2, r = 0, b = 5, l = 0)), 
        axis.text.y = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        panel.grid = element_blank(),
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 8), 
        legend.position = "bottom", 
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(1, "lines"),
        strip.background = element_blank())

dev.off()