############################################################################
#
# This script processes the results of the projections of the dewy-pine
# IBMs for natural and anthropogenicpopulations under current and future
# climatic conditions.
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




###########################################################################
#
# 2. Processing results ----
#
###########################################################################

n_sim = 500
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

ibm_results = expand.grid(year = seq(1, n_years),
                          simulation = seq(1, n_sim), 
                          climate = c("Control", 
                                      rep("Climate change", 
                                          times = length(climate_change_models))),
                          habitat = c(rep("Natural",
                                           times = length(populations_natural)), 
                                       rep("Anthropogenic", 
                                           times = length(populations_anthropogenic))))


# Add populations

ibm_results$population = c(rep(populations_natural, each = n_years * n_sim * (length(climate_change_models) + 1)),
                           rep(populations_anthropogenic, each = n_years * n_sim * (length(climate_change_models) + 1)))


# Add climate models

ibm_results$climate_model = NA
ibm_results$climate_model[which(ibm_results$climate == "Control")] = "Control"
ibm_results$climate_model[which(ibm_results$climate == "Climate change")] =
  rep(climate_change_models, each = n_years * n_sim)

ibm_results$treatment = paste(ibm_results$population,
                              ibm_results$climate_model,
                              ibm_results$simulation,
                              sep = "_")


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_Retin", list.files("Output/Projections/"))])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_Retin", list.files("Output/Projections/"))][i]))
}


# Add log lambda (per simulation = per row)

ibm_results$log_lambda = NA

ibm_results$log_lambda[which(ibm_results$population == "Retin")] = 
  c(t(rbind(ibm_retin_control[[1]]$log_lambda,
            ibm_retin_control[[2]]$log_lambda,
            ibm_retin_control[[3]]$log_lambda,
            ibm_retin_control[[4]]$log_lambda,
            ibm_retin_control[[5]]$log_lambda)),
    t(rbind(ibm_retin_canesm5[[1]]$log_lambda,
            ibm_retin_canesm5[[2]]$log_lambda,
            ibm_retin_canesm5[[3]]$log_lambda,
            ibm_retin_canesm5[[4]]$log_lambda,
            ibm_retin_canesm5[[5]]$log_lambda)),
    t(rbind(ibm_retin_ec_earth3[[1]]$log_lambda,
            ibm_retin_ec_earth3[[2]]$log_lambda,
            ibm_retin_ec_earth3[[3]]$log_lambda,
            ibm_retin_ec_earth3[[4]]$log_lambda,
            ibm_retin_ec_earth3[[5]]$log_lambda)),
    t(rbind(ibm_retin_fgoals_g3[[1]]$log_lambda,
            ibm_retin_fgoals_g3[[2]]$log_lambda,
            ibm_retin_fgoals_g3[[3]]$log_lambda,
            ibm_retin_fgoals_g3[[4]]$log_lambda,
            ibm_retin_fgoals_g3[[5]]$log_lambda)),
    t(rbind(ibm_retin_gfdl_esm4[[1]]$log_lambda,
            ibm_retin_gfdl_esm4[[2]]$log_lambda,
            ibm_retin_gfdl_esm4[[3]]$log_lambda,
            ibm_retin_gfdl_esm4[[4]]$log_lambda,
            ibm_retin_gfdl_esm4[[5]]$log_lambda)),
    t(rbind(ibm_retin_giss_e2_1_g[[1]]$log_lambda,
            ibm_retin_giss_e2_1_g[[2]]$log_lambda,
            ibm_retin_giss_e2_1_g[[3]]$log_lambda,
            ibm_retin_giss_e2_1_g[[4]]$log_lambda,
            ibm_retin_giss_e2_1_g[[5]]$log_lambda)),
    t(rbind(ibm_retin_inm_cm4_8[[1]]$log_lambda,
            ibm_retin_inm_cm4_8[[2]]$log_lambda,
            ibm_retin_inm_cm4_8[[3]]$log_lambda,
            ibm_retin_inm_cm4_8[[4]]$log_lambda,
            ibm_retin_inm_cm4_8[[5]]$log_lambda)),
    t(rbind(ibm_retin_ipsl_cm6a_lr[[1]]$log_lambda,
            ibm_retin_ipsl_cm6a_lr[[2]]$log_lambda,
            ibm_retin_ipsl_cm6a_lr[[3]]$log_lambda,
            ibm_retin_ipsl_cm6a_lr[[4]]$log_lambda,
            ibm_retin_ipsl_cm6a_lr[[5]]$log_lambda)),
    t(rbind(ibm_retin_miroc6[[1]]$log_lambda,
            ibm_retin_miroc6[[2]]$log_lambda,
            ibm_retin_miroc6[[3]]$log_lambda,
            ibm_retin_miroc6[[4]]$log_lambda,
            ibm_retin_miroc6[[5]]$log_lambda)),
    t(rbind(ibm_retin_mpi_esm1_2_lr[[1]]$log_lambda,
            ibm_retin_mpi_esm1_2_lr[[2]]$log_lambda,
            ibm_retin_mpi_esm1_2_lr[[3]]$log_lambda,
            ibm_retin_mpi_esm1_2_lr[[4]]$log_lambda,
            ibm_retin_mpi_esm1_2_lr[[5]]$log_lambda)),
    t(rbind(ibm_retin_mri_esm2_0[[1]]$log_lambda,
            ibm_retin_mri_esm2_0[[2]]$log_lambda,
            ibm_retin_mri_esm2_0[[3]]$log_lambda,
            ibm_retin_mri_esm2_0[[4]]$log_lambda,
            ibm_retin_mri_esm2_0[[5]]$log_lambda)),
    t(rbind(ibm_retin_noresm2_mm[[1]]$log_lambda,
            ibm_retin_noresm2_mm[[2]]$log_lambda,
            ibm_retin_noresm2_mm[[3]]$log_lambda,
            ibm_retin_noresm2_mm[[4]]$log_lambda,
            ibm_retin_noresm2_mm[[5]]$log_lambda)))


# Add extinction

ibm_results$extinction = NA

ibm_results$extinction[which(ibm_results$population == "Retin")] = 
  c(rep(c(ibm_retin_control[[1]]$extinction,
          ibm_retin_control[[2]]$extinction,
          ibm_retin_control[[3]]$extinction,
          ibm_retin_control[[4]]$extinction,
          ibm_retin_control[[5]]$extinction), each = n_years),
    rep(c(ibm_retin_canesm5[[1]]$extinction,
          ibm_retin_canesm5[[2]]$extinction,
          ibm_retin_canesm5[[3]]$extinction,
          ibm_retin_canesm5[[4]]$extinction,
          ibm_retin_canesm5[[5]]$extinction), each = n_years),
    rep(c(ibm_retin_ec_earth3[[1]]$extinction,
          ibm_retin_ec_earth3[[2]]$extinction,
          ibm_retin_ec_earth3[[3]]$extinction,
          ibm_retin_ec_earth3[[4]]$extinction,
          ibm_retin_ec_earth3[[5]]$extinction), each = n_years),
    rep(c(ibm_retin_fgoals_g3[[1]]$extinction,
          ibm_retin_fgoals_g3[[2]]$extinction,
          ibm_retin_fgoals_g3[[3]]$extinction,
          ibm_retin_fgoals_g3[[4]]$extinction,
          ibm_retin_fgoals_g3[[5]]$extinction), each = n_years),
    rep(c(ibm_retin_gfdl_esm4[[1]]$extinction,
          ibm_retin_gfdl_esm4[[2]]$extinction,
          ibm_retin_gfdl_esm4[[3]]$extinction,
          ibm_retin_gfdl_esm4[[4]]$extinction,
          ibm_retin_gfdl_esm4[[5]]$extinction), each = n_years),
    rep(c(ibm_retin_giss_e2_1_g[[1]]$extinction,
          ibm_retin_giss_e2_1_g[[2]]$extinction,
          ibm_retin_giss_e2_1_g[[3]]$extinction,
          ibm_retin_giss_e2_1_g[[4]]$extinction,
          ibm_retin_giss_e2_1_g[[5]]$extinction), each = n_years),
    rep(c(ibm_retin_inm_cm4_8[[1]]$extinction,
          ibm_retin_inm_cm4_8[[2]]$extinction,
          ibm_retin_inm_cm4_8[[3]]$extinction,
          ibm_retin_inm_cm4_8[[4]]$extinction,
          ibm_retin_inm_cm4_8[[5]]$extinction), each = n_years),
    rep(c(ibm_retin_ipsl_cm6a_lr[[1]]$extinction,
          ibm_retin_ipsl_cm6a_lr[[2]]$extinction,
          ibm_retin_ipsl_cm6a_lr[[3]]$extinction,
          ibm_retin_ipsl_cm6a_lr[[4]]$extinction,
          ibm_retin_ipsl_cm6a_lr[[5]]$extinction), each = n_years),
    rep(c(ibm_retin_miroc6[[1]]$extinction,
          ibm_retin_miroc6[[2]]$extinction,
          ibm_retin_miroc6[[3]]$extinction,
          ibm_retin_miroc6[[4]]$extinction,
          ibm_retin_miroc6[[5]]$extinction), each = n_years),
    rep(c(ibm_retin_mpi_esm1_2_lr[[1]]$extinction,
          ibm_retin_mpi_esm1_2_lr[[2]]$extinction,
          ibm_retin_mpi_esm1_2_lr[[3]]$extinction,
          ibm_retin_mpi_esm1_2_lr[[4]]$extinction,
          ibm_retin_mpi_esm1_2_lr[[5]]$extinction), each = n_years),
    rep(c(ibm_retin_mri_esm2_0[[1]]$extinction,
          ibm_retin_mri_esm2_0[[2]]$extinction,
          ibm_retin_mri_esm2_0[[3]]$extinction,
          ibm_retin_mri_esm2_0[[4]]$extinction,
          ibm_retin_mri_esm2_0[[5]]$extinction), each = n_years),
    rep(c(ibm_retin_noresm2_mm[[1]]$extinction,
          ibm_retin_noresm2_mm[[2]]$extinction,
          ibm_retin_noresm2_mm[[3]]$extinction,
          ibm_retin_noresm2_mm[[4]]$extinction,
          ibm_retin_noresm2_mm[[5]]$extinction), each = n_years))


rm(list = grep("ibm_retin", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_Prisoneros", list.files("Output/Projections/"))])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_Prisoneros", list.files("Output/Projections/"))][i]))
}


# Add log lambda (per simulation = per row)

ibm_results$log_lambda[which(ibm_results$population == "Prisoneros")] = 
  c(t(rbind(ibm_prisoneros_control[[1]]$log_lambda,
            ibm_prisoneros_control[[2]]$log_lambda,
            ibm_prisoneros_control[[3]]$log_lambda,
            ibm_prisoneros_control[[4]]$log_lambda,
            ibm_prisoneros_control[[5]]$log_lambda)),
    t(rbind(ibm_prisoneros_canesm5[[1]]$log_lambda,
            ibm_prisoneros_canesm5[[2]]$log_lambda,
            ibm_prisoneros_canesm5[[3]]$log_lambda,
            ibm_prisoneros_canesm5[[4]]$log_lambda,
            ibm_prisoneros_canesm5[[5]]$log_lambda)),
    t(rbind(ibm_prisoneros_ec_earth3[[1]]$log_lambda,
            ibm_prisoneros_ec_earth3[[2]]$log_lambda,
            ibm_prisoneros_ec_earth3[[3]]$log_lambda,
            ibm_prisoneros_ec_earth3[[4]]$log_lambda,
            ibm_prisoneros_ec_earth3[[5]]$log_lambda)),
    t(rbind(ibm_prisoneros_fgoals_g3[[1]]$log_lambda,
            ibm_prisoneros_fgoals_g3[[2]]$log_lambda,
            ibm_prisoneros_fgoals_g3[[3]]$log_lambda,
            ibm_prisoneros_fgoals_g3[[4]]$log_lambda,
            ibm_prisoneros_fgoals_g3[[5]]$log_lambda)),
    t(rbind(ibm_prisoneros_gfdl_esm4[[1]]$log_lambda,
            ibm_prisoneros_gfdl_esm4[[2]]$log_lambda,
            ibm_prisoneros_gfdl_esm4[[3]]$log_lambda,
            ibm_prisoneros_gfdl_esm4[[4]]$log_lambda,
            ibm_prisoneros_gfdl_esm4[[5]]$log_lambda)),
    t(rbind(ibm_prisoneros_giss_e2_1_g[[1]]$log_lambda,
            ibm_prisoneros_giss_e2_1_g[[2]]$log_lambda,
            ibm_prisoneros_giss_e2_1_g[[3]]$log_lambda,
            ibm_prisoneros_giss_e2_1_g[[4]]$log_lambda,
            ibm_prisoneros_giss_e2_1_g[[5]]$log_lambda)),
    t(rbind(ibm_prisoneros_inm_cm4_8[[1]]$log_lambda,
            ibm_prisoneros_inm_cm4_8[[2]]$log_lambda,
            ibm_prisoneros_inm_cm4_8[[3]]$log_lambda,
            ibm_prisoneros_inm_cm4_8[[4]]$log_lambda,
            ibm_prisoneros_inm_cm4_8[[5]]$log_lambda)),
    t(rbind(ibm_prisoneros_ipsl_cm6a_lr[[1]]$log_lambda,
            ibm_prisoneros_ipsl_cm6a_lr[[2]]$log_lambda,
            ibm_prisoneros_ipsl_cm6a_lr[[3]]$log_lambda,
            ibm_prisoneros_ipsl_cm6a_lr[[4]]$log_lambda,
            ibm_prisoneros_ipsl_cm6a_lr[[5]]$log_lambda)),
    t(rbind(ibm_prisoneros_miroc6[[1]]$log_lambda,
            ibm_prisoneros_miroc6[[2]]$log_lambda,
            ibm_prisoneros_miroc6[[3]]$log_lambda,
            ibm_prisoneros_miroc6[[4]]$log_lambda,
            ibm_prisoneros_miroc6[[5]]$log_lambda)),
    t(rbind(ibm_prisoneros_mpi_esm1_2_lr[[1]]$log_lambda,
            ibm_prisoneros_mpi_esm1_2_lr[[2]]$log_lambda,
            ibm_prisoneros_mpi_esm1_2_lr[[3]]$log_lambda,
            ibm_prisoneros_mpi_esm1_2_lr[[4]]$log_lambda,
            ibm_prisoneros_mpi_esm1_2_lr[[5]]$log_lambda)),
    t(rbind(ibm_prisoneros_mri_esm2_0[[1]]$log_lambda,
            ibm_prisoneros_mri_esm2_0[[2]]$log_lambda,
            ibm_prisoneros_mri_esm2_0[[3]]$log_lambda,
            ibm_prisoneros_mri_esm2_0[[4]]$log_lambda,
            ibm_prisoneros_mri_esm2_0[[5]]$log_lambda)),
    t(rbind(ibm_prisoneros_noresm2_mm[[1]]$log_lambda,
            ibm_prisoneros_noresm2_mm[[2]]$log_lambda,
            ibm_prisoneros_noresm2_mm[[3]]$log_lambda,
            ibm_prisoneros_noresm2_mm[[4]]$log_lambda,
            ibm_prisoneros_noresm2_mm[[5]]$log_lambda)))


# Add extinction

ibm_results$extinction[which(ibm_results$population == "Prisoneros")] = 
  c(rep(c(ibm_prisoneros_control[[1]]$extinction,
          ibm_prisoneros_control[[2]]$extinction,
          ibm_prisoneros_control[[3]]$extinction,
          ibm_prisoneros_control[[4]]$extinction,
          ibm_prisoneros_control[[5]]$extinction), each = n_years),
    rep(c(ibm_prisoneros_canesm5[[1]]$extinction,
          ibm_prisoneros_canesm5[[2]]$extinction,
          ibm_prisoneros_canesm5[[3]]$extinction,
          ibm_prisoneros_canesm5[[4]]$extinction,
          ibm_prisoneros_canesm5[[5]]$extinction), each = n_years),
    rep(c(ibm_prisoneros_ec_earth3[[1]]$extinction,
          ibm_prisoneros_ec_earth3[[2]]$extinction,
          ibm_prisoneros_ec_earth3[[3]]$extinction,
          ibm_prisoneros_ec_earth3[[4]]$extinction,
          ibm_prisoneros_ec_earth3[[5]]$extinction), each = n_years),
    rep(c(ibm_prisoneros_fgoals_g3[[1]]$extinction,
          ibm_prisoneros_fgoals_g3[[2]]$extinction,
          ibm_prisoneros_fgoals_g3[[3]]$extinction,
          ibm_prisoneros_fgoals_g3[[4]]$extinction,
          ibm_prisoneros_fgoals_g3[[5]]$extinction), each = n_years),
    rep(c(ibm_prisoneros_gfdl_esm4[[1]]$extinction,
          ibm_prisoneros_gfdl_esm4[[2]]$extinction,
          ibm_prisoneros_gfdl_esm4[[3]]$extinction,
          ibm_prisoneros_gfdl_esm4[[4]]$extinction,
          ibm_prisoneros_gfdl_esm4[[5]]$extinction), each = n_years),
    rep(c(ibm_prisoneros_giss_e2_1_g[[1]]$extinction,
          ibm_prisoneros_giss_e2_1_g[[2]]$extinction,
          ibm_prisoneros_giss_e2_1_g[[3]]$extinction,
          ibm_prisoneros_giss_e2_1_g[[4]]$extinction,
          ibm_prisoneros_giss_e2_1_g[[5]]$extinction), each = n_years),
    rep(c(ibm_prisoneros_inm_cm4_8[[1]]$extinction,
          ibm_prisoneros_inm_cm4_8[[2]]$extinction,
          ibm_prisoneros_inm_cm4_8[[3]]$extinction,
          ibm_prisoneros_inm_cm4_8[[4]]$extinction,
          ibm_prisoneros_inm_cm4_8[[5]]$extinction), each = n_years),
    rep(c(ibm_prisoneros_ipsl_cm6a_lr[[1]]$extinction,
          ibm_prisoneros_ipsl_cm6a_lr[[2]]$extinction,
          ibm_prisoneros_ipsl_cm6a_lr[[3]]$extinction,
          ibm_prisoneros_ipsl_cm6a_lr[[4]]$extinction,
          ibm_prisoneros_ipsl_cm6a_lr[[5]]$extinction), each = n_years),
    rep(c(ibm_prisoneros_miroc6[[1]]$extinction,
          ibm_prisoneros_miroc6[[2]]$extinction,
          ibm_prisoneros_miroc6[[3]]$extinction,
          ibm_prisoneros_miroc6[[4]]$extinction,
          ibm_prisoneros_miroc6[[5]]$extinction), each = n_years),
    rep(c(ibm_prisoneros_mpi_esm1_2_lr[[1]]$extinction,
          ibm_prisoneros_mpi_esm1_2_lr[[2]]$extinction,
          ibm_prisoneros_mpi_esm1_2_lr[[3]]$extinction,
          ibm_prisoneros_mpi_esm1_2_lr[[4]]$extinction,
          ibm_prisoneros_mpi_esm1_2_lr[[5]]$extinction), each = n_years),
    rep(c(ibm_prisoneros_mri_esm2_0[[1]]$extinction,
          ibm_prisoneros_mri_esm2_0[[2]]$extinction,
          ibm_prisoneros_mri_esm2_0[[3]]$extinction,
          ibm_prisoneros_mri_esm2_0[[4]]$extinction,
          ibm_prisoneros_mri_esm2_0[[5]]$extinction), each = n_years),
    rep(c(ibm_prisoneros_noresm2_mm[[1]]$extinction,
          ibm_prisoneros_noresm2_mm[[2]]$extinction,
          ibm_prisoneros_noresm2_mm[[3]]$extinction,
          ibm_prisoneros_noresm2_mm[[4]]$extinction,
          ibm_prisoneros_noresm2_mm[[5]]$extinction), each = n_years))


rm(list = grep("ibm_prisoneros", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_Bujeo", list.files("Output/Projections/"))])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_Bujeo", list.files("Output/Projections/"))][i]))
}


# Add log lambda (per simulation = per row)

ibm_results$log_lambda[which(ibm_results$population == "Bujeo")] = 
  c(t(rbind(ibm_bujeo_control[[1]]$log_lambda,
            ibm_bujeo_control[[2]]$log_lambda,
            ibm_bujeo_control[[3]]$log_lambda,
            ibm_bujeo_control[[4]]$log_lambda,
            ibm_bujeo_control[[5]]$log_lambda)),
    t(rbind(ibm_bujeo_canesm5[[1]]$log_lambda,
            ibm_bujeo_canesm5[[2]]$log_lambda,
            ibm_bujeo_canesm5[[3]]$log_lambda,
            ibm_bujeo_canesm5[[4]]$log_lambda,
            ibm_bujeo_canesm5[[5]]$log_lambda)),
    t(rbind(ibm_bujeo_ec_earth3[[1]]$log_lambda,
            ibm_bujeo_ec_earth3[[2]]$log_lambda,
            ibm_bujeo_ec_earth3[[3]]$log_lambda,
            ibm_bujeo_ec_earth3[[4]]$log_lambda,
            ibm_bujeo_ec_earth3[[5]]$log_lambda)),
    t(rbind(ibm_bujeo_fgoals_g3[[1]]$log_lambda,
            ibm_bujeo_fgoals_g3[[2]]$log_lambda,
            ibm_bujeo_fgoals_g3[[3]]$log_lambda,
            ibm_bujeo_fgoals_g3[[4]]$log_lambda,
            ibm_bujeo_fgoals_g3[[5]]$log_lambda)),
    t(rbind(ibm_bujeo_gfdl_esm4[[1]]$log_lambda,
            ibm_bujeo_gfdl_esm4[[2]]$log_lambda,
            ibm_bujeo_gfdl_esm4[[3]]$log_lambda,
            ibm_bujeo_gfdl_esm4[[4]]$log_lambda,
            ibm_bujeo_gfdl_esm4[[5]]$log_lambda)),
    t(rbind(ibm_bujeo_giss_e2_1_g[[1]]$log_lambda,
            ibm_bujeo_giss_e2_1_g[[2]]$log_lambda,
            ibm_bujeo_giss_e2_1_g[[3]]$log_lambda,
            ibm_bujeo_giss_e2_1_g[[4]]$log_lambda,
            ibm_bujeo_giss_e2_1_g[[5]]$log_lambda)),
    t(rbind(ibm_bujeo_inm_cm4_8[[1]]$log_lambda,
            ibm_bujeo_inm_cm4_8[[2]]$log_lambda,
            ibm_bujeo_inm_cm4_8[[3]]$log_lambda,
            ibm_bujeo_inm_cm4_8[[4]]$log_lambda,
            ibm_bujeo_inm_cm4_8[[5]]$log_lambda)),
    t(rbind(ibm_bujeo_ipsl_cm6a_lr[[1]]$log_lambda,
            ibm_bujeo_ipsl_cm6a_lr[[2]]$log_lambda,
            ibm_bujeo_ipsl_cm6a_lr[[3]]$log_lambda,
            ibm_bujeo_ipsl_cm6a_lr[[4]]$log_lambda,
            ibm_bujeo_ipsl_cm6a_lr[[5]]$log_lambda)),
    t(rbind(ibm_bujeo_miroc6[[1]]$log_lambda,
            ibm_bujeo_miroc6[[2]]$log_lambda,
            ibm_bujeo_miroc6[[3]]$log_lambda,
            ibm_bujeo_miroc6[[4]]$log_lambda,
            ibm_bujeo_miroc6[[5]]$log_lambda)),
    t(rbind(ibm_bujeo_mpi_esm1_2_lr[[1]]$log_lambda,
            ibm_bujeo_mpi_esm1_2_lr[[2]]$log_lambda,
            ibm_bujeo_mpi_esm1_2_lr[[3]]$log_lambda,
            ibm_bujeo_mpi_esm1_2_lr[[4]]$log_lambda,
            ibm_bujeo_mpi_esm1_2_lr[[5]]$log_lambda)),
    t(rbind(ibm_bujeo_mri_esm2_0[[1]]$log_lambda,
            ibm_bujeo_mri_esm2_0[[2]]$log_lambda,
            ibm_bujeo_mri_esm2_0[[3]]$log_lambda,
            ibm_bujeo_mri_esm2_0[[4]]$log_lambda,
            ibm_bujeo_mri_esm2_0[[5]]$log_lambda)),
    t(rbind(ibm_bujeo_noresm2_mm[[1]]$log_lambda,
            ibm_bujeo_noresm2_mm[[2]]$log_lambda,
            ibm_bujeo_noresm2_mm[[3]]$log_lambda,
            ibm_bujeo_noresm2_mm[[4]]$log_lambda,
            ibm_bujeo_noresm2_mm[[5]]$log_lambda)))


# Add extinction

ibm_results$extinction[which(ibm_results$population == "Bujeo")] = 
  c(rep(c(ibm_bujeo_control[[1]]$extinction,
          ibm_bujeo_control[[2]]$extinction,
          ibm_bujeo_control[[3]]$extinction,
          ibm_bujeo_control[[4]]$extinction,
          ibm_bujeo_control[[5]]$extinction), each = n_years),
    rep(c(ibm_bujeo_canesm5[[1]]$extinction,
          ibm_bujeo_canesm5[[2]]$extinction,
          ibm_bujeo_canesm5[[3]]$extinction,
          ibm_bujeo_canesm5[[4]]$extinction,
          ibm_bujeo_canesm5[[5]]$extinction), each = n_years),
    rep(c(ibm_bujeo_ec_earth3[[1]]$extinction,
          ibm_bujeo_ec_earth3[[2]]$extinction,
          ibm_bujeo_ec_earth3[[3]]$extinction,
          ibm_bujeo_ec_earth3[[4]]$extinction,
          ibm_bujeo_ec_earth3[[5]]$extinction), each = n_years),
    rep(c(ibm_bujeo_fgoals_g3[[1]]$extinction,
          ibm_bujeo_fgoals_g3[[2]]$extinction,
          ibm_bujeo_fgoals_g3[[3]]$extinction,
          ibm_bujeo_fgoals_g3[[4]]$extinction,
          ibm_bujeo_fgoals_g3[[5]]$extinction), each = n_years),
    rep(c(ibm_bujeo_gfdl_esm4[[1]]$extinction,
          ibm_bujeo_gfdl_esm4[[2]]$extinction,
          ibm_bujeo_gfdl_esm4[[3]]$extinction,
          ibm_bujeo_gfdl_esm4[[4]]$extinction,
          ibm_bujeo_gfdl_esm4[[5]]$extinction), each = n_years),
    rep(c(ibm_bujeo_giss_e2_1_g[[1]]$extinction,
          ibm_bujeo_giss_e2_1_g[[2]]$extinction,
          ibm_bujeo_giss_e2_1_g[[3]]$extinction,
          ibm_bujeo_giss_e2_1_g[[4]]$extinction,
          ibm_bujeo_giss_e2_1_g[[5]]$extinction), each = n_years),
    rep(c(ibm_bujeo_inm_cm4_8[[1]]$extinction,
          ibm_bujeo_inm_cm4_8[[2]]$extinction,
          ibm_bujeo_inm_cm4_8[[3]]$extinction,
          ibm_bujeo_inm_cm4_8[[4]]$extinction,
          ibm_bujeo_inm_cm4_8[[5]]$extinction), each = n_years),
    rep(c(ibm_bujeo_ipsl_cm6a_lr[[1]]$extinction,
          ibm_bujeo_ipsl_cm6a_lr[[2]]$extinction,
          ibm_bujeo_ipsl_cm6a_lr[[3]]$extinction,
          ibm_bujeo_ipsl_cm6a_lr[[4]]$extinction,
          ibm_bujeo_ipsl_cm6a_lr[[5]]$extinction), each = n_years),
    rep(c(ibm_bujeo_miroc6[[1]]$extinction,
          ibm_bujeo_miroc6[[2]]$extinction,
          ibm_bujeo_miroc6[[3]]$extinction,
          ibm_bujeo_miroc6[[4]]$extinction,
          ibm_bujeo_miroc6[[5]]$extinction), each = n_years),
    rep(c(ibm_bujeo_mpi_esm1_2_lr[[1]]$extinction,
          ibm_bujeo_mpi_esm1_2_lr[[2]]$extinction,
          ibm_bujeo_mpi_esm1_2_lr[[3]]$extinction,
          ibm_bujeo_mpi_esm1_2_lr[[4]]$extinction,
          ibm_bujeo_mpi_esm1_2_lr[[5]]$extinction), each = n_years),
    rep(c(ibm_bujeo_mri_esm2_0[[1]]$extinction,
          ibm_bujeo_mri_esm2_0[[2]]$extinction,
          ibm_bujeo_mri_esm2_0[[3]]$extinction,
          ibm_bujeo_mri_esm2_0[[4]]$extinction,
          ibm_bujeo_mri_esm2_0[[5]]$extinction), each = n_years),
    rep(c(ibm_bujeo_noresm2_mm[[1]]$extinction,
          ibm_bujeo_noresm2_mm[[2]]$extinction,
          ibm_bujeo_noresm2_mm[[3]]$extinction,
          ibm_bujeo_noresm2_mm[[4]]$extinction,
          ibm_bujeo_noresm2_mm[[5]]$extinction), each = n_years))


rm(list = grep("ibm_bujeo", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_MonteraTorero", list.files("Output/Projections/"))])){

  print(i)

  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_MonteraTorero", list.files("Output/Projections/"))][i]))
}


# Add log lambda (per simulation = per row)

ibm_results$log_lambda[which(ibm_results$population == "MonteraTorero")] =
  c(t(rbind(ibm_monteratorero_control[[1]]$log_lambda,
            ibm_monteratorero_control[[2]]$log_lambda,
            ibm_monteratorero_control[[3]]$log_lambda,
            ibm_monteratorero_control[[4]]$log_lambda,
            ibm_monteratorero_control[[5]]$log_lambda)),
    t(rbind(ibm_monteratorero_canesm5[[1]]$log_lambda,
            ibm_monteratorero_canesm5[[2]]$log_lambda,
            ibm_monteratorero_canesm5[[3]]$log_lambda,
            ibm_monteratorero_canesm5[[4]]$log_lambda,
            ibm_monteratorero_canesm5[[5]]$log_lambda)),
    t(rbind(ibm_monteratorero_ec_earth3[[1]]$log_lambda,
            ibm_monteratorero_ec_earth3[[2]]$log_lambda,
            ibm_monteratorero_ec_earth3[[3]]$log_lambda,
            ibm_monteratorero_ec_earth3[[4]]$log_lambda,
            ibm_monteratorero_ec_earth3[[5]]$log_lambda)),
    t(rbind(ibm_monteratorero_fgoals_g3[[1]]$log_lambda,
            ibm_monteratorero_fgoals_g3[[2]]$log_lambda,
            ibm_monteratorero_fgoals_g3[[3]]$log_lambda,
            ibm_monteratorero_fgoals_g3[[4]]$log_lambda,
            ibm_monteratorero_fgoals_g3[[5]]$log_lambda)),
    t(rbind(ibm_monteratorero_gfdl_esm4[[1]]$log_lambda,
            ibm_monteratorero_gfdl_esm4[[2]]$log_lambda,
            ibm_monteratorero_gfdl_esm4[[3]]$log_lambda,
            ibm_monteratorero_gfdl_esm4[[4]]$log_lambda,
            ibm_monteratorero_gfdl_esm4[[5]]$log_lambda)),
    t(rbind(ibm_monteratorero_giss_e2_1_g[[1]]$log_lambda,
            ibm_monteratorero_giss_e2_1_g[[2]]$log_lambda,
            ibm_monteratorero_giss_e2_1_g[[3]]$log_lambda,
            ibm_monteratorero_giss_e2_1_g[[4]]$log_lambda,
            ibm_monteratorero_giss_e2_1_g[[5]]$log_lambda)),
    t(rbind(ibm_monteratorero_inm_cm4_8[[1]]$log_lambda,
            ibm_monteratorero_inm_cm4_8[[2]]$log_lambda,
            ibm_monteratorero_inm_cm4_8[[3]]$log_lambda,
            ibm_monteratorero_inm_cm4_8[[4]]$log_lambda,
            ibm_monteratorero_inm_cm4_8[[5]]$log_lambda)),
    t(rbind(ibm_monteratorero_ipsl_cm6a_lr[[1]]$log_lambda,
            ibm_monteratorero_ipsl_cm6a_lr[[2]]$log_lambda,
            ibm_monteratorero_ipsl_cm6a_lr[[3]]$log_lambda,
            ibm_monteratorero_ipsl_cm6a_lr[[4]]$log_lambda,
            ibm_monteratorero_ipsl_cm6a_lr[[5]]$log_lambda)),
    t(rbind(ibm_monteratorero_miroc6[[1]]$log_lambda,
            ibm_monteratorero_miroc6[[2]]$log_lambda,
            ibm_monteratorero_miroc6[[3]]$log_lambda,
            ibm_monteratorero_miroc6[[4]]$log_lambda,
            ibm_monteratorero_miroc6[[5]]$log_lambda)),
    t(rbind(ibm_monteratorero_mpi_esm1_2_lr[[1]]$log_lambda,
            ibm_monteratorero_mpi_esm1_2_lr[[2]]$log_lambda,
            ibm_monteratorero_mpi_esm1_2_lr[[3]]$log_lambda,
            ibm_monteratorero_mpi_esm1_2_lr[[4]]$log_lambda,
            ibm_monteratorero_mpi_esm1_2_lr[[5]]$log_lambda)),
    t(rbind(ibm_monteratorero_mri_esm2_0[[1]]$log_lambda,
            ibm_monteratorero_mri_esm2_0[[2]]$log_lambda,
            ibm_monteratorero_mri_esm2_0[[3]]$log_lambda,
            ibm_monteratorero_mri_esm2_0[[4]]$log_lambda,
            ibm_monteratorero_mri_esm2_0[[5]]$log_lambda)),
    t(rbind(ibm_monteratorero_noresm2_mm[[1]]$log_lambda,
            ibm_monteratorero_noresm2_mm[[2]]$log_lambda,
            ibm_monteratorero_noresm2_mm[[3]]$log_lambda,
            ibm_monteratorero_noresm2_mm[[4]]$log_lambda,
            ibm_monteratorero_noresm2_mm[[5]]$log_lambda)))


# Add extinction

ibm_results$extinction[which(ibm_results$population == "MonteraTorero")] =
  c(rep(c(ibm_monteratorero_control[[1]]$extinction,
          ibm_monteratorero_control[[2]]$extinction,
          ibm_monteratorero_control[[3]]$extinction,
          ibm_monteratorero_control[[4]]$extinction,
          ibm_monteratorero_control[[5]]$extinction), each = n_years),
    rep(c(ibm_monteratorero_canesm5[[1]]$extinction,
          ibm_monteratorero_canesm5[[2]]$extinction,
          ibm_monteratorero_canesm5[[3]]$extinction,
          ibm_monteratorero_canesm5[[4]]$extinction,
          ibm_monteratorero_canesm5[[5]]$extinction), each = n_years),
    rep(c(ibm_monteratorero_ec_earth3[[1]]$extinction,
          ibm_monteratorero_ec_earth3[[2]]$extinction,
          ibm_monteratorero_ec_earth3[[3]]$extinction,
          ibm_monteratorero_ec_earth3[[4]]$extinction,
          ibm_monteratorero_ec_earth3[[5]]$extinction), each = n_years),
    rep(c(ibm_monteratorero_fgoals_g3[[1]]$extinction,
          ibm_monteratorero_fgoals_g3[[2]]$extinction,
          ibm_monteratorero_fgoals_g3[[3]]$extinction,
          ibm_monteratorero_fgoals_g3[[4]]$extinction,
          ibm_monteratorero_fgoals_g3[[5]]$extinction), each = n_years),
    rep(c(ibm_monteratorero_gfdl_esm4[[1]]$extinction,
          ibm_monteratorero_gfdl_esm4[[2]]$extinction,
          ibm_monteratorero_gfdl_esm4[[3]]$extinction,
          ibm_monteratorero_gfdl_esm4[[4]]$extinction,
          ibm_monteratorero_gfdl_esm4[[5]]$extinction), each = n_years),
    rep(c(ibm_monteratorero_giss_e2_1_g[[1]]$extinction,
          ibm_monteratorero_giss_e2_1_g[[2]]$extinction,
          ibm_monteratorero_giss_e2_1_g[[3]]$extinction,
          ibm_monteratorero_giss_e2_1_g[[4]]$extinction,
          ibm_monteratorero_giss_e2_1_g[[5]]$extinction), each = n_years),
    rep(c(ibm_monteratorero_inm_cm4_8[[1]]$extinction,
          ibm_monteratorero_inm_cm4_8[[2]]$extinction,
          ibm_monteratorero_inm_cm4_8[[3]]$extinction,
          ibm_monteratorero_inm_cm4_8[[4]]$extinction,
          ibm_monteratorero_inm_cm4_8[[5]]$extinction), each = n_years),
    rep(c(ibm_monteratorero_ipsl_cm6a_lr[[1]]$extinction,
          ibm_monteratorero_ipsl_cm6a_lr[[2]]$extinction,
          ibm_monteratorero_ipsl_cm6a_lr[[3]]$extinction,
          ibm_monteratorero_ipsl_cm6a_lr[[4]]$extinction,
          ibm_monteratorero_ipsl_cm6a_lr[[5]]$extinction), each = n_years),
    rep(c(ibm_monteratorero_miroc6[[1]]$extinction,
          ibm_monteratorero_miroc6[[2]]$extinction,
          ibm_monteratorero_miroc6[[3]]$extinction,
          ibm_monteratorero_miroc6[[4]]$extinction,
          ibm_monteratorero_miroc6[[5]]$extinction), each = n_years),
    rep(c(ibm_monteratorero_mpi_esm1_2_lr[[1]]$extinction,
          ibm_monteratorero_mpi_esm1_2_lr[[2]]$extinction,
          ibm_monteratorero_mpi_esm1_2_lr[[3]]$extinction,
          ibm_monteratorero_mpi_esm1_2_lr[[4]]$extinction,
          ibm_monteratorero_mpi_esm1_2_lr[[5]]$extinction), each = n_years),
    rep(c(ibm_monteratorero_mri_esm2_0[[1]]$extinction,
          ibm_monteratorero_mri_esm2_0[[2]]$extinction,
          ibm_monteratorero_mri_esm2_0[[3]]$extinction,
          ibm_monteratorero_mri_esm2_0[[4]]$extinction,
          ibm_monteratorero_mri_esm2_0[[5]]$extinction), each = n_years),
    rep(c(ibm_monteratorero_noresm2_mm[[1]]$extinction,
          ibm_monteratorero_noresm2_mm[[2]]$extinction,
          ibm_monteratorero_noresm2_mm[[3]]$extinction,
          ibm_monteratorero_noresm2_mm[[4]]$extinction,
          ibm_monteratorero_noresm2_mm[[5]]$extinction), each = n_years))


rm(list = grep("ibm_monteratorero", ls() , value = TRUE, invert = FALSE))


for(i in 1:length(list.files("Output/Projections/")[grep("Anthropogenic_SCarbDist", list.files("Output/Projections/"))])){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Anthropogenic_SCarbDist", list.files("Output/Projections/"))][i]))
}


# Add log lambda (per simulation = per row)

ibm_results$log_lambda[which(ibm_results$population == "SCarbDist")] = 
  c(t(rbind(ibm_scarbdist_control[[1]]$log_lambda,
            ibm_scarbdist_control[[2]]$log_lambda,
            ibm_scarbdist_control[[3]]$log_lambda,
            ibm_scarbdist_control[[4]]$log_lambda,
            ibm_scarbdist_control[[5]]$log_lambda)),
    t(rbind(ibm_scarbdist_canesm5[[1]]$log_lambda,
            ibm_scarbdist_canesm5[[2]]$log_lambda,
            ibm_scarbdist_canesm5[[3]]$log_lambda,
            ibm_scarbdist_canesm5[[4]]$log_lambda,
            ibm_scarbdist_canesm5[[5]]$log_lambda)),
    t(rbind(ibm_scarbdist_ec_earth3[[1]]$log_lambda,
            ibm_scarbdist_ec_earth3[[2]]$log_lambda,
            ibm_scarbdist_ec_earth3[[3]]$log_lambda,
            ibm_scarbdist_ec_earth3[[4]]$log_lambda,
            ibm_scarbdist_ec_earth3[[5]]$log_lambda)),
    t(rbind(ibm_scarbdist_fgoals_g3[[1]]$log_lambda,
            ibm_scarbdist_fgoals_g3[[2]]$log_lambda,
            ibm_scarbdist_fgoals_g3[[3]]$log_lambda,
            ibm_scarbdist_fgoals_g3[[4]]$log_lambda,
            ibm_scarbdist_fgoals_g3[[5]]$log_lambda)),
    t(rbind(ibm_scarbdist_gfdl_esm4[[1]]$log_lambda,
            ibm_scarbdist_gfdl_esm4[[2]]$log_lambda,
            ibm_scarbdist_gfdl_esm4[[3]]$log_lambda,
            ibm_scarbdist_gfdl_esm4[[4]]$log_lambda,
            ibm_scarbdist_gfdl_esm4[[5]]$log_lambda)),
    t(rbind(ibm_scarbdist_giss_e2_1_g[[1]]$log_lambda,
            ibm_scarbdist_giss_e2_1_g[[2]]$log_lambda,
            ibm_scarbdist_giss_e2_1_g[[3]]$log_lambda,
            ibm_scarbdist_giss_e2_1_g[[4]]$log_lambda,
            ibm_scarbdist_giss_e2_1_g[[5]]$log_lambda)),
    t(rbind(ibm_scarbdist_inm_cm4_8[[1]]$log_lambda,
            ibm_scarbdist_inm_cm4_8[[2]]$log_lambda,
            ibm_scarbdist_inm_cm4_8[[3]]$log_lambda,
            ibm_scarbdist_inm_cm4_8[[4]]$log_lambda,
            ibm_scarbdist_inm_cm4_8[[5]]$log_lambda)),
    t(rbind(ibm_scarbdist_ipsl_cm6a_lr[[1]]$log_lambda,
            ibm_scarbdist_ipsl_cm6a_lr[[2]]$log_lambda,
            ibm_scarbdist_ipsl_cm6a_lr[[3]]$log_lambda,
            ibm_scarbdist_ipsl_cm6a_lr[[4]]$log_lambda,
            ibm_scarbdist_ipsl_cm6a_lr[[5]]$log_lambda)),
    t(rbind(ibm_scarbdist_miroc6[[1]]$log_lambda,
            ibm_scarbdist_miroc6[[2]]$log_lambda,
            ibm_scarbdist_miroc6[[3]]$log_lambda,
            ibm_scarbdist_miroc6[[4]]$log_lambda,
            ibm_scarbdist_miroc6[[5]]$log_lambda)),
    t(rbind(ibm_scarbdist_mpi_esm1_2_lr[[1]]$log_lambda,
            ibm_scarbdist_mpi_esm1_2_lr[[2]]$log_lambda,
            ibm_scarbdist_mpi_esm1_2_lr[[3]]$log_lambda,
            ibm_scarbdist_mpi_esm1_2_lr[[4]]$log_lambda,
            ibm_scarbdist_mpi_esm1_2_lr[[5]]$log_lambda)),
    t(rbind(ibm_scarbdist_mri_esm2_0[[1]]$log_lambda,
            ibm_scarbdist_mri_esm2_0[[2]]$log_lambda,
            ibm_scarbdist_mri_esm2_0[[3]]$log_lambda,
            ibm_scarbdist_mri_esm2_0[[4]]$log_lambda,
            ibm_scarbdist_mri_esm2_0[[5]]$log_lambda)),
    t(rbind(ibm_scarbdist_noresm2_mm[[1]]$log_lambda,
            ibm_scarbdist_noresm2_mm[[2]]$log_lambda,
            ibm_scarbdist_noresm2_mm[[3]]$log_lambda,
            ibm_scarbdist_noresm2_mm[[4]]$log_lambda,
            ibm_scarbdist_noresm2_mm[[5]]$log_lambda)))


# Add extinction

ibm_results$extinction[which(ibm_results$population == "SCarbDist")] = 
  c(rep(c(ibm_scarbdist_control[[1]]$extinction,
          ibm_scarbdist_control[[2]]$extinction,
          ibm_scarbdist_control[[3]]$extinction,
          ibm_scarbdist_control[[4]]$extinction,
          ibm_scarbdist_control[[5]]$extinction), each = n_years),
    rep(c(ibm_scarbdist_canesm5[[1]]$extinction,
          ibm_scarbdist_canesm5[[2]]$extinction,
          ibm_scarbdist_canesm5[[3]]$extinction,
          ibm_scarbdist_canesm5[[4]]$extinction,
          ibm_scarbdist_canesm5[[5]]$extinction), each = n_years),
    rep(c(ibm_scarbdist_ec_earth3[[1]]$extinction,
          ibm_scarbdist_ec_earth3[[2]]$extinction,
          ibm_scarbdist_ec_earth3[[3]]$extinction,
          ibm_scarbdist_ec_earth3[[4]]$extinction,
          ibm_scarbdist_ec_earth3[[5]]$extinction), each = n_years),
    rep(c(ibm_scarbdist_fgoals_g3[[1]]$extinction,
          ibm_scarbdist_fgoals_g3[[2]]$extinction,
          ibm_scarbdist_fgoals_g3[[3]]$extinction,
          ibm_scarbdist_fgoals_g3[[4]]$extinction,
          ibm_scarbdist_fgoals_g3[[5]]$extinction), each = n_years),
    rep(c(ibm_scarbdist_gfdl_esm4[[1]]$extinction,
          ibm_scarbdist_gfdl_esm4[[2]]$extinction,
          ibm_scarbdist_gfdl_esm4[[3]]$extinction,
          ibm_scarbdist_gfdl_esm4[[4]]$extinction,
          ibm_scarbdist_gfdl_esm4[[5]]$extinction), each = n_years),
    rep(c(ibm_scarbdist_giss_e2_1_g[[1]]$extinction,
          ibm_scarbdist_giss_e2_1_g[[2]]$extinction,
          ibm_scarbdist_giss_e2_1_g[[3]]$extinction,
          ibm_scarbdist_giss_e2_1_g[[4]]$extinction,
          ibm_scarbdist_giss_e2_1_g[[5]]$extinction), each = n_years),
    rep(c(ibm_scarbdist_inm_cm4_8[[1]]$extinction,
          ibm_scarbdist_inm_cm4_8[[2]]$extinction,
          ibm_scarbdist_inm_cm4_8[[3]]$extinction,
          ibm_scarbdist_inm_cm4_8[[4]]$extinction,
          ibm_scarbdist_inm_cm4_8[[5]]$extinction), each = n_years),
    rep(c(ibm_scarbdist_ipsl_cm6a_lr[[1]]$extinction,
          ibm_scarbdist_ipsl_cm6a_lr[[2]]$extinction,
          ibm_scarbdist_ipsl_cm6a_lr[[3]]$extinction,
          ibm_scarbdist_ipsl_cm6a_lr[[4]]$extinction,
          ibm_scarbdist_ipsl_cm6a_lr[[5]]$extinction), each = n_years),
    rep(c(ibm_scarbdist_miroc6[[1]]$extinction,
          ibm_scarbdist_miroc6[[2]]$extinction,
          ibm_scarbdist_miroc6[[3]]$extinction,
          ibm_scarbdist_miroc6[[4]]$extinction,
          ibm_scarbdist_miroc6[[5]]$extinction), each = n_years),
    rep(c(ibm_scarbdist_mpi_esm1_2_lr[[1]]$extinction,
          ibm_scarbdist_mpi_esm1_2_lr[[2]]$extinction,
          ibm_scarbdist_mpi_esm1_2_lr[[3]]$extinction,
          ibm_scarbdist_mpi_esm1_2_lr[[4]]$extinction,
          ibm_scarbdist_mpi_esm1_2_lr[[5]]$extinction), each = n_years),
    rep(c(ibm_scarbdist_mri_esm2_0[[1]]$extinction,
          ibm_scarbdist_mri_esm2_0[[2]]$extinction,
          ibm_scarbdist_mri_esm2_0[[3]]$extinction,
          ibm_scarbdist_mri_esm2_0[[4]]$extinction,
          ibm_scarbdist_mri_esm2_0[[5]]$extinction), each = n_years),
    rep(c(ibm_scarbdist_noresm2_mm[[1]]$extinction,
          ibm_scarbdist_noresm2_mm[[2]]$extinction,
          ibm_scarbdist_noresm2_mm[[3]]$extinction,
          ibm_scarbdist_noresm2_mm[[4]]$extinction,
          ibm_scarbdist_noresm2_mm[[5]]$extinction), each = n_years))


rm(list = grep("ibm_scarbdist", ls() , value = TRUE, invert = FALSE))


for(i in 1:(length(list.files("Output/Projections/")[grep("Natural_SierraCarbonera", list.files("Output/Projections/"))]))){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SierraCarbonera", list.files("Output/Projections/"))][i]))
}


# Add log lambda (per simulation = per row)

ibm_results$log_lambda[which(ibm_results$population == "SierraCarboneraY5")] = 
  c(t(rbind(ibm_sierracarboneray5_control[[1]]$log_lambda,
            ibm_sierracarboneray5_control[[2]]$log_lambda,
            ibm_sierracarboneray5_control[[3]]$log_lambda,
            ibm_sierracarboneray5_control[[4]]$log_lambda,
            ibm_sierracarboneray5_control[[5]]$log_lambda)),
    t(rbind(ibm_sierracarboneray5_canesm5[[1]]$log_lambda,
            ibm_sierracarboneray5_canesm5[[2]]$log_lambda,
            ibm_sierracarboneray5_canesm5[[3]]$log_lambda,
            ibm_sierracarboneray5_canesm5[[4]]$log_lambda,
            ibm_sierracarboneray5_canesm5[[5]]$log_lambda)),
    t(rbind(ibm_sierracarboneray5_ec_earth3[[1]]$log_lambda,
            ibm_sierracarboneray5_ec_earth3[[2]]$log_lambda,
            ibm_sierracarboneray5_ec_earth3[[3]]$log_lambda,
            ibm_sierracarboneray5_ec_earth3[[4]]$log_lambda,
            ibm_sierracarboneray5_ec_earth3[[5]]$log_lambda)),
    t(rbind(ibm_sierracarboneray5_fgoals_g3[[1]]$log_lambda,
            ibm_sierracarboneray5_fgoals_g3[[2]]$log_lambda,
            ibm_sierracarboneray5_fgoals_g3[[3]]$log_lambda,
            ibm_sierracarboneray5_fgoals_g3[[4]]$log_lambda,
            ibm_sierracarboneray5_fgoals_g3[[5]]$log_lambda)),
    t(rbind(ibm_sierracarboneray5_gfdl_esm4[[1]]$log_lambda,
            ibm_sierracarboneray5_gfdl_esm4[[2]]$log_lambda,
            ibm_sierracarboneray5_gfdl_esm4[[3]]$log_lambda,
            ibm_sierracarboneray5_gfdl_esm4[[4]]$log_lambda,
            ibm_sierracarboneray5_gfdl_esm4[[5]]$log_lambda)),
    t(rbind(ibm_sierracarboneray5_giss_e2_1_g[[1]]$log_lambda,
            ibm_sierracarboneray5_giss_e2_1_g[[2]]$log_lambda,
            ibm_sierracarboneray5_giss_e2_1_g[[3]]$log_lambda,
            ibm_sierracarboneray5_giss_e2_1_g[[4]]$log_lambda,
            ibm_sierracarboneray5_giss_e2_1_g[[5]]$log_lambda)),
    t(rbind(ibm_sierracarboneray5_inm_cm4_8[[1]]$log_lambda,
            ibm_sierracarboneray5_inm_cm4_8[[2]]$log_lambda,
            ibm_sierracarboneray5_inm_cm4_8[[3]]$log_lambda,
            ibm_sierracarboneray5_inm_cm4_8[[4]]$log_lambda,
            ibm_sierracarboneray5_inm_cm4_8[[5]]$log_lambda)),
    t(rbind(ibm_sierracarboneray5_ipsl_cm6a_lr[[1]]$log_lambda,
            ibm_sierracarboneray5_ipsl_cm6a_lr[[2]]$log_lambda,
            ibm_sierracarboneray5_ipsl_cm6a_lr[[3]]$log_lambda,
            ibm_sierracarboneray5_ipsl_cm6a_lr[[4]]$log_lambda,
            ibm_sierracarboneray5_ipsl_cm6a_lr[[5]]$log_lambda)),
    t(rbind(ibm_sierracarboneray5_miroc6[[1]]$log_lambda,
            ibm_sierracarboneray5_miroc6[[2]]$log_lambda,
            ibm_sierracarboneray5_miroc6[[3]]$log_lambda,
            ibm_sierracarboneray5_miroc6[[4]]$log_lambda,
            ibm_sierracarboneray5_miroc6[[5]]$log_lambda)),
    t(rbind(ibm_sierracarboneray5_mpi_esm1_2_lr[[1]]$log_lambda,
            ibm_sierracarboneray5_mpi_esm1_2_lr[[2]]$log_lambda,
            ibm_sierracarboneray5_mpi_esm1_2_lr[[3]]$log_lambda,
            ibm_sierracarboneray5_mpi_esm1_2_lr[[4]]$log_lambda,
            ibm_sierracarboneray5_mpi_esm1_2_lr[[5]]$log_lambda)),
    t(rbind(ibm_sierracarboneray5_mri_esm2_0[[1]]$log_lambda,
            ibm_sierracarboneray5_mri_esm2_0[[2]]$log_lambda,
            ibm_sierracarboneray5_mri_esm2_0[[3]]$log_lambda,
            ibm_sierracarboneray5_mri_esm2_0[[4]]$log_lambda,
            ibm_sierracarboneray5_mri_esm2_0[[5]]$log_lambda)),
    t(rbind(ibm_sierracarboneray5_noresm2_mm[[1]]$log_lambda,
            ibm_sierracarboneray5_noresm2_mm[[2]]$log_lambda,
            ibm_sierracarboneray5_noresm2_mm[[3]]$log_lambda,
            ibm_sierracarboneray5_noresm2_mm[[4]]$log_lambda,
            ibm_sierracarboneray5_noresm2_mm[[5]]$log_lambda)))


# Add extinction

ibm_results$extinction[which(ibm_results$population == "SierraCarboneraY5")] = 
  c(rep(c(ibm_sierracarboneray5_control[[1]]$extinction,
          ibm_sierracarboneray5_control[[2]]$extinction,
          ibm_sierracarboneray5_control[[3]]$extinction,
          ibm_sierracarboneray5_control[[4]]$extinction,
          ibm_sierracarboneray5_control[[5]]$extinction), each = n_years),
    rep(c(ibm_sierracarboneray5_canesm5[[1]]$extinction,
          ibm_sierracarboneray5_canesm5[[2]]$extinction,
          ibm_sierracarboneray5_canesm5[[3]]$extinction,
          ibm_sierracarboneray5_canesm5[[4]]$extinction,
          ibm_sierracarboneray5_canesm5[[5]]$extinction), each = n_years),
    rep(c(ibm_sierracarboneray5_ec_earth3[[1]]$extinction,
          ibm_sierracarboneray5_ec_earth3[[2]]$extinction,
          ibm_sierracarboneray5_ec_earth3[[3]]$extinction,
          ibm_sierracarboneray5_ec_earth3[[4]]$extinction,
          ibm_sierracarboneray5_ec_earth3[[5]]$extinction), each = n_years),
    rep(c(ibm_sierracarboneray5_fgoals_g3[[1]]$extinction,
          ibm_sierracarboneray5_fgoals_g3[[2]]$extinction,
          ibm_sierracarboneray5_fgoals_g3[[3]]$extinction,
          ibm_sierracarboneray5_fgoals_g3[[4]]$extinction,
          ibm_sierracarboneray5_fgoals_g3[[5]]$extinction), each = n_years),
    rep(c(ibm_sierracarboneray5_gfdl_esm4[[1]]$extinction,
          ibm_sierracarboneray5_gfdl_esm4[[2]]$extinction,
          ibm_sierracarboneray5_gfdl_esm4[[3]]$extinction,
          ibm_sierracarboneray5_gfdl_esm4[[4]]$extinction,
          ibm_sierracarboneray5_gfdl_esm4[[5]]$extinction), each = n_years),
    rep(c(ibm_sierracarboneray5_giss_e2_1_g[[1]]$extinction,
          ibm_sierracarboneray5_giss_e2_1_g[[2]]$extinction,
          ibm_sierracarboneray5_giss_e2_1_g[[3]]$extinction,
          ibm_sierracarboneray5_giss_e2_1_g[[4]]$extinction,
          ibm_sierracarboneray5_giss_e2_1_g[[5]]$extinction), each = n_years),
    rep(c(ibm_sierracarboneray5_inm_cm4_8[[1]]$extinction,
          ibm_sierracarboneray5_inm_cm4_8[[2]]$extinction,
          ibm_sierracarboneray5_inm_cm4_8[[3]]$extinction,
          ibm_sierracarboneray5_inm_cm4_8[[4]]$extinction,
          ibm_sierracarboneray5_inm_cm4_8[[5]]$extinction), each = n_years),
    rep(c(ibm_sierracarboneray5_ipsl_cm6a_lr[[1]]$extinction,
          ibm_sierracarboneray5_ipsl_cm6a_lr[[2]]$extinction,
          ibm_sierracarboneray5_ipsl_cm6a_lr[[3]]$extinction,
          ibm_sierracarboneray5_ipsl_cm6a_lr[[4]]$extinction,
          ibm_sierracarboneray5_ipsl_cm6a_lr[[5]]$extinction), each = n_years),
    rep(c(ibm_sierracarboneray5_miroc6[[1]]$extinction,
          ibm_sierracarboneray5_miroc6[[2]]$extinction,
          ibm_sierracarboneray5_miroc6[[3]]$extinction,
          ibm_sierracarboneray5_miroc6[[4]]$extinction,
          ibm_sierracarboneray5_miroc6[[5]]$extinction), each = n_years),
    rep(c(ibm_sierracarboneray5_mpi_esm1_2_lr[[1]]$extinction,
          ibm_sierracarboneray5_mpi_esm1_2_lr[[2]]$extinction,
          ibm_sierracarboneray5_mpi_esm1_2_lr[[3]]$extinction,
          ibm_sierracarboneray5_mpi_esm1_2_lr[[4]]$extinction,
          ibm_sierracarboneray5_mpi_esm1_2_lr[[5]]$extinction), each = n_years),
    rep(c(ibm_sierracarboneray5_mri_esm2_0[[1]]$extinction,
          ibm_sierracarboneray5_mri_esm2_0[[2]]$extinction,
          ibm_sierracarboneray5_mri_esm2_0[[3]]$extinction,
          ibm_sierracarboneray5_mri_esm2_0[[4]]$extinction,
          ibm_sierracarboneray5_mri_esm2_0[[5]]$extinction), each = n_years),
    rep(c(ibm_sierracarboneray5_noresm2_mm[[1]]$extinction,
          ibm_sierracarboneray5_noresm2_mm[[2]]$extinction,
          ibm_sierracarboneray5_noresm2_mm[[3]]$extinction,
          ibm_sierracarboneray5_noresm2_mm[[4]]$extinction,
          ibm_sierracarboneray5_noresm2_mm[[5]]$extinction), each = n_years))


rm(list = grep("ibm_sierracarbonera", ls() , value = TRUE, invert = FALSE))



for(i in 1:(length(list.files("Output/Projections/")[grep("Natural_SierraRetin", list.files("Output/Projections/"))]))){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_SierraRetin", list.files("Output/Projections/"))][i]))
}


# Add log lambda (per simulation = per row)

ibm_results$log_lambda[which(ibm_results$population == "SierraRetinY5")] = 
  c(t(rbind(ibm_sierraretiny5_control[[1]]$log_lambda,
            ibm_sierraretiny5_control[[2]]$log_lambda,
            ibm_sierraretiny5_control[[3]]$log_lambda,
            ibm_sierraretiny5_control[[4]]$log_lambda,
            ibm_sierraretiny5_control[[5]]$log_lambda)),
    t(rbind(ibm_sierraretiny5_canesm5[[1]]$log_lambda,
            ibm_sierraretiny5_canesm5[[2]]$log_lambda,
            ibm_sierraretiny5_canesm5[[3]]$log_lambda,
            ibm_sierraretiny5_canesm5[[4]]$log_lambda,
            ibm_sierraretiny5_canesm5[[5]]$log_lambda)),
    t(rbind(ibm_sierraretiny5_ec_earth3[[1]]$log_lambda,
            ibm_sierraretiny5_ec_earth3[[2]]$log_lambda,
            ibm_sierraretiny5_ec_earth3[[3]]$log_lambda,
            ibm_sierraretiny5_ec_earth3[[4]]$log_lambda,
            ibm_sierraretiny5_ec_earth3[[5]]$log_lambda)),
    t(rbind(ibm_sierraretiny5_fgoals_g3[[1]]$log_lambda,
            ibm_sierraretiny5_fgoals_g3[[2]]$log_lambda,
            ibm_sierraretiny5_fgoals_g3[[3]]$log_lambda,
            ibm_sierraretiny5_fgoals_g3[[4]]$log_lambda,
            ibm_sierraretiny5_fgoals_g3[[5]]$log_lambda)),
    t(rbind(ibm_sierraretiny5_gfdl_esm4[[1]]$log_lambda,
            ibm_sierraretiny5_gfdl_esm4[[2]]$log_lambda,
            ibm_sierraretiny5_gfdl_esm4[[3]]$log_lambda,
            ibm_sierraretiny5_gfdl_esm4[[4]]$log_lambda,
            ibm_sierraretiny5_gfdl_esm4[[5]]$log_lambda)),
    t(rbind(ibm_sierraretiny5_giss_e2_1_g[[1]]$log_lambda,
            ibm_sierraretiny5_giss_e2_1_g[[2]]$log_lambda,
            ibm_sierraretiny5_giss_e2_1_g[[3]]$log_lambda,
            ibm_sierraretiny5_giss_e2_1_g[[4]]$log_lambda,
            ibm_sierraretiny5_giss_e2_1_g[[5]]$log_lambda)),
    t(rbind(ibm_sierraretiny5_inm_cm4_8[[1]]$log_lambda,
            ibm_sierraretiny5_inm_cm4_8[[2]]$log_lambda,
            ibm_sierraretiny5_inm_cm4_8[[3]]$log_lambda,
            ibm_sierraretiny5_inm_cm4_8[[4]]$log_lambda,
            ibm_sierraretiny5_inm_cm4_8[[5]]$log_lambda)),
    t(rbind(ibm_sierraretiny5_ipsl_cm6a_lr[[1]]$log_lambda,
            ibm_sierraretiny5_ipsl_cm6a_lr[[2]]$log_lambda,
            ibm_sierraretiny5_ipsl_cm6a_lr[[3]]$log_lambda,
            ibm_sierraretiny5_ipsl_cm6a_lr[[4]]$log_lambda,
            ibm_sierraretiny5_ipsl_cm6a_lr[[5]]$log_lambda)),
    t(rbind(ibm_sierraretiny5_miroc6[[1]]$log_lambda,
            ibm_sierraretiny5_miroc6[[2]]$log_lambda,
            ibm_sierraretiny5_miroc6[[3]]$log_lambda,
            ibm_sierraretiny5_miroc6[[4]]$log_lambda,
            ibm_sierraretiny5_miroc6[[5]]$log_lambda)),
    t(rbind(ibm_sierraretiny5_mpi_esm1_2_lr[[1]]$log_lambda,
            ibm_sierraretiny5_mpi_esm1_2_lr[[2]]$log_lambda,
            ibm_sierraretiny5_mpi_esm1_2_lr[[3]]$log_lambda,
            ibm_sierraretiny5_mpi_esm1_2_lr[[4]]$log_lambda,
            ibm_sierraretiny5_mpi_esm1_2_lr[[5]]$log_lambda)),
    t(rbind(ibm_sierraretiny5_mri_esm2_0[[1]]$log_lambda,
            ibm_sierraretiny5_mri_esm2_0[[2]]$log_lambda,
            ibm_sierraretiny5_mri_esm2_0[[3]]$log_lambda,
            ibm_sierraretiny5_mri_esm2_0[[4]]$log_lambda,
            ibm_sierraretiny5_mri_esm2_0[[5]]$log_lambda)),
    t(rbind(ibm_sierraretiny5_noresm2_mm[[1]]$log_lambda,
            ibm_sierraretiny5_noresm2_mm[[2]]$log_lambda,
            ibm_sierraretiny5_noresm2_mm[[3]]$log_lambda,
            ibm_sierraretiny5_noresm2_mm[[4]]$log_lambda,
            ibm_sierraretiny5_noresm2_mm[[5]]$log_lambda)))


# Add extinction

ibm_results$extinction[which(ibm_results$population == "SierraRetinY5")] = 
  c(rep(c(ibm_sierraretiny5_control[[1]]$extinction,
          ibm_sierraretiny5_control[[2]]$extinction,
          ibm_sierraretiny5_control[[3]]$extinction,
          ibm_sierraretiny5_control[[4]]$extinction,
          ibm_sierraretiny5_control[[5]]$extinction), each = n_years),
    rep(c(ibm_sierraretiny5_canesm5[[1]]$extinction,
          ibm_sierraretiny5_canesm5[[2]]$extinction,
          ibm_sierraretiny5_canesm5[[3]]$extinction,
          ibm_sierraretiny5_canesm5[[4]]$extinction,
          ibm_sierraretiny5_canesm5[[5]]$extinction), each = n_years),
    rep(c(ibm_sierraretiny5_ec_earth3[[1]]$extinction,
          ibm_sierraretiny5_ec_earth3[[2]]$extinction,
          ibm_sierraretiny5_ec_earth3[[3]]$extinction,
          ibm_sierraretiny5_ec_earth3[[4]]$extinction,
          ibm_sierraretiny5_ec_earth3[[5]]$extinction), each = n_years),
    rep(c(ibm_sierraretiny5_fgoals_g3[[1]]$extinction,
          ibm_sierraretiny5_fgoals_g3[[2]]$extinction,
          ibm_sierraretiny5_fgoals_g3[[3]]$extinction,
          ibm_sierraretiny5_fgoals_g3[[4]]$extinction,
          ibm_sierraretiny5_fgoals_g3[[5]]$extinction), each = n_years),
    rep(c(ibm_sierraretiny5_gfdl_esm4[[1]]$extinction,
          ibm_sierraretiny5_gfdl_esm4[[2]]$extinction,
          ibm_sierraretiny5_gfdl_esm4[[3]]$extinction,
          ibm_sierraretiny5_gfdl_esm4[[4]]$extinction,
          ibm_sierraretiny5_gfdl_esm4[[5]]$extinction), each = n_years),
    rep(c(ibm_sierraretiny5_giss_e2_1_g[[1]]$extinction,
          ibm_sierraretiny5_giss_e2_1_g[[2]]$extinction,
          ibm_sierraretiny5_giss_e2_1_g[[3]]$extinction,
          ibm_sierraretiny5_giss_e2_1_g[[4]]$extinction,
          ibm_sierraretiny5_giss_e2_1_g[[5]]$extinction), each = n_years),
    rep(c(ibm_sierraretiny5_inm_cm4_8[[1]]$extinction,
          ibm_sierraretiny5_inm_cm4_8[[2]]$extinction,
          ibm_sierraretiny5_inm_cm4_8[[3]]$extinction,
          ibm_sierraretiny5_inm_cm4_8[[4]]$extinction,
          ibm_sierraretiny5_inm_cm4_8[[5]]$extinction), each = n_years),
    rep(c(ibm_sierraretiny5_ipsl_cm6a_lr[[1]]$extinction,
          ibm_sierraretiny5_ipsl_cm6a_lr[[2]]$extinction,
          ibm_sierraretiny5_ipsl_cm6a_lr[[3]]$extinction,
          ibm_sierraretiny5_ipsl_cm6a_lr[[4]]$extinction,
          ibm_sierraretiny5_ipsl_cm6a_lr[[5]]$extinction), each = n_years),
    rep(c(ibm_sierraretiny5_miroc6[[1]]$extinction,
          ibm_sierraretiny5_miroc6[[2]]$extinction,
          ibm_sierraretiny5_miroc6[[3]]$extinction,
          ibm_sierraretiny5_miroc6[[4]]$extinction,
          ibm_sierraretiny5_miroc6[[5]]$extinction), each = n_years),
    rep(c(ibm_sierraretiny5_mpi_esm1_2_lr[[1]]$extinction,
          ibm_sierraretiny5_mpi_esm1_2_lr[[2]]$extinction,
          ibm_sierraretiny5_mpi_esm1_2_lr[[3]]$extinction,
          ibm_sierraretiny5_mpi_esm1_2_lr[[4]]$extinction,
          ibm_sierraretiny5_mpi_esm1_2_lr[[5]]$extinction), each = n_years),
    rep(c(ibm_sierraretiny5_mri_esm2_0[[1]]$extinction,
          ibm_sierraretiny5_mri_esm2_0[[2]]$extinction,
          ibm_sierraretiny5_mri_esm2_0[[3]]$extinction,
          ibm_sierraretiny5_mri_esm2_0[[4]]$extinction,
          ibm_sierraretiny5_mri_esm2_0[[5]]$extinction), each = n_years),
    rep(c(ibm_sierraretiny5_noresm2_mm[[1]]$extinction,
          ibm_sierraretiny5_noresm2_mm[[2]]$extinction,
          ibm_sierraretiny5_noresm2_mm[[3]]$extinction,
          ibm_sierraretiny5_noresm2_mm[[4]]$extinction,
          ibm_sierraretiny5_noresm2_mm[[5]]$extinction), each = n_years))


rm(list = grep("ibm_sierraretin", ls() , value = TRUE, invert = FALSE))


for(i in 1:(length(list.files("Output/Projections/")[grep("Natural_Vertedero", list.files("Output/Projections/"))]))){
  
  print(i)
  
  load(paste0("Output/Projections/", list.files("Output/Projections/")[grep("Natural_Vertedero", list.files("Output/Projections/"))][i]))
}


# Add log lambda (per simulation = per row)

ibm_results$log_lambda[which(ibm_results$population == "Vertedero")] = 
  c(t(rbind(ibm_vertedero_control[[1]]$log_lambda,
            ibm_vertedero_control[[2]]$log_lambda,
            ibm_vertedero_control[[3]]$log_lambda,
            ibm_vertedero_control[[4]]$log_lambda,
            ibm_vertedero_control[[5]]$log_lambda)),
    t(rbind(ibm_vertedero_canesm5[[1]]$log_lambda,
            ibm_vertedero_canesm5[[2]]$log_lambda,
            ibm_vertedero_canesm5[[3]]$log_lambda,
            ibm_vertedero_canesm5[[4]]$log_lambda,
            ibm_vertedero_canesm5[[5]]$log_lambda)),
    t(rbind(ibm_vertedero_ec_earth3[[1]]$log_lambda,
            ibm_vertedero_ec_earth3[[2]]$log_lambda,
            ibm_vertedero_ec_earth3[[3]]$log_lambda,
            ibm_vertedero_ec_earth3[[4]]$log_lambda,
            ibm_vertedero_ec_earth3[[5]]$log_lambda)),
    t(rbind(ibm_vertedero_fgoals_g3[[1]]$log_lambda,
            ibm_vertedero_fgoals_g3[[2]]$log_lambda,
            ibm_vertedero_fgoals_g3[[3]]$log_lambda,
            ibm_vertedero_fgoals_g3[[4]]$log_lambda,
            ibm_vertedero_fgoals_g3[[5]]$log_lambda)),
    t(rbind(ibm_vertedero_gfdl_esm4[[1]]$log_lambda,
            ibm_vertedero_gfdl_esm4[[2]]$log_lambda,
            ibm_vertedero_gfdl_esm4[[3]]$log_lambda,
            ibm_vertedero_gfdl_esm4[[4]]$log_lambda,
            ibm_vertedero_gfdl_esm4[[5]]$log_lambda)),
    t(rbind(ibm_vertedero_giss_e2_1_g[[1]]$log_lambda,
            ibm_vertedero_giss_e2_1_g[[2]]$log_lambda,
            ibm_vertedero_giss_e2_1_g[[3]]$log_lambda,
            ibm_vertedero_giss_e2_1_g[[4]]$log_lambda,
            ibm_vertedero_giss_e2_1_g[[5]]$log_lambda)),
    t(rbind(ibm_vertedero_inm_cm4_8[[1]]$log_lambda,
            ibm_vertedero_inm_cm4_8[[2]]$log_lambda,
            ibm_vertedero_inm_cm4_8[[3]]$log_lambda,
            ibm_vertedero_inm_cm4_8[[4]]$log_lambda,
            ibm_vertedero_inm_cm4_8[[5]]$log_lambda)),
    t(rbind(ibm_vertedero_ipsl_cm6a_lr[[1]]$log_lambda,
            ibm_vertedero_ipsl_cm6a_lr[[2]]$log_lambda,
            ibm_vertedero_ipsl_cm6a_lr[[3]]$log_lambda,
            ibm_vertedero_ipsl_cm6a_lr[[4]]$log_lambda,
            ibm_vertedero_ipsl_cm6a_lr[[5]]$log_lambda)),
    t(rbind(ibm_vertedero_miroc6[[1]]$log_lambda,
            ibm_vertedero_miroc6[[2]]$log_lambda,
            ibm_vertedero_miroc6[[3]]$log_lambda,
            ibm_vertedero_miroc6[[4]]$log_lambda,
            ibm_vertedero_miroc6[[5]]$log_lambda)),
    t(rbind(ibm_vertedero_mpi_esm1_2_lr[[1]]$log_lambda,
            ibm_vertedero_mpi_esm1_2_lr[[2]]$log_lambda,
            ibm_vertedero_mpi_esm1_2_lr[[3]]$log_lambda,
            ibm_vertedero_mpi_esm1_2_lr[[4]]$log_lambda,
            ibm_vertedero_mpi_esm1_2_lr[[5]]$log_lambda)),
    t(rbind(ibm_vertedero_mri_esm2_0[[1]]$log_lambda,
            ibm_vertedero_mri_esm2_0[[2]]$log_lambda,
            ibm_vertedero_mri_esm2_0[[3]]$log_lambda,
            ibm_vertedero_mri_esm2_0[[4]]$log_lambda,
            ibm_vertedero_mri_esm2_0[[5]]$log_lambda)),
    t(rbind(ibm_vertedero_noresm2_mm[[1]]$log_lambda,
            ibm_vertedero_noresm2_mm[[2]]$log_lambda,
            ibm_vertedero_noresm2_mm[[3]]$log_lambda,
            ibm_vertedero_noresm2_mm[[4]]$log_lambda,
            ibm_vertedero_noresm2_mm[[5]]$log_lambda)))


# Add extinction

ibm_results$extinction[which(ibm_results$population == "Vertedero")] = 
  c(rep(c(ibm_vertedero_control[[1]]$extinction,
          ibm_vertedero_control[[2]]$extinction,
          ibm_vertedero_control[[3]]$extinction,
          ibm_vertedero_control[[4]]$extinction,
          ibm_vertedero_control[[5]]$extinction), each = n_years),
    rep(c(ibm_vertedero_canesm5[[1]]$extinction,
          ibm_vertedero_canesm5[[2]]$extinction,
          ibm_vertedero_canesm5[[3]]$extinction,
          ibm_vertedero_canesm5[[4]]$extinction,
          ibm_vertedero_canesm5[[5]]$extinction), each = n_years),
    rep(c(ibm_vertedero_ec_earth3[[1]]$extinction,
          ibm_vertedero_ec_earth3[[2]]$extinction,
          ibm_vertedero_ec_earth3[[3]]$extinction,
          ibm_vertedero_ec_earth3[[4]]$extinction,
          ibm_vertedero_ec_earth3[[5]]$extinction), each = n_years),
    rep(c(ibm_vertedero_fgoals_g3[[1]]$extinction,
          ibm_vertedero_fgoals_g3[[2]]$extinction,
          ibm_vertedero_fgoals_g3[[3]]$extinction,
          ibm_vertedero_fgoals_g3[[4]]$extinction,
          ibm_vertedero_fgoals_g3[[5]]$extinction), each = n_years),
    rep(c(ibm_vertedero_gfdl_esm4[[1]]$extinction,
          ibm_vertedero_gfdl_esm4[[2]]$extinction,
          ibm_vertedero_gfdl_esm4[[3]]$extinction,
          ibm_vertedero_gfdl_esm4[[4]]$extinction,
          ibm_vertedero_gfdl_esm4[[5]]$extinction), each = n_years),
    rep(c(ibm_vertedero_giss_e2_1_g[[1]]$extinction,
          ibm_vertedero_giss_e2_1_g[[2]]$extinction,
          ibm_vertedero_giss_e2_1_g[[3]]$extinction,
          ibm_vertedero_giss_e2_1_g[[4]]$extinction,
          ibm_vertedero_giss_e2_1_g[[5]]$extinction), each = n_years),
    rep(c(ibm_vertedero_inm_cm4_8[[1]]$extinction,
          ibm_vertedero_inm_cm4_8[[2]]$extinction,
          ibm_vertedero_inm_cm4_8[[3]]$extinction,
          ibm_vertedero_inm_cm4_8[[4]]$extinction,
          ibm_vertedero_inm_cm4_8[[5]]$extinction), each = n_years),
    rep(c(ibm_vertedero_ipsl_cm6a_lr[[1]]$extinction,
          ibm_vertedero_ipsl_cm6a_lr[[2]]$extinction,
          ibm_vertedero_ipsl_cm6a_lr[[3]]$extinction,
          ibm_vertedero_ipsl_cm6a_lr[[4]]$extinction,
          ibm_vertedero_ipsl_cm6a_lr[[5]]$extinction), each = n_years),
    rep(c(ibm_vertedero_miroc6[[1]]$extinction,
          ibm_vertedero_miroc6[[2]]$extinction,
          ibm_vertedero_miroc6[[3]]$extinction,
          ibm_vertedero_miroc6[[4]]$extinction,
          ibm_vertedero_miroc6[[5]]$extinction), each = n_years),
    rep(c(ibm_vertedero_mpi_esm1_2_lr[[1]]$extinction,
          ibm_vertedero_mpi_esm1_2_lr[[2]]$extinction,
          ibm_vertedero_mpi_esm1_2_lr[[3]]$extinction,
          ibm_vertedero_mpi_esm1_2_lr[[4]]$extinction,
          ibm_vertedero_mpi_esm1_2_lr[[5]]$extinction), each = n_years),
    rep(c(ibm_vertedero_mri_esm2_0[[1]]$extinction,
          ibm_vertedero_mri_esm2_0[[2]]$extinction,
          ibm_vertedero_mri_esm2_0[[3]]$extinction,
          ibm_vertedero_mri_esm2_0[[4]]$extinction,
          ibm_vertedero_mri_esm2_0[[5]]$extinction), each = n_years),
    rep(c(ibm_vertedero_noresm2_mm[[1]]$extinction,
          ibm_vertedero_noresm2_mm[[2]]$extinction,
          ibm_vertedero_noresm2_mm[[3]]$extinction,
          ibm_vertedero_noresm2_mm[[4]]$extinction,
          ibm_vertedero_noresm2_mm[[5]]$extinction), each = n_years))


rm(list = grep("ibm_vertedero", ls() , value = TRUE, invert = FALSE))

write.csv(ibm_results, file = "Output/Projections/IBM_Results.csv", row.names = F)


## 2.2. Summarising results ----
# -------------------------

ibm_results$habitat = factor(ibm_results$habitat, 
                              levels = c("Natural", "Anthropogenic"))
ibm_results$climate = factor(ibm_results$climate,
                             levels = c("Control", "Climate change"))

mean_lambda = 
  aggregate(log_lambda ~ habitat + climate + population, 
            data = ibm_results[which(!is.infinite(ibm_results$log_lambda)), ],
            FUN = function(x) quantile(x, probs = c(0.025, 0.975), na.rm = T))

mean_lambda = rbind(mean_lambda[which(mean_lambda$habitat == "Natural"), ],
                    mean_lambda[which(mean_lambda$habitat == "Anthropogenic"), ])

ibm_results$population_full = ibm_results$population

ibm_results$population_full[which(ibm_results$population == "SierraCarboneraY5")] = "Sierra\nCarbonera\nYoung"
ibm_results$population_full[which(ibm_results$population == "SierraRetinY5")] = "Sierra del\nRetín\nYoung"
ibm_results$population_full[which(ibm_results$population == "MonteraTorero")] = "Montera\ndel\nTorero"
ibm_results$population_full[which(ibm_results$population == "Retin")] = "Sierra del\nRetín\nDisturbed"
ibm_results$population_full[which(ibm_results$population == "Prisoneros")] = "Prisioneros"
ibm_results$population_full[which(ibm_results$population == "SCarbDist")] = "Sierra\nCarbonera\nDisturbed"


mean_extinction = 
  aggregate(extinction ~ habitat + climate + population, 
            data = ibm_results,
            FUN = mean, na.rm = T)

mean_extinction_avg = 
  aggregate(extinction ~ habitat + climate, 
            data = mean_extinction,
            FUN = function(x) c(mean(x, na.rm = T), quantile(x, probs = c(0.025, 0.975), na.rm = T)))

mean_extinction_avg = cbind(mean_extinction_avg[, c("habitat", "climate")],
                            mean_extinction_avg$extinction)
colnames(mean_extinction_avg) = c("habitat", "climate", 
                                  "mean", "lwr", "upr")




###########################################################################
#
# 3. Plotting results ----
#
###########################################################################

cbbPalette <- c("#111111", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## 3.1. Log lambda ----
# ----------------

png(filename = "Output/Plots/LogLambda.png", 
    width = 16,
    height = 10,
    units = "cm",
    bg = "white",
    res = 600)

ggplot(ibm_results[-which(is.na(ibm_results$log_lambda)), ], aes(x = population_full, y = log_lambda, 
                                                                 fill = climate, colour = climate)) +
  facet_wrap(~ habitat, scales = "free") +
  geom_boxplot(ymin = mean_lambda$log_lambda[, 1],
               ymax = mean_lambda$log_lambda[, 2], outlier.size = 0.2, size = 0.2, alpha = 0.5) +
  stat_summary(fun = mean, colour = cbbPalette[5], aes(shape = "Mean", group = climate), geom = "point", shape = 17, size = 1, position = position_dodge(width = 0.75)) +
  labs(x = "Population", 
       y = "Stochastic log \u03BB") +
  scale_fill_manual(values = cbbPalette[1:2],
                    name = "Climate") +
  scale_colour_manual(values = cbbPalette[1:2],
                      name = "Climate") +
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


lambda_agg = aggregate(log_lambda ~ climate + habitat + year + simulation, ibm_results[-which(is.nan(ibm_results$log_lambda) | is.infinite(ibm_results$log_lambda)), ], mean)

mean_lambda_agg = 
  aggregate(log_lambda ~ climate + habitat, 
            data = lambda_agg,
            FUN = function(x) c(mean(x), quantile(x, probs = c(0.025, 0.975), na.rm = T)))


png(filename = "Output/Plots/LogLambda_Average.png", 
    width = 14,
    height = 10,
    units = "cm",
    bg = "white",
    res = 600)

lambda_average_plot = ggplot(lambda_agg, aes(x = habitat, y = log_lambda, 
                       fill = climate, colour = climate)) +
  geom_boxplot(ymin = mean_lambda_agg$log_lambda[, 2],
               ymax = mean_lambda_agg$log_lambda[, 3], outlier.size = 0.2, size = 0.2, alpha = 0.5) +
  stat_summary(fun = mean, colour = cbbPalette[5], aes(shape = "Mean", group = climate), geom = "point", shape = 17, size = 1, position = position_dodge(width = 0.75)) +
  labs(x = "Habitat type", 
       y = "Stochastic log \u03BB") +
  scale_fill_manual(values = cbbPalette[1:2],
                    name = "Climate scenario") +
  scale_colour_manual(values = cbbPalette[1:2],
                      name = "Climate scenario") +
  scale_shape_manual("", values = c("Mean" = 17)) +
  theme_minimal() %+replace%
  theme(axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title.x = element_text(size = 8, colour = "black", 
                                    margin = margin(t = 0, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 8, colour = "black", angle = 90, 
                                    margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 2, r = 0, b = 5, l = 0)), 
        axis.text.y = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        panel.grid = element_blank(),
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 8), 
        legend.position = "right", 
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(1, "lines"),
        strip.background = element_blank())

dev.off()


## 3.2. Extinction ----
# ----------------

png(filename = "Output/Plots/Extinction.png", 
    width = 14,
    height = 10,
    units = "cm",
    bg = "white",
    res = 600)

ggplot(mean_extinction, aes(x = population, y = extinction, 
                            fill = habitat, colour = climate)) +
  geom_point(shape = 21, size = 4, position = position_dodge(width = 0.75)) +
  labs(x = "Population", 
       y = "Extinction probability") +
  scale_fill_manual(values = cbbPalette[c(6, 4)],
                    name = "Habitat type") +
  scale_colour_manual(values = cbbPalette[1:2],
                      name = "Climate") +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title.x = element_text(size = 7, colour = "black", 
                                    margin = margin(t = 0, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 7, colour = "black", 
                                    margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust = 1, 
                                   margin = margin(t = 2, r = 0, b = 5, l = 0)), 
        axis.text.y = element_text(size = 6, colour = "black", 
                                   margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(face = "bold", size = 7), 
        legend.position = "right", 
        legend.key.size = unit(2, "lines"),
        strip.background = element_blank())

dev.off()


png(filename = "Output/Plots/Extinction_Average.png", 
    width = 14,
    height = 10,
    units = "cm",
    bg = "white",
    res = 600)

extinction_average_plot = ggplot(mean_extinction_avg, aes(x = habitat, y = mean, 
                                colour = climate)) +
  geom_point(shape = 20, size = 1.5, position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                width = 0.2,
                position = position_dodge(width = 0.75),
                linewidth = 0.4) + 
  labs(x = "Habitat type", 
       y = "Extinction probability") +
  scale_colour_manual(values = cbbPalette[1:2],
                      name = "Climate scenario") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title.x = element_text(size = 8, colour = "black", 
                                    margin = margin(t = 0, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 8, colour = "black", 
                                    margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 2, r = 0, b = 5, l = 0)), 
        axis.text.y = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 0, r = 2, b = 0, l = 0)),
        panel.grid = element_blank(),
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 8), 
        legend.position = "right", 
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(1, "lines"),
        strip.background = element_blank())

dev.off()


png(filename = "Output/Plots/Figure4.png", 
    width = 8, 
    height = 10, 
    units = "cm", 
    bg = "white", 
    res = 600, 
    type = "cairo")

lambda_average_plot + extinction_average_plot +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = 'a',
                  tag_prefix = "(",
                  tag_suffix = ")")  &
  theme(plot.tag = element_text(size = 10))

dev.off()