############################################################################
#
# This script projects the dewy-pine populations using an individual-based
# model (IBM). The goal is to assess the effect of climate change on 
# populations under natural fire regimes where seed dormancy follows a
# natural pattern and populations that do not burn anymore and have
# partly lost dormancy.
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

library(lubridate)
library(mgcv)
library(crch)
library(snowfall)


## 1.3. Loading data ----
# ------------------

# Dewy-pine data
droso_anthropogenic = read.csv("Data/droso_anthropogenic.csv")

droso_seedbank = read.csv("Data/droso_seedbank_anthropogenic.csv")

seeds_per_flower = 9.8


# Correction factors
seedSurv = 0.33


# Vital-rate models
load("Output/Models/Survival_GAM_Anthropogenic.RData")
load("Output/Models/Growth_GAM_Anthropogenic.RData")
load("Output/Models/FloweringProb_GAM_Anthropogenic.RData")
load("Output/Models/NbFlowers_GAM_Anthropogenic.RData")
load("Output/Models/SeedlingSize_GAM_Anthropogenic.RData")


# Average density per square for covariate standardization
nbSquares = aggregate(quadratID ~ time + site, 
                      data = droso_anthropogenic, 
                      function(x) length(unique(x)))
density_per_square = aggregate(abLarge_unscaled ~ quadratID + time + site, 
                               data = droso_anthropogenic, 
                               function(x) unique(x))
yearly_density_per_square = aggregate(abLarge_unscaled ~ time + site, 
                                      data = density_per_square, 
                                      function(x) sum(x))
yearly_density_per_square$abLarge_unscaled = yearly_density_per_square$abLarge_unscaled/nbSquares$quadratID


# Year- and population-specific climatic variables for covariate standardization
summerT_timeSeries = aggregate(summerT_unscaled ~ time + site, 
                               data = droso_anthropogenic, mean)
prevwinterT_timeSeries = aggregate(prevwinterT_unscaled ~ time + site, 
                                   data = droso_anthropogenic, mean)
fallR_timeSeries = aggregate(fallR_unscaled ~ time + site, 
                             data = droso_anthropogenic, mean)
prevfallR_timeSeries = aggregate(prevfallR_unscaled ~ time + site, 
                                 data = droso_anthropogenic, mean)
prevwinterR_timeSeries = aggregate(prevwinterR_unscaled ~ time + site, 
                                   data = droso_anthropogenic, mean)



###########################################################################
#
# 2. Building vital-rate functions ----
#
###########################################################################

## 2.1. Survival ----
# --------------

survival_function = function(size, abLarge, fallR, summerT, year, population){
  
  # Standardize covariates
  size_scaled = (size - mean(droso_anthropogenic$size_unscaled, na.rm = T)) / (2 * sd(droso_anthropogenic$size_unscaled, na.rm = T))
  abLarge_scaled = (abLarge - mean(yearly_density_per_square$abLarge_unscaled, na.rm = T)) / (2 * sd(yearly_density_per_square$abLarge_unscaled, na.rm = T))
  fallR_scaled = (fallR - mean(fallR_timeSeries$fallR_unscaled, na.rm = T)) / (2 * sd(fallR_timeSeries$fallR_unscaled, na.rm = T))
  summerT_scaled = (summerT - mean(summerT_timeSeries$summerT_unscaled, na.rm = T)) / (2 * sd(summerT_timeSeries$summerT_unscaled, na.rm = T))
  
  # Calculate survival
  survival = predict(surv_anthropogenic, 
                     newdata = data.frame(size = size_scaled,
                                          abLarge = abLarge_scaled,
                                          fallR = fallR_scaled,
                                          summerT = summerT_scaled,
                                          time = year,
                                          site = population), type = "response")
  
  return(survival)
}


## 2.2. Growth ----
# ------------

growth_function = function(size, abLarge, summerT, year, population){
  
  # Standardize covariates
  size_scaled = (size - mean(droso_anthropogenic$size_unscaled, na.rm = T)) / (2 * sd(droso_anthropogenic$size_unscaled, na.rm = T))
  abLarge_scaled = (abLarge - mean(yearly_density_per_square$abLarge_unscaled, na.rm = T)) / (2 * sd(yearly_density_per_square$abLarge_unscaled, na.rm = T))
  summerT_scaled = (summerT - mean(summerT_timeSeries$summerT_unscaled, na.rm = T)) / (2 * sd(summerT_timeSeries$summerT_unscaled, na.rm = T))
  
  # Get parameters from growth model (mean, sd, and degrees of freedom)
  growth_mean = predict(growth_anthropogenic, 
                        newdata = data.frame(size = size_scaled,
                                             abLarge = abLarge_scaled,
                                             summerT = summerT_scaled,
                                             time = year,
                                             site = population), type = "response")
  
  growth_sd = family(growth_anthropogenic)$getTheta(trans = T)[2] # Using trans = T, the second value of theta is sigma, the standard deviation (https://stats.stackexchange.com/questions/550339/extracting-the-degrees-of-freedom-of-t-distribution-of-a-gam)
  
  growth_df = family(growth_anthropogenic)$getTheta(trans = T)[1] # Using trans = T, the first value of theta is nu, the degrees of freedom (https://stats.stackexchange.com/questions/550339/extracting-the-degrees-of-freedom-of-t-distribution-of-a-gam)
  
  return(list(mean = growth_mean, 
              sd = growth_sd, 
              df = growth_df))
  
}


## 2.3. Seedling sizes ----
# --------------------

seedling_size_function = function(abLarge, prevwinterT, year, population){
  
  # Standardize covariates
  abLarge_scaled = (abLarge - mean(yearly_density_per_square$abLarge_unscaled, na.rm = T)) / (2 * sd(yearly_density_per_square$abLarge_unscaled, na.rm = T))
  prevwinterT_scaled = (prevwinterT - mean(prevwinterT_timeSeries$prevwinterT_unscaled, na.rm = T)) / (2 * sd(prevwinterT_timeSeries$prevwinterT_unscaled, na.rm = T))
  
  # Get parameters from seedling size model (mean, sd, and degrees of freedom)
  seedling_size_mean = as.numeric(predict(seedlingSize_anthropogenic, 
                                          newdata = data.frame(abLarge = abLarge_scaled,
                                                               prevwinterT = prevwinterT_scaled,
                                                               time = year,
                                                               site = population), type = "response"))
  
  seedling_size_sd = family(seedlingSize_anthropogenic)$getTheta(trans = T)[2] # Using trans = T, the second value of theta is sigma, the standard deviation (https://stats.stackexchange.com/questions/550339/extracting-the-degrees-of-freedom-of-t-distribution-of-a-gam)
  
  seedling_size_df = family(seedlingSize_anthropogenic)$getTheta(trans = T)[1] # Using trans = T, the first value of theta is nu, the degrees of freedom (https://stats.stackexchange.com/questions/550339/extracting-the-degrees-of-freedom-of-t-distribution-of-a-gam)
  
  return(list(mean = seedling_size_mean, 
              sd = seedling_size_sd, 
              df = seedling_size_df))
  
}


## 2.4. Flowering probability ----
# ---------------------------

flowering_function = function(size, abLarge, prevwinterR, year, population){
  
  # Standardize covariates
  size_scaled = (size - mean(droso_anthropogenic$size_unscaled, na.rm = T)) / (2 * sd(droso_anthropogenic$size_unscaled, na.rm = T))
  abLarge_scaled = (abLarge - mean(yearly_density_per_square$abLarge_unscaled, na.rm = T)) / (2 * sd(yearly_density_per_square$abLarge_unscaled, na.rm = T))
  prevwinterR_scaled = (prevwinterR - mean(prevwinterR_timeSeries$prevwinterR_unscaled, na.rm = T)) / (2 * sd(prevwinterR_timeSeries$prevwinterR_unscaled, na.rm = T))
  
  # Calculate flowering probability
  flowering = predict(flowering_anthropogenic, 
                      newdata = data.frame(size = size_scaled,
                                           abLarge = abLarge_scaled,
                                           prevwinterR = prevwinterR_scaled,
                                           time = year,
                                           site = population), type = "response")
  
  return(flowering)
}


## 2.5. Number of flowers ----
# -----------------------

nbFlowers_function = function(size, prevfallR, year, population){
  
  # Standardize covariates
  size_scaled = (size - mean(droso_anthropogenic$size_unscaled, na.rm = T)) / (2 * sd(droso_anthropogenic$size_unscaled, na.rm = T))
  prevfallR_scaled = (size - mean(prevfallR_timeSeries$prevfallR_unscaled, na.rm = T)) / (2 * sd(prevfallR_timeSeries$prevfallR_unscaled, na.rm = T))
  
  # Calculate number of flowers
  nbFlowers = predict(nbFlow_anthropogenic, 
                      newdata = data.frame(size = size_scaled,
                                           prevfallR = prevfallR_scaled,
                                           time = year,
                                           site = population), type = "response")
  
  return(nbFlowers)
}


## 2.6. Immediate germination (goCont) ----
# ------------------------------------

goCont_function = function(population){
  
  goCont = droso_seedbank$value[which(droso_seedbank$vital_rate == "goCont" &
                                        droso_seedbank$population == population)]
  
  if(population == "SCarbDist") goCont = goCont * 0.4
  
  return(goCont)
}


## 2.7. Staying in the seed bank (staySB) ----
# ---------------------------------------

staySB_function = function(population){
  
  staySB = droso_seedbank$value[which(droso_seedbank$vital_rate == "staySB" &
                                        droso_seedbank$population == population)]
  
  return(staySB)
}


## 2.8. Germinating out of the seed bank (outSB) ----
# ----------------------------------------------

outSB_function = function(population, seedSurv, seedlingSurv){
  
  outSB = droso_seedbank$value[which(droso_seedbank$vital_rate == "outSB" &
                                       droso_seedbank$population == population)] *
    seedSurv
  
  if(population == "SCarbDist") outSB = outSB * 0.4
  
  return(outSB)
}




###########################################################################
#
# 3. Climatic projection data ----
#
###########################################################################

CanESM5_climate = read.csv("Data/CanESM5_MonthlyClimateProjection.csv")
EC_Earth3_climate = read.csv("Data/EC_Earth3_MonthlyClimateProjection.csv")
FGOALS_G3_climate = read.csv("Data/FGOALS_G3_MonthlyClimateProjection.csv")
GFDL_ESM4_climate = read.csv("Data/GFDL_ESM4_MonthlyClimateProjection.csv")
GISS_E2_1_G_climate = read.csv("Data/GISS_E2_1_G_MonthlyClimateProjection.csv")
INM_CM4_8_climate = read.csv("Data/INM_CM4_8_MonthlyClimateProjection.csv")
IPSL_CM6A_LR_climate = read.csv("Data/IPSL_CM6A_LR_MonthlyClimateProjection.csv")
MIROC6_climate = read.csv("Data/MIROC6_MonthlyClimateProjection.csv")
MPI_ESM1_2_LR_climate = read.csv("Data/MPI_ESM1_2_LR_MonthlyClimateProjection.csv")
MRI_ESM2_0_climate = read.csv("Data/MRI_ESM2_0_MonthlyClimateProjection.csv")
NorESM2_MM_climate = read.csv("Data/NorESM2_MM_MonthlyClimateProjection.csv")


CanESM5_climate$date = as.Date(paste(CanESM5_climate$year, 
                                     CanESM5_climate$month, "1", 
                                     sep = "-"))
EC_Earth3_climate$date = as.Date(paste(EC_Earth3_climate$year, 
                                       EC_Earth3_climate$month, "1", 
                                       sep = "-"))
FGOALS_G3_climate$date = as.Date(paste(FGOALS_G3_climate$year, 
                                       FGOALS_G3_climate$month, "1", 
                                       sep = "-"))
GFDL_ESM4_climate$date = as.Date(paste(GFDL_ESM4_climate$year, 
                                       GFDL_ESM4_climate$month, "1", 
                                       sep = "-"))
GISS_E2_1_G_climate$date = as.Date(paste(GISS_E2_1_G_climate$year,
                                         GISS_E2_1_G_climate$month, "1",
                                         sep = "-"))
INM_CM4_8_climate$date = as.Date(paste(INM_CM4_8_climate$year,
                                       INM_CM4_8_climate$month, "1", 
                                       sep = "-"))
IPSL_CM6A_LR_climate$date = as.Date(paste(IPSL_CM6A_LR_climate$year, 
                                          IPSL_CM6A_LR_climate$month, "1", 
                                          sep = "-"))
MIROC6_climate$date = as.Date(paste(MIROC6_climate$year, 
                                    MIROC6_climate$month, "1", 
                                    sep = "-"))
MPI_ESM1_2_LR_climate$date = as.Date(paste(MPI_ESM1_2_LR_climate$year, 
                                           MPI_ESM1_2_LR_climate$month, "1", 
                                           sep = "-"))
MRI_ESM2_0_climate$date = as.Date(paste(MRI_ESM2_0_climate$year, 
                                        MRI_ESM2_0_climate$month, "1", 
                                        sep = "-"))
NorESM2_MM_climate$date = as.Date(paste(NorESM2_MM_climate$year, 
                                        NorESM2_MM_climate$month, "1", 
                                        sep = "-"))


proj_clim = expand.grid(year = seq(2021, 2100),
                        population = unique(droso_anthropogenic$site),
                        model = c("CanESM5", "EC_Earth3", "FGOALS_G3", "GFDL_ESM4",
                                  "GISS_E2_1_G", "INM_CM4_8", "IPSL_CM6A_LR", "MIROC6",
                                  "MPI_ESM1_2_LR", "MRI_ESM2_0", "NorESM2_MM"))
proj_clim$date = as.Date(paste(proj_clim$year, "5", "1", sep = "-"))
proj_clim$summerT_unscaled = NA
proj_clim$prevwinterT_unscaled = NA
proj_clim$fallR_unscaled = NA
proj_clim$prevfallR_unscaled = NA
proj_clim$prevwinterR_unscaled = NA

climate_cov = c("summerT", "prevwinterT", "fallR", "prevfallR", "prevwinterR")

for(cov in climate_cov){
  
  for(r in 1:nrow(proj_clim)){
    
    model_data = get(paste0(proj_clim$model[r], "_climate"))
    
    if(cov == "summerT"){
      
      start_date = (proj_clim$date[r]) %m+% months(0)
      end_date = (proj_clim$date[r]) %m+% months(5)
      
      proj_clim$summerT_unscaled[r] = mean(model_data$mean_max_temp[model_data$date >= start_date & 
                                                                      model_data$date < end_date & 
                                                                      model_data$pop %in% proj_clim$population[r]])
      
    }
    
    else if(cov == "prevwinterT"){
      
      start_date = (proj_clim$date[r]) %m+% months(-4)
      end_date = (proj_clim$date[r]) %m+% months(0)
      
      proj_clim$prevwinterT_unscaled[r] = mean(model_data$mean_max_temp[model_data$date >= start_date & 
                                                                          model_data$date < end_date & 
                                                                          model_data$pop %in% proj_clim$population[r]])
      
    }
    
    else if(cov == "fallR"){
      
      start_date = (proj_clim$date[r]) %m+% months(4)
      end_date = (proj_clim$date[r]) %m+% months(7)
      
      proj_clim$fallR_unscaled[r] = sum(model_data$cum_rain[model_data$date >= start_date & 
                                                              model_data$date < end_date & 
                                                              model_data$pop %in% proj_clim$population[r]])
    }
    
    else if(cov == "prevfallR"){
      
      start_date = (proj_clim$date[r]) %m+% months(-8)
      end_date = (proj_clim$date[r]) %m+% months(-5)
      
      proj_clim$prevfallR_unscaled[r] = sum(model_data$cum_rain[model_data$date >= start_date & 
                                                                  model_data$date < end_date & 
                                                                  model_data$pop %in% proj_clim$population[r]])
    }
    
    else if(cov == "prevwinterR"){
      
      start_date = (proj_clim$date[r]) %m+% months(-4)
      end_date = (proj_clim$date[r]) %m+% months(0)
      
      proj_clim$prevwinterR_unscaled[r] = sum(model_data$cum_rain[model_data$date >= start_date & 
                                                                    model_data$date < end_date & 
                                                                    model_data$pop %in% proj_clim$population[r]])
    }
  }
}


# Standardize climatic variables and remove values outside the
# observed range
proj_clim$summerT = proj_clim$prevwinterT = proj_clim$fallR = proj_clim$prevfallR = proj_clim$prevwinterR =  NA

# Next summer temperature
proj_clim$summerT_unscaled[which(proj_clim$summerT_unscaled > max(droso_anthropogenic$summerT_unscaled))] = max(droso_anthropogenic$summerT_unscaled)
proj_clim$summerT_unscaled[which(proj_clim$summerT_unscaled < min(droso_anthropogenic$summerT_unscaled))] = min(droso_anthropogenic$summerT_unscaled)

proj_clim$summerT = (proj_clim$summerT_unscaled - mean(summerT_timeSeries$summerT_unscaled, na.rm = T)) / (2 * sd(summerT_timeSeries$summerT_unscaled, na.rm = T))

# Previous winter temperature
proj_clim$prevwinterT_unscaled[which(proj_clim$prevwinterT_unscaled > max(droso_anthropogenic$prevwinterT_unscaled))] = max(droso_anthropogenic$prevwinterT_unscaled)
proj_clim$prevwinterT_unscaled[which(proj_clim$prevwinterT_unscaled < min(droso_anthropogenic$prevwinterT_unscaled))] = min(droso_anthropogenic$prevwinterT_unscaled)

proj_clim$prevwinterT = (proj_clim$prevwinterT_unscaled - mean(prevwinterT_timeSeries$prevwinterT_unscaled, na.rm = T)) / (2 * sd(prevwinterT_timeSeries$prevwinterT_unscaled, na.rm = T))

# Next fall rainfall
proj_clim$fallR_unscaled[which(proj_clim$fallR_unscaled > max(droso_anthropogenic$fallR_unscaled))] = max(droso_anthropogenic$fallR_unscaled)
proj_clim$fallR_unscaled[which(proj_clim$fallR_unscaled < min(droso_anthropogenic$fallR_unscaled))] = min(droso_anthropogenic$fallR_unscaled)

proj_clim$fallR = (proj_clim$fallR_unscaled - mean(fallR_timeSeries$fallR_unscaled, na.rm = T)) / (2 * sd(fallR_timeSeries$fallR_unscaled, na.rm = T))

# Previous fall rainfall
proj_clim$prevfallR_unscaled[which(proj_clim$prevfallR_unscaled > max(droso_anthropogenic$prevfallR_unscaled))] = max(droso_anthropogenic$prevfallR_unscaled)
proj_clim$prevfallR_unscaled[which(proj_clim$prevfallR_unscaled < min(droso_anthropogenic$prevfallR_unscaled))] = min(droso_anthropogenic$prevfallR_unscaled)

proj_clim$prevfallR = (proj_clim$prevfallR_unscaled - mean(prevfallR_timeSeries$prevfallR_unscaled, na.rm = T)) / (2 * sd(prevfallR_timeSeries$prevfallR_unscaled, na.rm = T))

# Previous winter rainfall
proj_clim$prevwinterR_unscaled[which(proj_clim$prevwinterR_unscaled > max(droso_anthropogenic$prevwinterR_unscaled))] = max(droso_anthropogenic$prevwinterR_unscaled)
proj_clim$prevwinterR_unscaled[which(proj_clim$prevwinterR_unscaled < min(droso_anthropogenic$prevwinterR_unscaled))] = min(droso_anthropogenic$prevwinterR_unscaled)

proj_clim$prevwinterR = (proj_clim$prevwinterR_unscaled - mean(prevwinterR_timeSeries$prevwinterR_unscaled, na.rm = T)) / (2 * sd(prevwinterR_timeSeries$prevwinterR_unscaled, na.rm = T))




###########################################################################
#
# 4. Individual-based model function ----
#
###########################################################################

## 4.1. Timestep projection function ----
# ----------------------------------

ibm_sim = function(n_sim,             # Number of simulations
                   n_years,           # Number of years per simulation
                   sim,               # Current simulation number
                   first_year,        # Year of start of simulation
                   years_RE,          # Observed years available for random year effect
                   data_initial,      # Initial dataset
                   seedbank_size,     # Initial seedbank size
                   recruitCap,        # Maximum number of recruits per quadrat
                   population,        # Population ID
                   climate_model,     # Name of climate model
                   ibm_data,          # Previously stored projection data
                   seedbank_initial){ # Initial seedbank data
  
  # Year sequence for climatic variable predictions and random effects
  
  years_obs = c(first_year, sample(years_RE, n_years-1, replace = T)) # Random effects
  
  if(climate_model == "Control") years = years_obs # Climatic variable predictions
  # same as RE in the control scenario
  
  else years = c(first_year, seq(first_year + 1, first_year + n_years - 2)) # Climatic variable predictions
  # sequence of 30 years from 2021
  # for the climate-change scenario
  
  
  # Empty files to hold results 
  
  log_meanChangeAboveground = rep(NA, n_years)
  log_lambda = rep(NA, n_years)
  seedbank_size_vec = rep(NA, n_years)
  seedbank_size_vec[1] = seedbank_size
  pop_density = vector(mode = "list", length = n_years)
  extinction = 0
  
  sim_data = NULL # Full individual data across whole simulation
  
  data_droso = data_initial[which(data_initial$time == years_obs[1]), ] # Yearly individual data
  
  # Add density data 
  data_droso$abLarge_unscaled = apply(data_droso, 1, function(x) nrow(data_droso[which(data_droso$quadratID == x[2] & data_droso$size_unscaled > 4.5), ]))
  pop_density[[1]] = aggregate(abLarge_unscaled ~ quadratID, data = data_droso, mean, na.rm = T)$abLarge_unscaled
  
  # Seedbank seeds data
  data_SB = seedbank_initial
  
  # Assign seeds in initial seedbank to quadrats
  data_SB$quadratID = sample(unique(data_droso$quadratID), size = seedbank_size, 
                             replace = T, 
                             prob = as.numeric(table(data_droso$quadratID)/sum(table(data_droso$quadratID))))
  
  # Add density data to seedbank
  data_SB$abLarge_unscaled = apply(data_SB, 1, function(x) nrow(data_droso[which(data_droso$quadratID == x[2] & data_droso$size_unscaled > 4.5), ]))
  
  # Highest seed ID (format Seed_XXX) to give names to new seeds 
  max_seed_ID = max(as.numeric(unlist(lapply(strsplit(data_SB$ID, split = "_"), function(x) x[2]))))
  
  seed_produced = data_droso[0, ] # Data on seeds produced by reproducing plants
  
  time_sim = 1 # Timestep
  sim_data = rbind(sim_data, cbind(data_droso, time_sim)) # Merge full individual 
  # data and yearly individual 
  # data with timestep info
  
  
  # Project the population
  for(i in 2:n_years){
    
    if(ncol(seed_produced) == 0) seed_produced = data_droso[0, ]
    
    # Get new years for projected climatic variables and random effects
    year = years[i-1]
    year_obs = years_obs[i-1]
    
    if(nrow(data_droso) > 0) data_droso$time = (year + 1)
    
    
    ### CLIMATIC VARIABLES ###
    
    # Climatic variables for the control scenario
    if(climate_model == "Control"){ 
      
      summerT = unique(droso_anthropogenic$summerT_unscaled[which(droso_anthropogenic$site == population &
                                                                    droso_anthropogenic$time == year)])
      prevwinterT = unique(droso_anthropogenic$prevwinterT_unscaled[which(droso_anthropogenic$site == population &
                                                                            droso_anthropogenic$time == year)])
      
      fallR = unique(droso_anthropogenic$fallR_unscaled[which(droso_anthropogenic$site == population &
                                                                droso_anthropogenic$time == year)])
      prevfallR = unique(droso_anthropogenic$prevfallR_unscaled[which(droso_anthropogenic$site == population &
                                                                        droso_anthropogenic$time == year)])
      prevwinterR = unique(droso_anthropogenic$prevwinterR_unscaled[which(droso_anthropogenic$site == population &
                                                                            droso_anthropogenic$time == year)])
    }
    
    # Climatic variables for the climate-change scenario
    else{ 
      
      summerT = unique(proj_clim$summerT_unscaled[which(proj_clim$population == population &
                                                          proj_clim$year == year &
                                                          proj_clim$model == climate_model)])
      prevwinterT = unique(proj_clim$prevwinterT_unscaled[which(proj_clim$population == population &
                                                                  proj_clim$year == year &
                                                                  proj_clim$model == climate_model)])
      fallR = unique(proj_clim$fallR_unscaled[which(proj_clim$population == population &
                                                      proj_clim$year == year &
                                                      proj_clim$model == climate_model)])
      prevfallR = unique(proj_clim$prevfallR_unscaled[which(proj_clim$population == population &
                                                              proj_clim$year == year &
                                                              proj_clim$model == climate_model)])
      prevwinterR = unique(proj_clim$prevwinterR_unscaled[which(proj_clim$population == population &
                                                                  proj_clim$year == year &
                                                                  proj_clim$model == climate_model)])
    }
    
    
    ### CORRECTION FACTORS AND NUMBER OF SQUARES ###
    
    # print("CORRECTION FACTORS AND NUMBER OF SQUARES")
    
    corr_seed_surv = seedSurv # Survival of seeds above the ground, from Paniw et al. 2017 (J. Appl. Ecol.)
    
    
    if(nrow(data_droso) > 0){
      
      ### FLOWERING ###
      
      data_droso$flowering = NA
      
      data_droso$flowering = rbinom(n = nrow(data_droso), 
                                    flowering_function(size = data_droso$size_unscaled,
                                                       abLarge = data_droso$abLarge_unscaled,
                                                       prevwinterR = prevwinterR,
                                                       year = year_obs,
                                                       population = population), 
                                    size = 1)
      
      
      ### RECRUITMENT (number of flowers) ###
      
      # Number of flowers per individual reproducing
      
      data_droso$nbFlowers = NA
      
      if(length(data_droso$nbFlowers[which(data_droso$flowering == 1)]) > 0){
        
        data_droso$nbFlowers[which(data_droso$flowering == 1)] = rnbinom(n = nrow(data_droso[which(data_droso$flowering == 1), ]),
                                                                         mu = nbFlowers_function(size = data_droso[which(data_droso$flowering == 1), ]$size_unscaled,
                                                                                                 prevfallR = prevfallR,
                                                                                                 year = year_obs,
                                                                                                 population = population),
                                                                         size = 1)
        
        
        # Number of seeds per individual reproducing
        
        data_droso_sub = data_droso[which(data_droso$flowering == 1), ]
        
        if(any(data_droso_sub$nbFlowers > 0)){
          
          nb_seeds = rpois(nrow(data_droso_sub), lambda = seeds_per_flower) # Number of seeds per flower per individual
          
          seed_produced = data_droso_sub[rep(row.names(data_droso_sub), 
                                             data_droso_sub$nbFlowers * nb_seeds), ] # Dataset of seeds produced
          
          # Format dataset to match other datasets
          seed_produced$size_unscaled = seed_produced$sizeNext = seed_produced$flowering = seed_produced$nbFlowers = seed_produced$survival = NA
          rownames(seed_produced) = seq(1, nrow(seed_produced))
          
          
          # Seeds germinating directly
          
          # print("SEEDS GERMINATING")
          
          seed_produced$goCont = rbinom(n = nrow(seed_produced),
                                        prob = goCont_function(population = population) * corr_seed_surv,
                                        size = 1)
          
          
          # Seeds going to the seedbank
          
          seed_produced$goSB = NA
          
          if(any(seed_produced$goCont == 0)) seed_produced$goSB[which(seed_produced$goCont == 0)] = rbinom(n = nrow(seed_produced[which(seed_produced$goCont == 0), ]),
                                                                                                           prob = (1 - goCont_function(population = population)) * corr_seed_surv,
                                                                                                           size = 1)
          
          
          # Size of the germinated seedlings
          
          # Seedling size parameters (mean, sd, and degrees of freedom)
          seedling_size_parameters = seedling_size_function(prevwinterT = prevwinterT,
                                                            abLarge = seed_produced$abLarge_unscaled,
                                                            year = year_obs, 
                                                            population = population)
          
          # Get seedling size by sampling a truncated Student-t distribution
          # to match the model family, using the corresponding mean, sd, and df
          seed_produced$rownb = seq(1, nrow(seed_produced))
          
          seed_produced$size_unscaled[which(seed_produced$goCont == 1)] = apply(seed_produced[which(seed_produced$goCont == 1), ], 
                                                                                1, 
                                                                                function(x) rtt(1, 
                                                                                                location = seedling_size_parameters$mean[as.numeric(x[grep("rownb", names(x))])], 
                                                                                                scale = seedling_size_parameters$sd, 
                                                                                                df = seedling_size_parameters$df, 
                                                                                                left = 0))
          # Assign zero to individuals with negative size
          if(any(seed_produced$size_unscaled[which(seed_produced$goCont == 1)] == Inf, na.rm = T)) seed_produced$size_unscaled[which(seed_produced$goCont == 1 & seed_produced$size_unscaled == Inf)] = max(seed_produced$size_unscaled[which(seed_produced$goCont == 1 & seed_produced$size_unscaled != Inf)], na.rm = T)
          
          if(any(seed_produced$size_unscaled[which(seed_produced$goCont == 1)] < 0, na.rm = T)) seed_produced$size_unscaled[which(seed_produced$goCont == 1 & seed_produced$size_unscaled < 0)] = 0
          
          
          # Seedling ID
          
          seed_produced$ID = paste("Seed", 
                                   seq(max_seed_ID + 1, max_seed_ID + 1 + nrow(seed_produced) - 1), 
                                   sep = "_")
          
          seed_produced = seed_produced[, colnames(seed_produced)[-which(colnames(seed_produced) == "rownb")]]
          
          # Update max seed ID
          max_seed_ID = max(as.numeric(unlist(lapply(strsplit(seed_produced$ID, split = "_"), function(x) x[2]))))
        }
      } 
      
      
      ### SURVIVAL ###
      
      data_droso$survival = rbinom(n = nrow(data_droso),
                                   prob = survival_function(size = data_droso$size_unscaled,
                                                            abLarge = data_droso$abLarge_unscaled,
                                                            fallR = fallR,
                                                            summerT = summerT,
                                                            year = year_obs,
                                                            population = population),
                                   size = 1)
      
      
      ### GROWTH ###
      
      data_droso$sizeNext = NA
      
      # Growth parameters (mean, sd, and degrees of freedom)
      growth_parameters = growth_function(size = data_droso$size_unscaled, 
                                          summerT = summerT, 
                                          abLarge = data_droso$abLarge_unscaled,
                                          year = year_obs, 
                                          population = population)
      
      # Get size at next timestep by sampling a truncated Student-t distribution
      # to match the model family, using the corresponding mean, sd, and df
      data_droso$rownb = seq(1, nrow(data_droso))
      
      data_droso$sizeNext[which(data_droso$survival == 1)] = apply(data_droso[which(data_droso$survival == 1), ], 
                                                                   1, 
                                                                   function(x) rtt(1, 
                                                                                   location = growth_parameters$mean[as.numeric(x[grep("rownb", names(x))])], 
                                                                                   scale = growth_parameters$sd, 
                                                                                   df = growth_parameters$df, 
                                                                                   left = 0))
      
      # Assign max size in current dataset to individual with infinite size
      if(any(data_droso$sizeNext[which(data_droso$survival == 1)] == Inf, na.rm = T)) data_droso$sizeNext[which(data_droso$survival == 1 & data_droso$sizeNext == Inf)] = max(data_droso$sizeNext[which(data_droso$survival == 1 & data_droso$sizeNext != Inf)], na.rm = T)
      
      if(any(data_droso$sizeNext[which(data_droso$survival == 1)] < 0, na.rm = T)) data_droso$sizeNext[which(data_droso$survival == 1 & data_droso$sizeNext < 0)] = 0                                                          
      
      # Format dataset
      data_droso = data_droso[which(data_droso$survival == 1), colnames(data_droso)[-which(colnames(data_droso) %in% c("rownb"))]]
    }
    
    
    if(nrow(data_SB) > 0){
      
      ### SEEDBANK ###
      
      # Seeds germinating from the seedbank 
      
      data_SB$outSB = rbinom(n = nrow(data_SB),
                             prob = outSB_function(population = population,
                                                   seedSurv = corr_seed_surv),
                             size = 1)
      
      
      # Assign a size to the germinated seeds, add them to the dataset 
      # of aboveground individuals and remove them from the seedbank data
      
      # Seedling size parameters (mean, sd, and degrees of freedom)
      seedling_size_parameters = seedling_size_function(prevwinterT = prevwinterT,
                                                        abLarge = data_SB$abLarge_unscaled,
                                                        year = year_obs, 
                                                        population = population)
      
      # Get seedling size by sampling a truncated Student-t distribution
      # to match the model family, using the corresponding mean, sd, and df
      data_SB$rownb = seq(1, nrow(data_SB))
      
      if(length(which(data_SB$outSB == 1)) > 0){
        
        data_SB$size_unscaled[which(data_SB$outSB == 1)] = apply(data_SB[which(data_SB$outSB == 1), ],
                                                                 1, 
                                                                 function(x) rtt(1, 
                                                                                 location = seedling_size_parameters$mean[as.numeric(x[grep("rownb", names(x))])], 
                                                                                 scale = seedling_size_parameters$sd, 
                                                                                 df = seedling_size_parameters$df, 
                                                                                 left = 0))
        
        # Assign zero in current dataset to individual with negative size
        if(any(data_SB$size_unscaled[which(data_SB$outSB == 1)] < 0, na.rm = T)) data_SB$size_unscaled[which(data_SB$outSB == 1 & data_SB$size_unscaled < 0)] = 0
        
        if(any(data_SB$size_unscaled[which(data_SB$outSB == 1)] == Inf, na.rm = T)) data_SB$size_unscaled[which(data_SB$outSB == 1 & data_SB$size_unscaled == Inf)] = max(data_SB$size_unscaled[which(data_SB$outSB == 1 & data_SB$size_unscaled != Inf)], na.rm = T)                                                                       
      }
      
      # Preparing the seedbank data to merge with the continuous germination data
      
      data_SB$goCont = data_SB$outSB
      data_SB$goSB = NA
      
      # Merge datasets 
      seed_produced = rbind(seed_produced, data_SB[which(data_SB$outSB == 1), colnames(data_SB)[which(colnames(data_SB) %in% colnames(seed_produced))]])
      
      
      
      # Remove germinated seeds from seedbank data
      data_SB = data_SB[which(data_SB$outSB == 0), ]
      
      data_SB = data_SB[, colnames(data_SB)[-which(colnames(data_SB) %in% c("rownb", "outSB"))]]
      
      
      # Seeds staying in and going to the seedbank
      
      data_SB$staySB = rbinom(n = nrow(data_SB),
                              prob = staySB_function(population = population) * corr_seed_surv,
                              size = 1)
      
      # Keep only seeds staying in the seedbank
      data_SB = data_SB[which(data_SB$staySB == 1), ]
      
    }
    
    # Merge seeds going to SB to seedbank data
    seed_produced$staySB = seed_produced$goSB
    
    if(nrow(seed_produced[which(seed_produced$staySB == 1), ]) > 0){
      
      data_SB = rbind(data_SB, seed_produced[which(seed_produced$staySB == 1), ])
    }
    
    if(nrow(data_SB) > 0){
      
      # Format seedbank data to match other datasets
      data_SB = data_SB[, colnames(data_SB)[-which(colnames(data_SB) %in% c("staySB", "outSB", "rownb", "goCont", "goSB"))]]
    }
    
    seedbank_size_vec[i] = nrow(data_SB)
    
    if(nrow(seed_produced) > 0){
      
      # Keep only seeds germinating 
      seed_produced = seed_produced[which(seed_produced$goCont == 1), colnames(seed_produced)[-which(colnames(seed_produced) %in% c("goCont", "goSB", "staySB"))]]
      
      
      # Cap the number of recruits if needed
      if(!is.null(recruitCap) & nrow(seed_produced) > 0){
        
        # Get number of seedlings per quadrat
        quadratsAboveMaxSeedlings = aggregate(ID ~ quadratID, 
                                              data = seed_produced,
                                              FUN = function(x) length(x))
        
        # Keep quadrats where the number of seedlings is above the threshold
        quadratsAboveMaxSeedlings = quadratsAboveMaxSeedlings[which(quadratsAboveMaxSeedlings$ID > recruitCap), ]
        
        # If quadrats are above the threshold, sample the seedlings that
        # will be kept
        if(nrow(quadratsAboveMaxSeedlings) > 0){
          
          for(quadrat in quadratsAboveMaxSeedlings$quadratID){
            
            recruitsKept_ID = sample(seq(1, nrow(seed_produced[which(seed_produced$quadratID == quadrat), ])), 
                                     size = recruitCap, replace = F)
            
            recruitsKept = seed_produced[which(seed_produced$quadratID == quadrat)[recruitsKept_ID], ]
            
            seed_produced = seed_produced[-which(seed_produced$quadratID == quadrat), ]
            seed_produced = rbind(seed_produced, recruitsKept)
          }
        }
      }
      
      
      # Adding new seedlings to the population
      
      data_droso = rbind(data_droso, seed_produced)
    }
    
    if(nrow(data_droso) > 0) data_droso$time = (year + 1)
    
    log_meanChangeAboveground[i] = log(nrow(data_droso)/nrow(sim_data[which(sim_data$time_sim == i-1), ])) # Calculate mean change in aboveground population abundance
    log_lambda[i] = log((nrow(data_droso) + nrow(data_SB))/(nrow(sim_data[which(sim_data$time_sim == i-1), ]) + seedbank_size_vec[i-1])) # Calculate log lambda (with seedbank)
    
    
    ### SAVE DATA ###
    
    time_sim = i # Update timestep
    
    # Merge yearly individual data with full individual data
    if(nrow(data_droso) > 0) data_droso = cbind(data_droso, time_sim) 
    sim_data = rbind(sim_data, data_droso[, colnames(data_droso)[which(colnames(data_droso) %in% colnames(sim_data))]])
    
    # Format yearly data
    data_droso = data_droso[, -ncol(data_droso)]
    
    # Assess population quasi-extinction
    if(nrow(data_droso) < 5 & nrow(data_SB) < 50){
      
      extinction = 1
      break
      
    }
    
    
    ##### NEW  DATA FOR T + 1 #####
    
    data_droso$size_unscaled[which(!is.na(data_droso$sizeNext))] = data_droso$sizeNext[which(!is.na(data_droso$sizeNext))] # Assign new size to individuals
    
    ### DENSITY ###
    
    # Get quadrat-specific density from number of individuals with size > 4.5 
    data_droso$abLarge_unscaled = as.numeric(apply(data_droso, 1, FUN = function(x) nrow(data_droso[which(data_droso$quadratID == x[2] & data_droso$size_unscaled > 4.5), ])))
    if(nrow(data_droso) > 0) pop_density[[i]] = aggregate(abLarge_unscaled ~ quadratID, data = data_droso, mean, na.rm = T)$abLarge_unscaled
    
    
  }
  
  # Create summary data
  data_agg = aggregate(ID ~ time_sim, data = sim_data, function(x) length(x))
  
  data_agg$run = sim
  
  ibm_data = rbind(ibm_data, data_agg)
  
  
  return(list(pop_data = sim_data,
              pop_size = ibm_data,
              pop_density = pop_density,
              log_meanChangeAboveground = log_meanChangeAboveground,
              log_lambda = log_lambda,
              seedbank_size_vec = seedbank_size_vec,
              extinction = extinction))
}


## 4.2. Simulation function ----
# -------------------------

ibm_anthropogenic = function(n_sim = 1000,                # Number of simulations
                             n_years = 50,                # Number of years per simulation
                             first_year = 2021,           # Year of start of simulation
                             years_RE = seq(2016, 2021),  # Observed years available for random year effect
                             data_initial,                # Initial dataset
                             seedbank_size = 3000,        # Initial seedbank size
                             recruitCap,                  # Maximum number of recruits per quadrat
                             population,                  # Population ID
                             climate_model){              # Name of climate model
  
  
  # Initialization - Prepare storing objects
  
  extinction_vector = rep(0, n_sim)
  log_meanChangeAboveground_array = array(NA, dim = c(n_sim, n_years))
  log_lambda_array = array(NA, dim = c(n_sim, n_years))
  seedbank_size_array = array(NA, dim = c(n_sim, n_years))
  pop_size_list = vector(mode = "list", length = n_sim)
  pop_data_list = vector(mode = "list", length = n_sim)
  pop_density_list = vector(mode = "list", length = n_sim)
  
  ibm_data = NULL
  
  # Calculate number of flowers
  data_initial$nbFlowers = data_initial$fs * data_initial$fps
  
  # Format initial data
  data_initial = data_initial[, c("site", "quadratID", "ID", "size_unscaled", "sizeNext", 
                                  "fl", "nbFlowers", 
                                  "surv", "time", "abLarge_unscaled")]
  
  colnames(data_initial) = c("site", "quadratID", "ID", "size_unscaled", "sizeNext", "flowering", "nbFlowers", 
                             "survival", "time", "abLarge_unscaled")
  
  # Prepare initial seedbank data
  seedbank_initial = data.frame(site = population, quadratID = NA, ID = paste("Seed", seq(1, seedbank_size), sep = "_"),
                                size_unscaled = NA, sizeNext = NA, flowering = NA, nbFlowers = NA,
                                survival = NA, time = NA, abLarge_unscaled = NA)
  
  
  # Projection
  
  # For each simulation, project the population
  # If there is an error (exploding population), restart the simulation
  # to have another sequence of years
  
  for(sim in 1:n_sim){
    
    print(paste("Iteration", sim))
    
    while(TRUE){
      
      ibm_sim_result =
        try(ibm_sim(n_sim = n_sim,
                    n_years = n_years,
                    sim = sim,
                    first_year = first_year,
                    years_RE = years_RE,
                    seedbank_size = seedbank_size,
                    recruitCap = recruitCap,
                    population = population,
                    climate_model = climate_model,
                    ibm_data = ibm_data,
                    data_initial = data_initial,
                    seedbank_initial = seedbank_initial),
            silent = TRUE)
      
      if(!is(ibm_sim_result, 'try-error')) break
    }
    
    # Fill in result objects
    log_meanChangeAboveground_array[sim, ] = ibm_sim_result$log_meanChangeAboveground
    log_lambda_array[sim, ] = ibm_sim_result$log_lambda
    seedbank_size_array[sim, ] = ibm_sim_result$seedbank_size_vec
    extinction_vector[sim] = ibm_sim_result$extinction
    pop_size_list[[sim]] = ibm_sim_result$pop_size
    pop_data_list[[sim]] = ibm_sim_result$pop_data
    pop_density_list[[sim]] = ibm_sim_result$pop_density
    
  }
  
  return(list(pop_data = pop_data_list,
              pop_size = pop_size_list,
              pop_density = pop_density_list,
              log_meanChangeAboveground = log_meanChangeAboveground_array,
              log_lambda = log_lambda_array,
              seedbank_size = seedbank_size_array,
              extinction = extinction_vector))
}




###########################################################################
#
# 5. Projections ----
#
###########################################################################

n_sim = 500
n_years = 30

## 5.1. Retin ----
# -----------

population = "Retin"

recruitCap = aggregate(ID ~ quadratID + time,
                       data = droso_anthropogenic[which(droso_anthropogenic$site == population &
                                                          droso_anthropogenic$stage == "SD"), ], FUN = function(x) length(x))

recruitCap = round(max(recruitCap$ID))


# Function to run the projections in parallel
ncpus <- 5

startIBM <- function(sim){
  
  if (sim > 1) {      #slave, give it other random seed,
    rm(".Random.seed", envir = .GlobalEnv); runif(1)
  }
  
  ibm_anthropogenic(n_sim = n_sim/ncpus,
                    n_years = n_years,
                    first_year = 2021, 
                    data_initial = droso_anthropogenic[which(droso_anthropogenic$site == population), ], 
                    recruitCap = recruitCap, 
                    population = population, 
                    climate_model = climate_model)
  
}


## 5.1.1. Control ----
# ---------------

climate_model = "Control"


# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, 
       slaveOutfile = paste0("Output/Projections/IBM_Progress_", 
                             population, "_", climate_model, ".txt"))

# Export data and libraries to each parallel core
sfExport(list = c(ls(), ".Random.seed")) 

sfLibrary("mgcv", character.only = TRUE)
sfLibrary("crch", character.only = TRUE)
sfLibrary("snowfall", character.only = TRUE)

# Run projections in parallel
ibm_retin_control = sfClusterApplyLB(1:ncpus, startIBM)

# Save results and stop cluster
save(ibm_retin_control, 
     file = paste0("Output/Projections/IBM_Anthropogenic_", 
                   population, "_", climate_model, ".RData"))
sfStop()

rm(ibm_retin_control)


## 5.1.2. CanESM5 ----
# ---------------

climate_model = "CanESM5"

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, 
       slaveOutfile = paste0("Output/Projections/IBM_Progress_", 
                             population, "_", climate_model, ".txt"))

# Export data and libraries to each parallel core
sfExport(list = c(ls(), ".Random.seed")) 

sfLibrary("mgcv", character.only = TRUE)
sfLibrary("crch", character.only = TRUE)
sfLibrary("snowfall", character.only = TRUE)

# Run projections in parallel
ibm_retin_canesm5 = sfClusterApplyLB(1:ncpus, startIBM)

# Save results and stop cluster
save(ibm_retin_canesm5, 
     file = paste0("Output/Projections/IBM_Anthropogenic_", 
                   population, "_", climate_model, ".RData"))
sfStop()

rm(ibm_retin_canesm5)


## 5.1.3. EC_Earth3 ----
# -----------------

climate_model = "EC_Earth3"

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, 
       slaveOutfile = paste0("Output/Projections/IBM_Progress_", 
                             population, "_", climate_model, ".txt"))

# Export data and libraries to each parallel core
sfExport(list = c(ls(), ".Random.seed")) 

sfLibrary("mgcv", character.only = TRUE)
sfLibrary("crch", character.only = TRUE)
sfLibrary("snowfall", character.only = TRUE)

# Run projections in parallel
ibm_retin_ec_earth3 = sfClusterApplyLB(1:ncpus, startIBM)

# Save results and stop cluster
save(ibm_retin_ec_earth3,
     file = paste0("Output/Projections/IBM_Anthropogenic_", 
                   population, "_", climate_model, ".RData"))
sfStop()

rm(ibm_retin_ec_earth3)


## 5.1.4. FGOALS_G3 ----
# -----------------

climate_model = "FGOALS_G3"

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, 
       slaveOutfile = paste0("Output/Projections/IBM_Progress_", 
                             population, "_", climate_model, ".txt"))

# Export data and libraries to each parallel core
sfExport(list = c(ls(), ".Random.seed")) 

sfLibrary("mgcv", character.only = TRUE)
sfLibrary("crch", character.only = TRUE)
sfLibrary("snowfall", character.only = TRUE)

# Run projections in parallel
ibm_retin_fgoals_g3 = sfClusterApplyLB(1:ncpus, startIBM)

# Save results and stop cluster
save(ibm_retin_fgoals_g3, 
     file = paste0("Output/Projections/IBM_Anthropogenic_",
                   population, "_", climate_model, ".RData"))
sfStop()

rm(ibm_retin_fgoals_g3)


## 5.1.5. GFDL_ESM4 ----
# -----------------

climate_model = "GFDL_ESM4"

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, 
       slaveOutfile = paste0("Output/Projections/IBM_Progress_", 
                             population, "_", climate_model, ".txt"))

# Export data and libraries to each parallel core
sfExport(list = c(ls(), ".Random.seed")) 

sfLibrary("mgcv", character.only = TRUE)
sfLibrary("crch", character.only = TRUE)
sfLibrary("snowfall", character.only = TRUE)

# Run projections in parallel
ibm_retin_gfdl_esm4 = sfClusterApplyLB(1:ncpus, startIBM)

# Save results and stop cluster
save(ibm_retin_gfdl_esm4, 
     file = paste0("Output/Projections/IBM_Anthropogenic_", 
                   population, "_", climate_model, ".RData"))
sfStop()

rm(ibm_retin_gfdl_esm4)


## 5.1.6. GISS_E2_1_G ----
# -------------------

climate_model = "GISS_E2_1_G"

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, 
       slaveOutfile = paste0("Output/Projections/IBM_Progress_", 
                             population, "_", climate_model, ".txt"))

# Export data and libraries to each parallel core
sfExport(list = c(ls(), ".Random.seed")) 

sfLibrary("mgcv", character.only = TRUE)
sfLibrary("crch", character.only = TRUE)
sfLibrary("snowfall", character.only = TRUE)

# Run projections in parallel
ibm_retin_giss_e2_1_g = sfClusterApplyLB(1:ncpus, startIBM)

# Save results and stop cluster
save(ibm_retin_giss_e2_1_g, 
     file = paste0("Output/Projections/IBM_Anthropogenic_", 
                   population, "_", climate_model, ".RData"))
sfStop()

rm(ibm_retin_giss_e2_1_g)


## 5.1.7. INM_CM4_8 ----
# -----------------

climate_model = "INM_CM4_8"

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, 
       slaveOutfile = paste0("Output/Projections/IBM_Progress_", 
                             population, "_", climate_model, ".txt"))

# Export data and libraries to each parallel core
sfExport(list = c(ls(), ".Random.seed")) 

sfLibrary("mgcv", character.only = TRUE)
sfLibrary("crch", character.only = TRUE)
sfLibrary("snowfall", character.only = TRUE)

# Run projections in parallel
ibm_retin_inm_cm4_8 = sfClusterApplyLB(1:ncpus, startIBM)

# Save results and stop cluster
save(ibm_retin_inm_cm4_8, 
     file = paste0("Output/Projections/IBM_Anthropogenic_", 
                   population, "_", climate_model, ".RData"))
sfStop()

rm(ibm_retin_inm_cm4_8)


## 5.1.8. IPSL_CM6A_LR ----
# --------------------

climate_model = "IPSL_CM6A_LR"

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, 
       slaveOutfile = paste0("Output/Projections/IBM_Progress_", 
                             population, "_", climate_model, ".txt"))

# Export data and libraries to each parallel core
sfExport(list = c(ls(), ".Random.seed")) 

sfLibrary("mgcv", character.only = TRUE)
sfLibrary("crch", character.only = TRUE)
sfLibrary("snowfall", character.only = TRUE)

# Run projections in parallel
ibm_retin_ipsl_cm6a_lr = sfClusterApplyLB(1:ncpus, startIBM)

# Save results and stop cluster
save(ibm_retin_ipsl_cm6a_lr, 
     file = paste0("Output/Projections/IBM_Anthropogenic_", 
                   population, "_", climate_model, ".RData"))
sfStop()

rm(ibm_retin_ipsl_cm6a_lr)


## 5.1.9. MIROC6 ----
# --------------

climate_model = "MIROC6"

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, 
       slaveOutfile = paste0("Output/Projections/IBM_Progress_", 
                             population, "_", climate_model, ".txt"))

# Export data and libraries to each parallel core
sfExport(list = c(ls(), ".Random.seed")) 

sfLibrary("mgcv", character.only = TRUE)
sfLibrary("crch", character.only = TRUE)
sfLibrary("snowfall", character.only = TRUE)

# Run projections in parallel
ibm_retin_miroc6 = sfClusterApplyLB(1:ncpus, startIBM)

# Save results and stop cluster
save(ibm_retin_miroc6, 
     file = paste0("Output/Projections/IBM_Anthropogenic_", 
                   population, "_", climate_model, ".RData"))
sfStop()

rm(ibm_retin_miroc6)


## 5.1.10. MPI_ESM1_2_LR ----
# ----------------------

climate_model = "MPI_ESM1_2_LR"

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, 
       slaveOutfile = paste0("Output/Projections/IBM_Progress_", 
                             population, "_", climate_model, ".txt"))

# Export data and libraries to each parallel core
sfExport(list = c(ls(), ".Random.seed")) 

sfLibrary("mgcv", character.only = TRUE)
sfLibrary("crch", character.only = TRUE)
sfLibrary("snowfall", character.only = TRUE)

# Run projections in parallel
ibm_retin_mpi_esm1_2_lr = sfClusterApplyLB(1:ncpus, startIBM)

# Save results and stop cluster
save(ibm_retin_mpi_esm1_2_lr, 
     file = paste0("Output/Projections/IBM_Anthropogenic_", 
                   population, "_", climate_model, ".RData"))
sfStop()

rm(ibm_retin_mpi_esm1_2_lr)


## 5.1.11. MRI_ESM2_0 ----
# -------------------

climate_model = "MRI_ESM2_0"

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, 
       slaveOutfile = paste0("Output/Projections/IBM_Progress_", 
                             population, "_", climate_model, ".txt"))

# Export data and libraries to each parallel core
sfExport(list = c(ls(), ".Random.seed")) 

sfLibrary("mgcv", character.only = TRUE)
sfLibrary("crch", character.only = TRUE)
sfLibrary("snowfall", character.only = TRUE)

# Run projections in parallel
ibm_retin_mri_esm2_0 = sfClusterApplyLB(1:ncpus, startIBM)

# Save results and stop cluster
save(ibm_retin_mri_esm2_0, 
     file = paste0("Output/Projections/IBM_Anthropogenic_",
                   population, "_", climate_model, ".RData"))
sfStop()

rm(ibm_retin_mri_esm2_0)


## 5.1.12. NorESM2_MM ----
# -------------------

climate_model = "NorESM2_MM"

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, 
       slaveOutfile = paste0("Output/Projections/IBM_Progress_", 
                             population, "_", climate_model, ".txt"))

# Export data and libraries to each parallel core
sfExport(list = c(ls(), ".Random.seed")) 

sfLibrary("mgcv", character.only = TRUE)
sfLibrary("crch", character.only = TRUE)
sfLibrary("snowfall", character.only = TRUE)

# Run projections in parallel
ibm_retin_noresm2_mm = sfClusterApplyLB(1:ncpus, startIBM)

# Save results and stop cluster
save(ibm_retin_noresm2_mm, 
     file = paste0("Output/Projections/IBM_Anthropogenic_", 
                   population, "_", climate_model, ".RData"))
sfStop()

rm(ibm_retin_noresm2_mm)