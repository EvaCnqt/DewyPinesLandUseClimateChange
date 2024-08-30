############################################################################
#
# This script uses the observed years not used for the models
# to validate the vital-rate models and the IBM.
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


## 1.2. Loading data ----
# ------------------

library(crch)


## 1.3. Loading data ----
# ------------------

# Dewy-pine data
droso_natural = read.csv("Data/droso_natural.csv")
droso_natural_full = read.csv("Data/droso_natural_full.csv")

droso_seedbank = read.csv("Data/droso_seedbank_natural.csv")

seeds_per_flower = 9.8


# Vital-rate models
load("Output/Models/Survival_GAM_Natural.RData")
load("Output/Models/Growth_GAM_Natural.RData")
load("Output/Models/FloweringProb_GAM_Natural.RData")
load("Output/Models/NbFlowers_GAM_Natural.RData")
load("Output/Models/SeedlingSize_GAM_Natural.RData")


# Average density per square for covariate standardization
nbSquares = aggregate(quadratID ~ time + site, 
                      data = droso_natural, 
                      function(x) length(unique(x)))
density_per_square = aggregate(abLarge_unscaled ~ quadratID + time + site, 
                               data = droso_natural, 
                               function(x) unique(x))
yearly_density_per_square = aggregate(abLarge_unscaled ~ time + site, 
                                      data = density_per_square, 
                                      function(x) sum(x))
yearly_density_per_square$abLarge_unscaled = yearly_density_per_square$abLarge_unscaled/nbSquares$quadratID


# Year- and population-specific climatic variables for covariate standardization
summerT_timeSeries = aggregate(summerT_unscaled ~ time + site, 
                               data = droso_natural, mean)
prevwinterT_timeSeries = aggregate(prevwinterT_unscaled ~ time + site, 
                                   data = droso_natural, mean)
fallR_timeSeries = aggregate(fallR_unscaled ~ time + site, 
                             data = droso_natural, mean)
prevfallR_timeSeries = aggregate(prevfallR_unscaled ~ time + site, 
                                 data = droso_natural, mean)




###########################################################################
#
# 2. Building vital-rate functions ----
#
###########################################################################

## 2.1. Survival ----
# --------------

survival_function = function(size, abLarge, TSFcont, fallR, summerT, year, population){
  
  # Standardize covariates
  size_scaled = (size - mean(droso_natural$size_unscaled, na.rm = T)) / (2 * sd(droso_natural$size_unscaled, na.rm = T))
  abLarge_scaled = (abLarge - mean(yearly_density_per_square$abLarge_unscaled, na.rm = T)) / (2 * sd(yearly_density_per_square$abLarge_unscaled, na.rm = T))
  TSFcont_scaled = (TSFcont - mean(droso_natural$TSFcont_unscaled, na.rm = T)) / (2 * sd(droso_natural$TSFcont_unscaled, na.rm = T))
  fallR_scaled = (fallR - mean(fallR_timeSeries$fallR_unscaled, na.rm = T)) / (2 * sd(fallR_timeSeries$fallR_unscaled, na.rm = T))
  summerT_scaled = (summerT - mean(summerT_timeSeries$summerT_unscaled, na.rm = T)) / (2 * sd(summerT_timeSeries$summerT_unscaled, na.rm = T))

  # Calculate survival
  survival = predict(surv_natural, newdata = data.frame(size = size_scaled,
                                                   abLarge = abLarge_scaled,
                                                   TSFcont = TSFcont_scaled,
                                                   fallR = fallR_scaled,
                                                   summerT = summerT_scaled,
                                                   time = year,
                                                   site = population), type = "response")
  
  return(survival)
}


## 2.2. Growth ----
# ------------

growth_function = function(size, abLarge, TSFcont, fallR, year, population){
  
  # Standardize covariates
  size_scaled = (size - mean(droso_natural$size_unscaled, na.rm = T)) / (2 * sd(droso_natural$size_unscaled, na.rm = T))
  abLarge_scaled = (abLarge - mean(yearly_density_per_square$abLarge_unscaled, na.rm = T)) / (2 * sd(yearly_density_per_square$abLarge_unscaled, na.rm = T))
  TSFcont_scaled = (TSFcont - mean(droso_natural$TSFcont_unscaled, na.rm = T)) / (2 * sd(droso_natural$TSFcont_unscaled, na.rm = T))
  fallR_scaled = (fallR - mean(fallR_timeSeries$fallR_unscaled, na.rm = T)) / (2 * sd(fallR_timeSeries$fallR_unscaled, na.rm = T))

  # Get parameters from growth model (mean, sd, and degrees of freedom)
  growth_mean = predict(growth_natural, newdata = data.frame(size = size_scaled,
                                                        abLarge = abLarge_scaled,
                                                        TSFcont = TSFcont_scaled,
                                                        fallR = fallR_scaled,
                                                        time = year,
                                                        site = population), type = "response")
  
  growth_sd = family(growth_natural)$getTheta(trans = T)[2] # Using trans = T, the second value of theta is sigma, the standard deviation (https://stats.stackexchange.com/questions/550339/extracting-the-degrees-of-freedom-of-t-distribution-of-a-gam)
  
  growth_df = family(growth_natural)$getTheta(trans = T)[1] # Using trans = T, the first value of theta is nu, the degrees of freedom (https://stats.stackexchange.com/questions/550339/extracting-the-degrees-of-freedom-of-t-distribution-of-a-gam)
  
  return(list(mean = growth_mean, 
              sd = growth_sd, 
              df = growth_df))
}


## 2.3. Seedling sizes ----
# --------------------

seedling_size_function = function(abLarge, prevwinterT, TSFcont, year, population){
  
  # Standardize covariates
  abLarge_scaled = (abLarge - mean(yearly_density_per_square$abLarge_unscaled, na.rm = T)) / (2 * sd(yearly_density_per_square$abLarge_unscaled, na.rm = T))
  prevwinterT_scaled = (prevwinterT - mean(prevwinterT_timeSeries$prevwinterT_unscaled, na.rm = T)) / (2 * sd(prevwinterT_timeSeries$prevwinterT_unscaled, na.rm = T))
  TSFcont_scaled = (TSFcont - mean(droso_natural$TSFcont_unscaled, na.rm = T)) / (2 * sd(droso_natural$TSFcont_unscaled, na.rm = T))
  
  # Get parameters from seedling size model (mean, sd, and degrees of freedom)
  seedling_size_mean = as.numeric(predict(seedlingSize_natural, newdata = data.frame(abLarge = abLarge_scaled,
                                                                                prevwinterT = prevwinterT_scaled,
                                                                                TSFcont = TSFcont_scaled,
                                                                                time = year,
                                                                                site = population), type = "response"))
  
  seedling_size_sd = family(seedlingSize_natural)$getTheta(trans = T)[2] # Using trans = T, the second value of theta is sigma, the standard deviation (https://stats.stackexchange.com/questions/550339/extracting-the-degrees-of-freedom-of-t-distribution-of-a-gam)
  
  seedling_size_df = family(seedlingSize_natural)$getTheta(trans = T)[1] # Using trans = T, the first value of theta is nu, the degrees of freedom (https://stats.stackexchange.com/questions/550339/extracting-the-degrees-of-freedom-of-t-distribution-of-a-gam)
  
  return(list(mean = seedling_size_mean, 
              sd = seedling_size_sd, 
              df = seedling_size_df))
  
}


## 2.4. Flowering probability ----
# ---------------------------

flowering_function = function(size, abLarge, prevwinterT, prevfallR, TSFcont, year, population){
  
  # Standardize covariates
  size_scaled = (size - mean(droso_natural$size_unscaled, na.rm = T)) / (2 * sd(droso_natural$size_unscaled, na.rm = T))
  abLarge_scaled = (abLarge - mean(yearly_density_per_square$abLarge_unscaled, na.rm = T)) / (2 * sd(yearly_density_per_square$abLarge_unscaled, na.rm = T))
  prevwinterT_scaled = (prevwinterT - mean(prevwinterT_timeSeries$prevwinterT_unscaled, na.rm = T)) / (2 * sd(prevwinterT_timeSeries$prevwinterT_unscaled, na.rm = T))
  prevfallR_scaled = (prevfallR - mean(prevfallR_timeSeries$prevfallR_unscaled, na.rm = T)) / (2 * sd(prevfallR_timeSeries$prevfallR_unscaled, na.rm = T))
  TSFcont_scaled = (TSFcont - mean(droso_natural$TSFcont_unscaled, na.rm = T)) / (2 * sd(droso_natural$TSFcont_unscaled, na.rm = T))
  
  # Calculate flowering probability
  flowering = predict(flowering_natural, newdata = data.frame(size = size_scaled,
                                                         abLarge = abLarge_scaled,
                                                         prevwinterT = prevwinterT_scaled,
                                                         prevfallR = prevfallR_scaled,
                                                         TSFcont = TSFcont_scaled,
                                                         time = year,
                                                         site = population), type = "response")
  
  return(flowering)
}


## 2.5. Number of flowers ----
# -----------------------

nbFlowers_function = function(size, prevwinterT, TSFcont, year, population){
  
  # Standardize covariates
  size_scaled = (size - mean(droso_natural$size_unscaled, na.rm = T)) / (2 * sd(droso_natural$size_unscaled, na.rm = T))
  prevwinterT_scaled = (prevwinterT - mean(prevwinterT_timeSeries$prevwinterT_unscaled, na.rm = T)) / (2 * sd(prevwinterT_timeSeries$prevwinterT_unscaled, na.rm = T))
  TSFcont_scaled = (TSFcont - mean(droso_natural$TSFcont_unscaled, na.rm = T)) / (2 * sd(droso_natural$TSFcont_unscaled, na.rm = T))
  
  # Calculate number of flowers
  nbFlowers = predict(nbFlow_natural, newdata = data.frame(size = size_scaled,
                                                      prevwinterT = prevwinterT_scaled,
                                                      TSFcont = TSFcont_scaled,
                                                      time = year,
                                                      site = population), type = "response")
  
  return(nbFlowers)
}


## 2.6. Immediate germination (goCont) ----
# ------------------------------------

goCont_function = function(TSF){
  
  if(TSF %in% seq(0, 5)){
    
    goCont = droso_seedbank$value[which(droso_seedbank$vital_rate == "goCont" &
                                          droso_seedbank$TSF == TSF)]
  }
  
  else{
    
    goCont = droso_seedbank$value[which(droso_seedbank$vital_rate == "goCont" &
                                          droso_seedbank$TSF == 5)]
    
  }
  
  return(goCont)
}


## 2.7. Staying in the seed bank (staySB) ----
# ---------------------------------------

staySB_function = function(TSF){
  
  if(TSF %in% seq(0, 5)){
    
    staySB = droso_seedbank$value[which(droso_seedbank$vital_rate == "staySB" &
                                        droso_seedbank$TSF == TSF)]
    
  }
  
  else{
    
    staySB = droso_seedbank$value[which(droso_seedbank$vital_rate == "staySB" &
                                        droso_seedbank$TSF == 5)]
    
  }
  
  return(staySB)
}

## 2.8. Germinating out of the seed bank (outSB) ----
# ----------------------------------------------

outSB_function = function(TSF){
  
  if(TSF %in% seq(0, 5)){
    
    outSB = droso_seedbank$value[which(droso_seedbank$vital_rate == "outSB" &
                                         droso_seedbank$TSF == TSF)]
    
  }
  
  else{
    
    outSB = droso_seedbank$value[which(droso_seedbank$vital_rate == "outSB" &
                                         droso_seedbank$TSF == 5)]
  }
  
  return(outSB)
}




###########################################################################
#
# 3. Individual-based model function ----
#
###########################################################################

## 3.1. Timestep projection function ----
# ----------------------------------

ibm_sim = function(n_sim,              # Number of simulations
                   n_years,            # Number of years per simulation
                   sim = sim,          # Current simulation number
                   years_RE,           # Observed years available for random year effect
                   seedbank_size,      # Initial seedbank size
                   recruitCap,         # Maximum number of recruits per quadrat
                   max_nbFlowers,      # Maximum number of flowers per individual
                   min_seedlingSize,
                   population,         # Population ID
                   ibm_data,           # Previously stored projection data
                   data_initial,       # Initial dataset
                   seedbank_initial){  # Initial seedbank data
  
  # Year sequence for climatic variable predictions and random effects
  
  years = years_obs = years_RE # Random effects
  
  # Empty files to hold results 
  
  log_meanChangeAboveground = rep(NA, n_years)
  log_lambda = rep(NA, n_years)
  seedbank_size_vec = rep(NA, n_years)
  seedbank_size_vec[1] = seedbank_size
  pop_density = vector(mode = "list", length = n_years)
  
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
  
  time_sim = 1 # Timestep
  sim_data = rbind(sim_data, cbind(data_droso, time_sim)) # Merge full individual 
  # data and yearly individual 
  # data with timestep info
  
  # Building the post-fire environmental states succession vector
  states = unique(data_initial$TSFcont_unscaled)
  states[which(states > 4)] = 4
  TSFcont = unique(data_initial$TSFcont_unscaled)
  
  
  # Project the population
  for(i in 2:n_years){
    
    seed_produced = data_droso[0, ] # Data on seeds produced by reproducing plants
    seed_produced[1, ] = NA
    seed_produced$goCont = seed_produced$goSB = seed_produced$staySB = seed_produced$outSB = seed_produced$rownb = NA
    seed_produced = seed_produced[0, ]
    
    
    # Get new years for projected climatic variables and random effects
    year = years[i-1]
    year_obs = years_obs[i-1]
    
    if(nrow(data_droso) > 0) data_droso$time = (year + 1)
    
    
    ### CLIMATIC VARIABLES ###
    
    # Climatic variables for the control scenario
    summerT = unique(droso_natural_full$summerT_unscaled[which(droso_natural_full$site == population &
                                                               droso_natural_full$time == year)])
    prevwinterT = unique(droso_natural_full$prevwinterT_unscaled[which(droso_natural_full$site == population &
                                                                       droso_natural_full$time == year)])
    
    fallR = unique(droso_natural_full$fallR_unscaled[which(droso_natural_full$site == population &
                                                           droso_natural_full$time == year)])
    prevfallR = unique(droso_natural_full$prevfallR_unscaled[which(droso_natural_full$site == population &
                                                                   droso_natural_full$time == year)])
    
    
    ### TSF, CORRECTION FACTORS AND NUMBER OF SQUARES ###
    
    TSFcat = states[i-1]
    TSF = TSFcont[i-1]
    if(nrow(data_droso) > 0) data_droso$TSFcont_unscaled = TSF
    
    if(states[i-1] != 0){
      
      if(states[i-1] > 1){
        
        if(nrow(data_droso) > 0){
          
          ### FLOWERING ###
          
          data_droso$flowering = NA
          
          data_droso$flowering = rbinom(n = nrow(data_droso), 
                                        flowering_function(size = data_droso$size_unscaled,
                                                           abLarge = data_droso$abLarge_unscaled,
                                                           prevwinterT = prevwinterT,
                                                           prevfallR = prevfallR,
                                                           TSFcont = TSF,
                                                           year = year_obs,
                                                           population = population), 
                                        size = 1)
          

          ### RECRUITMENT (number of flowers) ###
          
          # Number of flowers per individual reproducing
          
          data_droso$nbFlowers = NA
          
          if(length(data_droso$nbFlowers[which(data_droso$flowering == 1)]) > 0){
            
            data_droso$nbFlowers[which(data_droso$flowering == 1)] = rnbinom(n = nrow(data_droso[which(data_droso$flowering == 1), ]),
                                                                             mu = nbFlowers_function(size = data_droso[which(data_droso$flowering == 1), ]$size_unscaled,
                                                                                                     prevwinterT = prevwinterT,
                                                                                                     TSFcont = TSF,
                                                                                                     year = year_obs,
                                                                                                     population = population),
                                                                             size = 1)
            
            
            # Cap the number of flowers
            if(!is.null(max_nbFlowers)){
              
              if(any(data_droso$nbFlowers[which(data_droso$flowering == 1)] > max_nbFlowers$nbFlowers[which(max_nbFlowers$TSFcont == TSFcat)], na.rm = T)) data_droso$nbFlowers[which(data_droso$flowering == 1 & data_droso$nbFlowers > max_nbFlowers$nbFlowers[which(max_nbFlowers$TSFcont == TSFcat)])] = max_nbFlowers$nbFlowers[which(max_nbFlowers$TSFcont == TSFcat)]
            }
            
            
            # Number of seeds per individual reproducing
            
            data_droso_sub = data_droso[which(data_droso$flowering == 1 & 
                                              data_droso$nbFlowers > 0), ]
            
            if(any(data_droso_sub$nbFlowers > 0)){
              
              nb_seeds = rpois(nrow(data_droso_sub), lambda = seeds_per_flower) # Number of seeds per flower per individual
              
              seed_produced = data_droso_sub[rep(row.names(data_droso_sub), 
                                                 data_droso_sub$nbFlowers * nb_seeds), ] # Dataset of seeds produced
              
              # Format dataset to match other datasets
              seed_produced$size_unscaled = seed_produced$sizeNext = seed_produced$flowering = seed_produced$nbFlowers = seed_produced$survival = NA
              rownames(seed_produced) = seq(1, nrow(seed_produced))
              
              
              # Seeds germinating directly
              
              seed_produced$goCont = rbinom(n = nrow(seed_produced),
                                            prob = goCont_function(TSF = TSF),
                                            size = 1)
              
              
              # Seeds going to the seedbank
              
              seed_produced$goSB = NA
              
              if(any(seed_produced$goCont == 0)) seed_produced$goSB[which(seed_produced$goCont == 0)] = rbinom(n = nrow(seed_produced[which(seed_produced$goCont == 0), ]),
                                                                                                               prob = (1 - goCont_function(TSF = TSF)),
                                                                                                               size = 1)
              
              
              # Size of the germinated seedlings
              
              # Seedling size parameters (mean, sd, and degrees of freedom)
              seedling_size_parameters = seedling_size_function(abLarge = seed_produced$abLarge_unscaled,
                                                                prevwinterT = prevwinterT,
                                                                TSFcont = TSF,
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
              
              if(any(is.na(seed_produced$size_unscaled[which(seed_produced$goCont == 1)]), na.rm = T)) seed_produced$size_unscaled[which(seed_produced$goCont == 1 & seed_produced$size_unscaled < 0)] = 0
              
              
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
                                                                TSFcont = TSF,
                                                                fallR = fallR,
                                                                summerT = summerT,
                                                                year = year_obs,
                                                                population = population),
                                       size = 1)
          
          
          ### GROWTH ###
          
          data_droso$sizeNext = NA
          
          # Growth parameters (mean, sd, and degrees of freedom)
          growth_parameters = growth_function(size = data_droso$size_unscaled, 
                                              fallR = fallR, 
                                              abLarge = data_droso$abLarge_unscaled,
                                              TSFcont = TSF,
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
                                                                                       left = 0, 
                                                                                       right = max(data_droso$size_unscaled)))
          
          
          # Assign max size in current dataset to individual with infinite size
          if(any(data_droso$sizeNext[which(data_droso$survival == 1)] == Inf, na.rm = T)) data_droso$sizeNext[which(data_droso$survival == 1 & data_droso$sizeNext == Inf)] = max(data_droso$sizeNext[which(data_droso$survival == 1 & data_droso$sizeNext != Inf)], na.rm = T)
          
          if(any(is.na(data_droso$sizeNext[which(data_droso$survival == 1)]), na.rm = T)) data_droso$sizeNext[which(data_droso$survival == 1 & data_droso$sizeNext < 0)] = 0
          
          # Format dataset
          data_droso = data_droso[which(data_droso$survival == 1), colnames(data_droso)[-which(colnames(data_droso) %in% c("rownb"))]]
        } 
      }
    }
    
    
    if(nrow(data_SB) > 0){
      
      ### SEEDBANK ###
      
      # Seeds germinating from the seedbank 
      
      data_SB$outSB = rbinom(n = nrow(data_SB),
                             prob = outSB_function(TSF = TSF),
                             size = 1)
      
      
      # Assign a size to the germinated seeds, add them to the dataset 
      # of aboveground individuals and remove them from the seedbank data
      
      # Seedling size parameters (mean, sd, and degrees of freedom)
      seedling_size_parameters = seedling_size_function(abLarge = data_SB$abLarge_unscaled,
                                                        prevwinterT = prevwinterT,
                                                        TSFcont = TSF,
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
        if(any(is.na(data_SB$size_unscaled[which(data_SB$outSB == 1)]), na.rm = T)) data_SB$size_unscaled[which(data_SB$outSB == 1 & data_SB$size_unscaled < 0)] = 0
        
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
                              prob = staySB_function(TSF = TSF),
                              size = 1)
      
      # Keep only seeds staying in the seedbank
      data_SB = data_SB[which(data_SB$staySB == 1), ]
    }
    
    # Merge seeds going to SB to seedbank data
    seed_produced$staySB = seed_produced$goSB
    
    if(nrow(seed_produced[which(seed_produced$staySB == 1), ]) > 0){
      
      data_SB = rbind(data_SB, seed_produced[which(seed_produced$staySB == 1), ])
    }
    
    
    # Format seedbank data to match other datasets
    if(nrow(data_SB) > 0) data_SB = data_SB[, colnames(data_SB)[-which(colnames(data_SB) %in% c("staySB", "outSB", "rownb", "goCont", "goSB"))]]
    
    seedbank_size_vec[i] = nrow(data_SB)
    
    if(nrow(seed_produced) > 0){
      
      # Keep only seeds germinating 
      seed_produced = seed_produced[which(seed_produced$goCont == 1), colnames(seed_produced)[-which(colnames(seed_produced) %in% c("goCont", "goSB", "staySB", "outSB", "rownb"))]]
      
      # Cap the number of recruits if needed
      if(!is.null(recruitCap) & nrow(seed_produced) > 0){
        
        # Get number of seedlings per quadrat
        quadratsAboveMaxSeedlings = aggregate(ID ~ quadratID, 
                                              data = seed_produced,
                                              FUN = function(x) length(x))
        
        # Keep quadrats where the number of seedlings is above the threshold
        quadratsAboveMaxSeedlings = quadratsAboveMaxSeedlings[which(quadratsAboveMaxSeedlings$ID > recruitCap$ID[which(recruitCap$TSFcont == TSFcat)]), ]
        
        # If quadrats are above the threshold, sample the seedlings that
        # will be kept
        if(nrow(quadratsAboveMaxSeedlings) > 0){
          
          for(quadrat in quadratsAboveMaxSeedlings$quadratID){
            
            recruitsKept_ID = sample(seq(1, nrow(seed_produced[which(seed_produced$quadratID == quadrat), ])), 
                                     size = recruitCap$ID[which(recruitCap$TSFcont == TSFcat)], replace = F)
            
            recruitsKept = seed_produced[which(seed_produced$quadratID == quadrat)[recruitsKept_ID], ]
            
            seed_produced = seed_produced[-which(seed_produced$quadratID == quadrat), ]
            seed_produced = rbind(seed_produced, recruitsKept)
          }
        }
      }
      
      
      # Adding new seedlings to the population
      
      if(any(is.na(data_droso$ID))) data_droso = rbind(data_droso[which(is.na(data_droso$ID)), ], seed_produced)
      
      else data_droso = rbind(data_droso, seed_produced)
    }
    
    if(nrow(data_droso) > 0){
      
      data_droso$TSFcont_unscaled[which(is.na(data_droso$TSFcont_unscaled))] = TSF
      data_droso$time = (year + 1)
    }
    
    if(is.infinite(log(nrow(data_droso)/nrow(sim_data[which(sim_data$time_sim == i-1), ])))) break
    
    log_meanChangeAboveground[i] = log(nrow(data_droso)/nrow(sim_data[which(sim_data$time_sim == i-1), ])) # Calculate mean change in aboveground population abundance
    log_lambda[i] = log((nrow(data_droso) + nrow(data_SB))/(nrow(sim_data[which(sim_data$time_sim == i-1), ]) + seedbank_size_vec[i-1])) # Calculate log lambda (with seedbank)
    
    
    ### SAVE DATA ###
    
    time_sim = i # Update timestep
    
    # Merge yearly individual data with full individual data
    if(nrow(data_droso) > 0) data_droso = cbind(data_droso, time_sim)
    sim_data = rbind(sim_data, data_droso[, colnames(data_droso)[which(colnames(data_droso) %in% colnames(sim_data))]])
    
    # Format yearly data
    data_droso = data_droso[, -ncol(data_droso)]
    
    
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
              log_lambda = log_lambda))
}


## 3.2. Simulation function ----
# -------------------------

ibm_natural = function(n_sim = 1000,               # Number of simulations
                       n_years = 50,               # Number of years per simulation
                       years_RE,                   # Observed years available for random year effect
                       data_initial,               # Initial dataset
                       seedbank_size = 10000,      # Initial seedbank size
                       recruitCap,                 # Maximum number of recruits per quadrat
                       max_nbFlowers,
                       population){                # Population ID
  
  # Initialization - Prepare storing objects
  
  log_meanChangeAboveground_array = array(NA, dim = c(n_sim, n_years))
  log_lambda_array = array(NA, dim = c(n_sim, n_years))
  pop_size_list = vector(mode = "list", length = n_sim)
  pop_data_list = vector(mode = "list", length = n_sim)
  pop_density_list = vector(mode = "list", length = n_sim)
  
  ibm_data = NULL
  
  # Calculate number of flowers
  data_initial$nbFlowers = data_initial$fs * data_initial$fps
  
  # Format initial data
  data_initial = data_initial[, c("site", "quadratID", "ID", "TSFcont_unscaled", "size_unscaled", "sizeNext", 
                                  "fl", "nbFlowers", 
                                  "surv", "time", "abLarge_unscaled")]
  
  colnames(data_initial) = c("site", "quadratID", "ID", "TSFcont_unscaled", "size_unscaled", "sizeNext", 
                             "flowering", "nbFlowers", 
                             "survival", "time", "abLarge_unscaled")
  
  # Prepare initial seedbank data
  seedbank_initial = data.frame(site = population, quadratID = NA, ID = paste("Seed", seq(1, seedbank_size), sep = "_"),
                                TSFcont_unscaled = NA, size_unscaled = NA, sizeNext = NA, 
                                flowering = NA, nbFlowers = NA,
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
                    years_RE = years_RE,
                    seedbank_size = seedbank_size,
                    recruitCap = recruitCap,
                    max_nbFlowers = max_nbFlowers,
                    population = population,
                    ibm_data = ibm_data,
                    data_initial = data_initial,
                    seedbank_initial = seedbank_initial),
            silent = TRUE)

      if(!is(ibm_sim_result, 'try-error')) break
    }
    
    # Fill in result objects
    log_meanChangeAboveground_array[sim, ] = ibm_sim_result$log_meanChangeAboveground
    log_lambda_array[sim, ] = ibm_sim_result$log_lambda
    pop_size_list[[sim]] = ibm_sim_result$pop_size
    pop_data_list[[sim]] = ibm_sim_result$pop_data
    pop_density_list[[sim]] = ibm_sim_result$pop_density
    
  }
  
  return(list(pop_data = pop_data_list,
              pop_size = pop_size_list,
              pop_density = pop_density_list,
              log_meanChangeAboveground = log_meanChangeAboveground_array,
              log_lambda = log_lambda_array))
}




###########################################################################
#
# 4. Projections ----
#
###########################################################################

n_sim = 5

## 4.1. SierraCarboneraY5 ----
# -----------------------

population = "SierraCarboneraY5"
droso_pop = droso_natural_full[which(droso_natural_full$site == population), ]

recruitCap = aggregate(ID ~ quadratID + time + TSFcont_unscaled,
                       data = droso_pop[which(droso_pop$stage == "SD"), ], 
                       FUN = function(x) length(x))

colnames(recruitCap)[which(colnames(recruitCap) == "TSFcont_unscaled")] = "TSFcont"

recruitCap$TSFcont[which(recruitCap$TSFcont > 4)] = 4

recruitCap = aggregate(ID ~ TSFcont,
                       data = recruitCap, FUN = max)

recruitCap = rbind(recruitCap, 
                   data.frame(TSFcont = 0,
                              ID = round(max(recruitCap$ID) * 1.5)))

recruitCap = recruitCap[order(recruitCap$TSFcont), ]


max_nbFlowers = aggregate(fs * fps ~ TSFcont_unscaled,
                          data = droso_pop, FUN = max)
colnames(max_nbFlowers) = c("TSFcont", "nbFlowers")

max_nbFlowers$TSFcont[which(max_nbFlowers$TSFcont > 4)] = 4

max_nbFlowers = aggregate(nbFlowers ~ TSFcont,
                          data = max_nbFlowers, FUN = max)


n_years = length(unique(droso_pop$time))

ibm_sierracarboneray5_validation = ibm_natural(n_sim = n_sim,
                                               n_years = n_years,
                                               years_RE = unique(droso_pop$time),
                                               data_initial = droso_pop,
                                               recruitCap = recruitCap,
                                               max_nbFlowers = max_nbFlowers,
                                               population = population)

# Saving results
save(ibm_sierracarboneray5_validation, 
     file = paste0("Output/Projections/IBM_Natural_Validation_", 
                   population, ".RData"))


## 4.2. SierraRetinY5 ----
# -------------------

population = "SierraRetinY5"
droso_pop = droso_natural_full[which(droso_natural_full$site == population), ]

recruitCap = aggregate(ID ~ quadratID + time + TSFcont_unscaled,
                       data = droso_pop[which(droso_pop$stage == "SD"), ], 
                       FUN = function(x) length(x))
colnames(recruitCap)[which(colnames(recruitCap) == "TSFcont_unscaled")] = "TSFcont"

recruitCap$TSFcont[which(recruitCap$TSFcont > 4)] = 4

recruitCap = aggregate(ID ~ TSFcont,
                       data = recruitCap, FUN = max)

recruitCap$ID = round(recruitCap$ID)

recruitCap = rbind(recruitCap, 
                   data.frame(TSFcont = c(0, 2),
                              ID = c(round(max(recruitCap$ID) * 1.5), round(mean(recruitCap$ID[which(recruitCap$TSFcont != 1)])))))

recruitCap = recruitCap[order(recruitCap$TSFcont), ]


max_nbFlowers = aggregate(fs * fps ~ TSFcont_unscaled,
                          data = droso_pop, FUN = max)
colnames(max_nbFlowers) = c("TSFcont", "nbFlowers")

max_nbFlowers$TSFcont[which(max_nbFlowers$TSFcont > 4)] = 4

max_nbFlowers = aggregate(nbFlowers ~ TSFcont,
                          data = max_nbFlowers, FUN = max)

max_nbFlowers = rbind(max_nbFlowers, 
                      data.frame(TSFcont = 3,
                                 nbFlowers = round(mean(max_nbFlowers$nbFlowers))))

max_nbFlowers = max_nbFlowers[order(max_nbFlowers$TSFcont), ]


n_years = length(unique(droso_pop$time))


ibm_sierraretiny5_validation = ibm_natural(n_sim = n_sim,
                                           n_years = n_years,
                                           years_RE = unique(droso_pop$time),
                                           data_initial = droso_pop,
                                           recruitCap = recruitCap,
                                           max_nbFlowers = max_nbFlowers,
                                           population = population)

# Saving results
save(ibm_sierraretiny5_validation, 
     file = paste0("Output/Projections/IBM_Natural_Validation_", 
                   population, ".RData"))


## 4.3. Vertedero ----
# ---------------

population = "Vertedero"
droso_pop = droso_natural_full[which(droso_natural_full$site == population), ]

recruitCap = aggregate(ID ~ quadratID + time + TSFcont_unscaled,
                       data = droso_pop[which(droso_pop$stage == "SD"), ], FUN = function(x) length(x))
colnames(recruitCap)[which(colnames(recruitCap) == "TSFcont_unscaled")] = "TSFcont"

recruitCap$TSFcont[which(recruitCap$TSFcont > 4)] = 4

recruitCap = aggregate(ID ~ TSFcont,
                       data = recruitCap, FUN = max)

recruitCap$ID = round(recruitCap$ID)

recruitCap = rbind(recruitCap, 
                   data.frame(TSFcont = c(0, 1, 2),
                              ID = c(round(max(recruitCap$ID) * 1.5), max(recruitCap$ID), round(mean(recruitCap$ID[which(recruitCap$TSFcont != 1)])))))

recruitCap = recruitCap[order(recruitCap$TSFcont), ]


max_nbFlowers = aggregate(fs * fps ~ TSFcont_unscaled,
                          data = droso_pop, FUN = max)
colnames(max_nbFlowers) = c("TSFcont", "nbFlowers")

max_nbFlowers$TSFcont[which(max_nbFlowers$TSFcont > 4)] = 4

max_nbFlowers = aggregate(nbFlowers ~ TSFcont,
                          data = max_nbFlowers, FUN = max)


n_years = length(unique(droso_pop$time))


ibm_vertedero_validation = ibm_natural(n_sim = n_sim,
                                       n_years = n_years,
                                       years_RE = unique(droso_pop$time),
                                       data_initial = droso_pop,
                                       recruitCap = recruitCap,
                                       max_nbFlowers = max_nbFlowers,
                                       population = population)

# Saving results
save(ibm_vertedero_validation, 
     file = paste0("Output/Projections/IBM_Natural_Validation_", 
                   population, ".RData"))




###########################################################################
#
# 5. Observed and projected population metrics ----
#
###########################################################################

# Results data frame

ibm_results = expand.grid(year = unique(droso_natural_full$time),
                          simulation = seq(1, n_sim), 
                          dormancy = "Natural",
                          source = "Projected",
                          population = c("SierraCarboneraY5", "SierraRetinY5"
                                         ,
                                         "Vertedero"
                                         ))

ibm_results = rbind(ibm_results, expand.grid(year = unique(droso_natural_full$time),
                                             simulation = NA,
                                             dormancy = "Natural",
                                             source = "Observed",
                                             population = c("SierraCarboneraY5", "SierraRetinY5"
                                                            ,
                                                            "Vertedero"
                                                            )))

ibm_results$pop_size = NA
ibm_results$pop_size_repro = NA


## 5.1. Functions to calculate population size ----
# --------------------------------------------

# Get population size

addPopSize = function(dataset, population){
  
  demo_data = dataset$pop_data[[length(dataset$pop_data)]]
  n_years = length(unique(droso_natural_full$time[which(droso_natural_full$site == population)]))

  pop_size = lapply(dataset$pop_size, 
                    FUN = function(x){
                      
                      sim = as.numeric(unique(x[3]))
                      
                      max_time = max(x[1])
                      
                      time_sim_to_add = which(!(seq(1, n_years) %in% x[1]$time_sim))
                      
                      if(length(time_sim_to_add) > 0){
                        
                        x = rbind(data.frame(time_sim = time_sim_to_add,
                                             ID = NA,
                                             run = sim),
                                  x)
                        
                        x = x[order(x[1]$time_sim), ]
                      }
                      
                      return(x[2])
                    })
  
  pop_size_length = unlist(lapply(pop_size, 
                                  FUN = function(x) nrow(x)))
  
  return(as.numeric(unlist(pop_size)))
}


## 5.2. Observed population metrics ----
# ---------------------------------

pop_size = aggregate(ID ~ time + site,
                     data = droso_natural_full, 
                     FUN = function(x) length(x))


for(pop in unique(pop_size$site)){
  
  for(yr in unique(pop_size$time[which(pop_size$site == pop)])){
    
    ibm_results$pop_size[which(ibm_results$source == "Observed" &
                               ibm_results$population == pop &
                               ibm_results$year == yr)] = pop_size$ID[which(pop_size$site == pop &
                                                                            pop_size$time == yr)]
    
  }
}

ibm_results$log_meanChangeAboveground = NA

for(pop in unique(ibm_results$population)){
  
  for(yr in unique(ibm_results$year[which(ibm_results$population == pop)])){
    
    if(yr == unique(ibm_results$year[which(ibm_results$population == pop)])[1]) next
    
    ibm_results$log_meanChangeAboveground[which(ibm_results$year == yr & ibm_results$population == pop & ibm_results$source == "Observed")] = log(nrow(droso_natural_full[which(droso_natural_full$site == pop & droso_natural_full$time == yr), ])/nrow(droso_natural_full[which(droso_natural_full$site == pop & droso_natural_full$time == (yr-1)), ]))
  }
}


## 6.3. Projected population metrics ----
# ----------------------------------

ibm_results$pop_size[which(ibm_results$source == "Projected" &
                           ibm_results$population == "SierraCarboneraY5" &
                           ibm_results$year %in% unique(droso_natural_full$time[which(droso_natural_full$site == "SierraCarboneraY5")]))] = addPopSize(ibm_sierracarboneray5_validation, "SierraCarboneraY5")


ibm_results$pop_size[which(ibm_results$source == "Projected" &
                           ibm_results$population == "SierraRetinY5" &
                           ibm_results$year %in% unique(droso_natural_full$time[which(droso_natural_full$site == "SierraRetinY5")]))] = addPopSize(ibm_sierraretiny5_validation, "SierraRetinY5")


ibm_results$pop_size[which(ibm_results$source == "Projected" &
                           ibm_results$population == "Vertedero" &
                           ibm_results$year %in% unique(droso_natural_full$time[which(droso_natural_full$site == "Vertedero")]))] = addPopSize(ibm_vertedero_validation, "Vertedero")


ibm_results$log_lambda = NA
ibm_results$log_meanChangeAboveground[which(ibm_results$population == "SierraCarboneraY5" &
                                            ibm_results$source == "Projected" &
                                            ibm_results$year %in% unique(droso_natural_full$time[which(droso_natural_full$site == "SierraCarboneraY5")]))] =
  c(t(ibm_sierracarboneray5_validation$log_meanChangeAboveground))
ibm_results$log_lambda[which(ibm_results$population == "SierraCarboneraY5" &
                             ibm_results$source == "Projected" &
                             ibm_results$year %in% unique(droso_natural_full$time[which(droso_natural_full$site == "SierraCarboneraY5")]))] =
  c(t(ibm_sierracarboneray5_validation$log_lambda))

ibm_results$log_meanChangeAboveground[which(ibm_results$population == "SierraRetinY5" &
                                            ibm_results$source == "Projected" &
                                            ibm_results$year %in% unique(droso_natural_full$time[which(droso_natural_full$site == "SierraRetinY5")]))] =
  c(t(ibm_sierraretiny5_validation$log_meanChangeAboveground))
ibm_results$log_lambda[which(ibm_results$population == "SierraRetinY5" &
                             ibm_results$source == "Projected" &
                             ibm_results$year %in% unique(droso_natural_full$time[which(droso_natural_full$site == "SierraRetinY5")]))] =
  c(t(ibm_sierraretiny5_validation$log_lambda))

ibm_results$log_meanChangeAboveground[which(ibm_results$population == "Vertedero" &
                                            ibm_results$source == "Projected" &
                                            ibm_results$year %in% unique(droso_natural_full$time[which(droso_natural_full$site == "Vertedero")]))] =
  c(t(ibm_vertedero_validation$log_meanChangeAboveground))
ibm_results$log_lambda[which(ibm_results$population == "Vertedero" &
                             ibm_results$source == "Projected" &
                             ibm_results$year %in% unique(droso_natural_full$time[which(droso_natural_full$site == "Vertedero")]))] =
  c(t(ibm_vertedero_validation$log_lambda))


ibm_results_natural = ibm_results[which(!is.na(ibm_results$pop_size)), ]


size_dist_natural = ibm_sierracarboneray5_validation$pop_data[[1]][0, ]

for(i in 1:n_sim){
  
  size_dist_natural = rbind(size_dist_natural, ibm_sierracarboneray5_validation$pop_data[[i]])
}

for(i in 1:n_sim){
  
  size_dist_natural = rbind(size_dist_natural, ibm_sierraretiny5_validation$pop_data[[i]])
}

for(i in 1:n_sim){
  
  size_dist_natural = rbind(size_dist_natural, ibm_vertedero_validation$pop_data[[i]])
}

size_dist_natural = size_dist_natural[, which(colnames(size_dist_natural) %in% colnames(droso_natural_full))]
size_dist_natural$srce = "Projected"

size_dist_natural = rbind(size_dist_natural,
                  cbind(droso_natural_full[, which(colnames(droso_natural_full) %in% colnames(size_dist_natural))], data.frame(srce = "Observed")))


write.csv(ibm_results_natural, file = "Output/Projections/ValidationResults_Natural.csv", row.names = F)
write.csv(size_dist_natural, file = "Output/Projections/SizeDistribution_Natural.csv", row.names = F)