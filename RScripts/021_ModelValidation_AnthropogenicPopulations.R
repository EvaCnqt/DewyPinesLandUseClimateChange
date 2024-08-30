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


## 1.2. Loading libraries ----
# -----------------------

library(crch)


## 1.3. Loading data ----
# ------------------

# Dewy-pine data
droso_anthropogenic = read.csv("Data/droso_anthropogenic.csv")
droso_anthropogenic_full = read.csv("Data/droso_anthropogenic_full.csv")

droso_seedbank = read.csv("Data/droso_seedbank_anthropogenic.csv")

seeds_per_flower = 9.8


# Correction factors (seed survival)
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
  survival = predict(surv_anthropogenic, newdata = data.frame(size = size_scaled,
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
  growth_mean = predict(growth_anthropogenic, newdata = data.frame(size = size_scaled,
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
  seedling_size_mean = as.numeric(predict(seedlingSize_anthropogenic, newdata = data.frame(abLarge = abLarge_scaled,
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
  flowering = predict(flowering_anthropogenic, newdata = data.frame(size = size_scaled,
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
  nbFlowers = predict(nbFlow_anthropogenic, newdata = data.frame(size = size_scaled,
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
# 3. Individual-based model function ----
#
###########################################################################

## 3.1. Timestep projection function ----
# ----------------------------------

ibm_sim = function(n_sim,             # Number of simulations
                   n_years,           # Number of years per simulation
                   sim,               # Current simulation number
                   years_RE,          # Observed years available for random year effect
                   data_initial,      # Initial dataset
                   seedbank_size,     # Initial seedbank size
                   recruitCap,        # Maximum number of recruits per quadrat
                   population,        # Population ID
                   ibm_data,          # Previously stored projection data
                   seedbank_initial){ # Initial seedbank data
  
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
  data_droso$abLarge_unscaled = apply(data_droso, 1, 
                                      function(x) nrow(data_droso[which(data_droso$quadratID == x[2] & data_droso$size_unscaled > 4.5), ]))
  pop_density[[1]] = aggregate(abLarge_unscaled ~ quadratID, 
                               data = data_droso, mean, 
                               na.rm = T)$abLarge_unscaled
  
  # Seedbank seeds data
  data_SB = seedbank_initial
  
  # Assign seeds in initial seedbank to quadrats
  data_SB$quadratID = sample(unique(data_droso$quadratID), size = seedbank_size, 
                             replace = T, 
                             prob = as.numeric(table(data_droso$quadratID)/sum(table(data_droso$quadratID))))
  
  # Add density data to seedbank
  data_SB$abLarge_unscaled = apply(data_SB, 1, 
                                   function(x) nrow(data_droso[which(data_droso$quadratID == x[2] & data_droso$size_unscaled > 4.5), ]))
  
  # Highest seed ID (format Seed_XXX) to give names to new seeds 
  max_seed_ID = max(as.numeric(unlist(lapply(strsplit(data_SB$ID, split = "_"), 
                                             function(x) x[2]))))
  
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
    
    summerT = unique(droso_anthropogenic_full$summerT_unscaled[which(droso_anthropogenic_full$site == population &
                                                                   droso_anthropogenic_full$time == year)])
    prevwinterT = unique(droso_anthropogenic_full$prevwinterT_unscaled[which(droso_anthropogenic_full$site == population &
                                                                         droso_anthropogenic_full$time == year)])
    
    fallR = unique(droso_anthropogenic_full$fallR_unscaled[which(droso_anthropogenic_full$site == population &
                                                             droso_anthropogenic_full$time == year)])
    prevfallR = unique(droso_anthropogenic_full$prevfallR_unscaled[which(droso_anthropogenic_full$site == population &
                                                                     droso_anthropogenic_full$time == year)])
    prevwinterR = unique(droso_anthropogenic_full$prevwinterR_unscaled[which(droso_anthropogenic_full$site == population &
                                                                         droso_anthropogenic_full$time == year)])

    
    ### CORRECTION FACTORS AND NUMBER OF SQUARES ###
    
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

    #
    
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
    
    log_meanChangeAboveground[i] = log(nrow(data_droso)/nrow(sim_data[which(sim_data$time_sim == i-1), ])) # Calculate mean change in aboveground population abundance
    log_lambda[i] = log((nrow(data_droso) + nrow(data_SB))/(nrow(sim_data[which(sim_data$time_sim == i-1), ]) + seedbank_size_vec[i-1])) # Calculate log lambda (with seedbank)
    
    
    ### SAVE DATA ###
    
    time_sim = i # Update timestep
    
    # Merge yearly individual data with full individual data
    if(nrow(data_droso) > 0){
      
      data_droso = cbind(data_droso, time_sim)
      data_droso$time = (year + 1)
    
    }
      
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

ibm_anthropogenic = function(n_sim = 1000,                # Number of simulations
                             n_years = 50,                # Number of years per simulation
                             years_RE,                    # Observed years available for random year effect
                             data_initial,                # Initial dataset
                             seedbank_size = 3000,        # Initial seedbank size
                             recruitCap,                  # Maximum number of recruits per quadrat
                             population){                 # Population ID
  
  
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
                    years_RE = years_RE,
                    seedbank_size = seedbank_size,
                    recruitCap = recruitCap,
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

n_sim = 500

## 4.1. Retin ----
# -----------

population = "Retin"
droso_pop = droso_anthropogenic_full[which(droso_anthropogenic_full$site == population), ]

recruitCap = NULL

n_years = length(unique(droso_pop$time[which(droso_pop$site == population)]))


ibm_retin_validation = ibm_anthropogenic(n_sim = n_sim, 
                                         n_years = n_years, 
                                         years_RE = unique(droso_pop$time[which(droso_pop$site == population)]),
                                         data_initial = droso_pop, 
                                         recruitCap = recruitCap,
                                         population = population)

# Saving results
save(ibm_retin_validation,
     file = paste0("Output/Projections/IBM_Anthropogenic_Validation_",
                   population, ".RData"))


## 4.2. Prisoneros ----
# ----------------

population = "Prisoneros"
droso_pop = droso_anthropogenic_full[which(droso_anthropogenic_full$site == population), ]

recruitCap = NULL

n_years = length(unique(droso_pop$time[which(droso_pop$site == population)]))


ibm_prisoneros_validation = ibm_anthropogenic(n_sim = n_sim, 
                                              n_years = n_years, 
                                              years_RE = unique(droso_pop$time[which(droso_pop$site == population)]),
                                              data_initial = droso_pop, 
                                              recruitCap = recruitCap,
                                              population = population)

# Saving results
save(ibm_prisoneros_validation,
     file = paste0("Output/Projections/IBM_Anthropogenic_Validation_",
                   population, ".RData"))


## 4.3. Bujeo ----
# -----------

population = "Bujeo"
droso_pop = droso_anthropogenic_full[which(droso_anthropogenic_full$site == population), ]

recruitCap = aggregate(ID ~ quadratID + time,
                       data = droso_anthropogenic[which(droso_anthropogenic$site == population &
                                                      droso_anthropogenic$stage == "SD"), ], FUN = function(x) length(x))

recruitCap = round(max(recruitCap$ID))

n_years = length(unique(droso_pop$time[which(droso_pop$site == population)]))


ibm_bujeo_validation = ibm_anthropogenic(n_sim = n_sim, 
                                         n_years = n_years, 
                                         years_RE = unique(droso_pop$time[which(droso_pop$site == population)]),
                                         data_initial = droso_pop, 
                                         recruitCap = recruitCap,
                                         population = population)

# Saving results
save(ibm_bujeo_validation,
     file = paste0("Output/Projections/IBM_Anthropogenic_Validation_",
                   population, ".RData"))


## 4.4. MonteraTorero ----
# -------------------

population = "MonteraTorero"
droso_pop = droso_anthropogenic_full[which(droso_anthropogenic_full$site == population), ]

recruitCap = NULL

n_years = length(unique(droso_pop$time[which(droso_pop$site == population)]))


ibm_monteratorero_validation = ibm_anthropogenic(n_sim = n_sim, 
                                                 n_years = n_years, 
                                                 years_RE = unique(droso_pop$time[which(droso_pop$site == population)]),
                                                 data_initial = droso_pop, 
                                                 recruitCap = recruitCap,
                                                 population = population)

# Saving results
save(ibm_monteratorero_validation,
     file = paste0("Output/Projections/IBM_Anthropogenic_Validation_",
                   population, ".RData"))


## 4.5. SCarbDist ----
# ---------------

population = "SCarbDist"
droso_pop = droso_anthropogenic_full[which(droso_anthropogenic_full$site == population), ]

recruitCap = aggregate(ID ~ quadratID + time,
                       data = droso_anthropogenic[which(droso_anthropogenic$site == population &
                                                      droso_anthropogenic$stage == "SD"), ], FUN = function(x) length(x))

recruitCap = round(max(recruitCap$ID))

n_years = length(unique(droso_pop$time[which(droso_pop$site == population)]))


ibm_scarbdist_validation = ibm_anthropogenic(n_sim = n_sim, 
                                             n_years = n_years, 
                                             years_RE = unique(droso_pop$time[which(droso_pop$site == population)]),
                                             data_initial = droso_pop, 
                                             recruitCap = recruitCap,
                                             population = population)

# Saving results
save(ibm_scarbdist_validation,
     file = paste0("Output/Projections/IBM_Anthropogenic_Validation_",
                   population, ".RData"))




###########################################################################
#
# 5. Observed and projected population metrics ----
#
###########################################################################

# Results data frame

ibm_results = expand.grid(year = unique(droso_anthropogenic_full$time),
                          simulation = seq(1, n_sim), 
                          dormancy = "Anthropogenic",
                          source = "Projected",
                          population = c("Retin", "Prisoneros", "Bujeo",
                                         "MonteraTorero", "SCarbDist"))

ibm_results = rbind(ibm_results, expand.grid(year = unique(droso_anthropogenic_full$time),
                                             simulation = NA,
                                             dormancy = "Anthropogenic",
                                             source = "Observed",
                                             population = c("Retin", "Prisoneros", 
                                                            "Bujeo", "MonteraTorero", 
                                                            "SCarbDist")))

ibm_results$pop_size = NA


## 5.1. Functions to calculate population size ----
# --------------------------------------------

# Get population size

addPopSize = function(dataset){
  
  demo_data = dataset$pop_data[[length(dataset$pop_data)]]
  n_years = length(unique(demo_data$time[which(!is.na(demo_data$time))]))
  
  pop_size = lapply(dataset$pop_size, 
                    FUN = function(x){
                      
                      max_time = max(x[1])
                      
                      if(max_time < n_years){
                        
                        sim = unique(x[3])
                        complement_df = data.frame(time_sim = seq(max_time + 1, n_years),
                                                   ID = NA,
                                                   run = sim)
                        x = rbind(x, complement_df)
                      }
                      
                      return(x[2])
                    })
  
  pop_size_length = unlist(lapply(pop_size, 
                                  FUN = function(x) nrow(x)))
  
  if(any(pop_size_length > n_years)){
    
    pop_size[which(pop_size_length > n_years)] = lapply(pop_size[which(pop_size_length > n_years)],
                                                   FUN = function(x){
                                                     
                                                     return(x[-c(n_years+1:nrow(x)), ])
                                                   })
  }
  
  if(any(pop_size_length < n_years)){
    
    pop_size[which(pop_size_length < n_years)] = lapply(pop_size[which(pop_size_length < n_years)],
                                                              FUN = function(x){
                                                                
                                                                x = rbind(x, rep(NA, n_years-nrow(x)))
                                                                
                                                                return(x)
                                                              })
  }
  
  return(as.numeric(unlist(pop_size)))
}


## 5.2. Observed population metrics ----
# ---------------------------------

pop_size = aggregate(ID ~ time + site,
                     data = droso_anthropogenic_full, 
                     FUN = function(x) length(x))


for(pop in unique(pop_size$site)){
  
  for(yr in unique(pop_size$time[which(pop_size$site == pop)])){
    
    ibm_results$pop_size[which(ibm_results$source == "Observed" &
                               ibm_results$population == pop &
                               ibm_results$year == yr)] = pop_size$ID[which(pop_size$site == pop &
                                                                            pop_size$time == yr)]
    
  }
}


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
    
    ibm_results$log_meanChangeAboveground[which(ibm_results$year == yr & ibm_results$population == pop & ibm_results$source == "Observed")] = log(nrow(droso_anthropogenic_full[which(droso_anthropogenic_full$site == pop & droso_anthropogenic_full$time == yr), ])/nrow(droso_anthropogenic_full[which(droso_anthropogenic_full$site == pop & droso_anthropogenic_full$time == (yr-1)), ]))
  }
}

ibm_results = ibm_results[-which(ibm_results$population != "Retin" &
                                ibm_results$year < 2016), ]


## 5.3. Projected population metrics ----
# ----------------------------------

ibm_results$pop_size[which(ibm_results$source == "Projected" &
                           ibm_results$population == "Retin" &
                           ibm_results$year %in% unique(droso_anthropogenic_full$time[which(droso_anthropogenic_full$site == "Retin")]))] = addPopSize(ibm_retin_validation)


ibm_results$pop_size[which(ibm_results$source == "Projected" &
                           ibm_results$population == "Prisoneros" &
                           ibm_results$year %in% unique(droso_anthropogenic_full$time[which(droso_anthropogenic_full$site == "Prisoneros")]))] = addPopSize(ibm_prisoneros_validation)


ibm_results$pop_size[which(ibm_results$source == "Projected" &
                           ibm_results$population == "Bujeo" &
                           ibm_results$year %in% unique(droso_anthropogenic_full$time[which(droso_anthropogenic_full$site == "Bujeo")]))] = addPopSize(ibm_bujeo_validation)


ibm_results$pop_size[which(ibm_results$source == "Projected" &
                           ibm_results$population == "MonteraTorero" &
                           ibm_results$year %in% unique(droso_anthropogenic_full$time[which(droso_anthropogenic_full$site == "MonteraTorero")]))] = addPopSize(ibm_monteratorero_validation)


ibm_results$pop_size[which(ibm_results$source == "Projected" &
                           ibm_results$population == "SCarbDist" &
                           ibm_results$year %in% unique(droso_anthropogenic_full$time[which(droso_anthropogenic_full$site == "SCarbDist")]))] = addPopSize(ibm_scarbdist_validation)


ibm_results$log_lambda = NA

ibm_results$log_meanChangeAboveground[which(ibm_results$population == "Retin" &
                                            ibm_results$source == "Projected" &
                                            ibm_results$year %in% unique(droso_anthropogenic_full$time[which(droso_anthropogenic_full$site == "Retin")]))] = 
  c(t(ibm_retin_validation$log_meanChangeAboveground))
ibm_results$log_lambda[which(ibm_results$population == "Retin" &
                               ibm_results$source == "Projected" &
                               ibm_results$year %in% unique(droso_anthropogenic_full$time[which(droso_anthropogenic_full$site == "Retin")]))] = 
  c(t(ibm_retin_validation$log_lambda))

ibm_results$log_meanChangeAboveground[which(ibm_results$population == "Prisoneros" &
                                              ibm_results$source == "Projected" &
                                              ibm_results$year %in% unique(droso_anthropogenic_full$time[which(droso_anthropogenic_full$site == "Prisoneros")]))] = 
  c(t(ibm_prisoneros_validation$log_meanChangeAboveground))
ibm_results$log_lambda[which(ibm_results$population == "Prisoneros" &
                               ibm_results$source == "Projected" &
                               ibm_results$year %in% unique(droso_anthropogenic_full$time[which(droso_anthropogenic_full$site == "Prisoneros")]))] = 
  c(t(ibm_prisoneros_validation$log_lambda))

ibm_results$log_meanChangeAboveground[which(ibm_results$population == "Bujeo" &
                                              ibm_results$source == "Projected" &
                                              ibm_results$year %in% unique(droso_anthropogenic_full$time[which(droso_anthropogenic_full$site == "Bujeo")]))] = 
  c(t(ibm_bujeo_validation$log_meanChangeAboveground))
ibm_results$log_lambda[which(ibm_results$population == "Bujeo" &
                               ibm_results$source == "Projected" &
                               ibm_results$year %in% unique(droso_anthropogenic_full$time[which(droso_anthropogenic_full$site == "Bujeo")]))] = 
  c(t(ibm_bujeo_validation$log_lambda))

ibm_results$log_meanChangeAboveground[which(ibm_results$population == "MonteraTorero" &
                                              ibm_results$source == "Projected" &
                                              ibm_results$year %in% unique(droso_anthropogenic_full$time[which(droso_anthropogenic_full$site == "MonteraTorero")]))] = 
  c(t(ibm_monteratorero_validation$log_meanChangeAboveground))
ibm_results$log_lambda[which(ibm_results$population == "MonteraTorero" &
                               ibm_results$source == "Projected" &
                               ibm_results$year %in% unique(droso_anthropogenic_full$time[which(droso_anthropogenic_full$site == "MonteraTorero")]))] = 
  c(t(ibm_monteratorero_validation$log_lambda))

ibm_results$log_meanChangeAboveground[which(ibm_results$population == "SCarbDist" &
                                              ibm_results$source == "Projected" &
                                              ibm_results$year %in% unique(droso_anthropogenic_full$time[which(droso_anthropogenic_full$site == "SCarbDist")]))] = 
  c(t(ibm_scarbdist_validation$log_meanChangeAboveground))
ibm_results$log_lambda[which(ibm_results$population == "SCarbDist" &
                               ibm_results$source == "Projected" &
                               ibm_results$year %in% unique(droso_anthropogenic_full$time[which(droso_anthropogenic_full$site == "SCarbDist")]))] = 
  c(t(ibm_scarbdist_validation$log_lambda))

ibm_results_anthropogenic = ibm_results[which(!is.na(ibm_results$pop_size)), ]


# Size distributions

size_dist_anthropogenic = ibm_bujeo_validation$pop_data[[1]][0, ]

for(i in 1:n_sim){
  
  size_dist_anthropogenic = rbind(size_dist_anthropogenic, ibm_bujeo_validation$pop_data[[i]])
}

for(i in 1:n_sim){
  
  size_dist_anthropogenic = rbind(size_dist_anthropogenic, ibm_monteratorero_validation$pop_data[[i]])
}

for(i in 1:n_sim){
  
  size_dist_anthropogenic = rbind(size_dist_anthropogenic, ibm_prisoneros_validation$pop_data[[i]])
}

for(i in 1:n_sim){
  
  size_dist_anthropogenic = rbind(size_dist_anthropogenic, ibm_retin_validation$pop_data[[i]])
}

for(i in 1:n_sim){
  
  size_dist_anthropogenic = rbind(size_dist_anthropogenic, ibm_scarbdist_validation$pop_data[[i]])
}

size_dist_anthropogenic = size_dist_anthropogenic[, which(colnames(size_dist_anthropogenic) %in% colnames(droso_anthropogenic_full))]
size_dist_anthropogenic$srce = "Projected"

size_dist_anthropogenic = rbind(size_dist_anthropogenic,
                                cbind(droso_anthropogenic_full[, which(colnames(droso_anthropogenic_full) %in% colnames(size_dist_anthropogenic))], data.frame(srce = "Observed")))


write.csv(ibm_results_anthropogenic, file = "Output/Projections/ValidationResults_Anthropogenic.csv", row.names = F)
write.csv(size_dist_anthropogenic, file = "Output/Projections/SizeDistribution_Anthropogenic.csv", row.names = F)