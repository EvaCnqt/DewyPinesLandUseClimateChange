############################################################################
#
# The aim of this script is to fit GAMs to estimate dewy-pine seedling
# size as a function of size, population density, temperature, 
# and rainfall in various populations.
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

library(bbmle)
library(mgcv)
library(MuMIn)


## 1.3. Loading data ----
# ------------------

droso_anthropogenic = read.csv("Data/droso_anthropogenic.csv")




###########################################################################
#
# 2. Preparing data ----
#
###########################################################################

## 2.1. Calculating mean density per quadrat for density standardization ----
# ----------------------------------------------------------------------

nbSquares = aggregate(quadratID ~ time + site, 
                      data = droso_anthropogenic,
                      FUN = function(x) length(unique(x)))
density_per_square = aggregate(abLarge_unscaled ~ quadratID + time + site, 
                               data = droso_anthropogenic, 
                               function(x) unique(x))
yearly_density_per_square = aggregate(abLarge_unscaled ~ time + site, 
                                      data = density_per_square, 
                                      function(x) sum(x))
yearly_density_per_square$abLarge_unscaled = yearly_density_per_square$abLarge_unscaled/nbSquares$quadratID


droso_anthropogenic$time = factor(droso_anthropogenic$time)
droso_anthropogenic$site = factor(droso_anthropogenic$site)


## 2.2. Subsetting data for seedling size only ----
# --------------------------------------------

droso_seedlingSize = droso_anthropogenic[which(droso_anthropogenic$stage == "SD"), ] # Seedlings




###########################################################################
#
# 3. Modelling seedling size ----
#
###########################################################################

# Null model
seedlingSize_null = gam(size_unscaled ~ s(time, bs = "re") + s(site, bs = "re"),
                        data = droso_seedlingSize, family = "scat", method = "REML",
                        gamma = 1.4, select = T)


## 3.1. Best model with climate only ----
# ----------------------------------

# Find the best lagged window:

# Previous winter temperature (mean max temperature from January to April of the census year, i.e., prior to the census)
# January = Census date - 4
# April = Census date - 1

seedlingSize_prevwinterT = update(seedlingSize_null,    ~ . + s(prevwinterT, bs = "cr", k = 3), 
                                   data = droso_seedlingSize)

# Previous fall rainfall (cumulative rainfall from September to November of the year prior to the census year)
# September = Census date - 8
# November = Census date - 6

seedlingSize_prevfallR = update(seedlingSize_null,    ~ . + s(prevfallR, bs = "cr", k = 3), 
                                 data = droso_seedlingSize)

# Previous winter rainfall (cumulative rainfall from January to April of the census year, i.e., prior to the census)
# January = Census date - 4
# April = Census date - 1

seedlingSize_prevwinterR = update(seedlingSize_null,    ~ . + s(prevwinterR, bs = "cr", k = 3), 
                                   data = droso_seedlingSize)


# Put models in two lists (temperature and rainfall)
mod_listT = list(seedlingSize_null, seedlingSize_prevwinterT)
mod_listR = list(seedlingSize_null, seedlingSize_prevfallR, seedlingSize_prevwinterR)


# Model names
names(mod_listT) = c("null", "previous_winterT")
names(mod_listR) = c("null", "previous_fallR", "previous_winterR")

model.sel(mod_listR)

model.sel(mod_listT)

AICtab(mod_listR[["null"]], mod_listT[["previous_winterT"]])

# No effect of rainfall
# Mean max. temperature in previous winter (January to April)


# Checking correlation
cor(droso_seedlingSize$prevwinterT, droso_seedlingSize$abLarge) # No
cor(droso_seedlingSize$prevwinterT, droso_seedlingSize$TSFcont) # No

cor(droso_seedlingSize$abLarge, droso_seedlingSize$TSFcont) # No


# Model selection 

RT1a = gam(size_unscaled ~ prevwinterT + 
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT1b = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") + 
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(seedlingSize_null, RT1a, RT1b, base = T)


## 3.2. Best model including abundance ----
# ------------------------------------

RT2a = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           abLarge +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT2b = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT1b, RT2a, RT2b, base = T)


RT3a = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           abLarge +
                           s(prevwinterT, site, bs = "re") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT3b = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           abLarge +
                           s(abLarge, site, bs = "re") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT2b, RT3a, RT3b, base = T)


RT4 = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                          abLarge +
                          s(prevwinterT, site, bs = "re") +
                          s(abLarge, site, bs = "re") +
                          s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT3a, RT4, base = T)


RT5a = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           abLarge +
                           s(prevwinterT, site, bs = "re") +
                           s(abLarge, site, bs = "re") +
                           ti(prevwinterT, abLarge, k = 3, bs = "cr") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT5b = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           abLarge +
                           s(prevwinterT, site, bs = "re") +
                           s(abLarge, site, bs = "re") +
                           prevwinterT:abLarge +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT4, RT5a, RT5b, base = T)

seedlingSize_anthropogenic = RT5a

save(seedlingSize_anthropogenic, file = "Output/Models/SeedlingSize_GAM_Anthropogenic.RData")