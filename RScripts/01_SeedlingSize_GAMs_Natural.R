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

droso_natural = read.csv("Data/droso_natural.csv")




###########################################################################
#
# 2. Preparing data ----
#
###########################################################################

## 2.1. Calculating mean density per quadrat for density standardization ----
# ----------------------------------------------------------------------

nbSquares = aggregate(quadratID ~ time + site, 
                      data = droso_natural, 
                      FUN = function(x) length(unique(x)))
density_per_square = aggregate(abLarge_unscaled ~ quadratID + time + site,
                               data = droso_natural, 
                               function(x) unique(x))
yearly_density_per_square = aggregate(abLarge_unscaled ~ time + site,
                                      data = density_per_square, 
                                      function(x) sum(x))
yearly_density_per_square$abLarge_unscaled = yearly_density_per_square$abLarge_unscaled/nbSquares$quadratID


droso_natural$time = factor(droso_natural$time)
droso_natural$site = factor(droso_natural$site)


## 2.2. Subsetting data for seedling size only ----
# --------------------------------------------

droso_seedlingSize = droso_natural[which(droso_natural$stage == "SD"), ] # Seedlings




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

RT2c = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           TSFcont +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT2d = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(TSFcont, k = 3, bs = "cr") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT1b, RT2a, RT2b, RT2c, RT2d, base = T)


RT3a = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT3b = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           s(TSFcont, k = 3, bs = "cr") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT2b, RT3a, RT3b, base = T)


RT4a = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(prevwinterT, site, bs = "re") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT4b = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(abLarge, site, bs = "re") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT4c = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(TSFcont, site, bs = "re") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT3a, RT4a, RT4b, RT4c, base = T)


RT5a = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(abLarge, site, bs = "re") +
                           s(prevwinterT, site, bs = "re") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT5b = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(abLarge, site, bs = "re") +
                           s(TSFcont, site, bs = "re") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT4b, RT5a, RT5b, base = T)


RT6a = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(abLarge, site, bs = "re") +
                           s(TSFcont, site, bs = "re") +
                           prevwinterT:abLarge +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT6b = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(abLarge, site, bs = "re") +
                           s(TSFcont, site, bs = "re") +
                           ti(prevwinterT, abLarge, k = 3, bs = "cr") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT6c = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(abLarge, site, bs = "re") +
                           s(TSFcont, site, bs = "re") +
                           prevwinterT:TSFcont +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT6d = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(abLarge, site, bs = "re") +
                           s(TSFcont, site, bs = "re") +
                           ti(prevwinterT, TSFcont, k = 3, bs = "cr") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT6e = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(abLarge, site, bs = "re") +
                           s(TSFcont, site, bs = "re") +
                           abLarge:TSFcont +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT6f = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(abLarge, site, bs = "re") +
                           s(TSFcont, site, bs = "re") +
                           ti(abLarge, TSFcont, k = 3, bs = "cr") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT5b, RT6a, RT6b, RT6c, RT6d, RT6e, RT6f, base = T)


RT7a = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(abLarge, site, bs = "re") +
                           s(TSFcont, site, bs = "re") +
                           ti(prevwinterT, abLarge, k = 3, bs = "cr") +
                           prevwinterT:TSFcont +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT7b = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(abLarge, site, bs = "re") +
                           s(TSFcont, site, bs = "re") +
                           ti(prevwinterT, abLarge, k = 3, bs = "cr") +
                           ti(prevwinterT, TSFcont, k = 3, bs = "cr") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT7c = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(abLarge, site, bs = "re") +
                           s(TSFcont, site, bs = "re") +
                           ti(prevwinterT, abLarge, k = 3, bs = "cr") +
                           abLarge:TSFcont +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT7d = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(abLarge, site, bs = "re") +
                           s(TSFcont, site, bs = "re") +
                           ti(prevwinterT, abLarge, k = 3, bs = "cr") +
                           ti(abLarge, TSFcont, k = 3, bs = "cr") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT6b, RT7a, RT7b, RT7c, RT7d, base = T)


RT8a = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(abLarge, site, bs = "re") +
                           s(TSFcont, site, bs = "re") +
                           ti(prevwinterT, abLarge, k = 3, bs = "cr") +
                           ti(prevwinterT, TSFcont, k = 3, bs = "cr") +
                           abLarge:TSFcont +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT8b = gam(size_unscaled ~ s(prevwinterT, k = 3, bs = "cr") +
                           s(abLarge, k = 3, bs = "cr") +
                           TSFcont +
                           s(abLarge, site, bs = "re") +
                           s(TSFcont, site, bs = "re") +
                           ti(prevwinterT, abLarge, k = 3, bs = "cr") +
                           ti(prevwinterT, TSFcont, k = 3, bs = "cr") +
                           ti(abLarge, TSFcont, k = 3, bs = "cr") +
                           s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_seedlingSize, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT7b, RT8a, RT8b, base = T)

seedlingSize_natural = RT8b

save(seedlingSize_natural, file = "Output/Models/SeedlingSize_GAM_Natural.RData")