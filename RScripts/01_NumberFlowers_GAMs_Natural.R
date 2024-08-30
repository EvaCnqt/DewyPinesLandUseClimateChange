############################################################################
#
# The aim of this script is to fit GAMs to estimate dewy-pine number of
# flowers as a function of size, population density, temperature, and
# rainfall in various populations.
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
                      function(x) length(unique(x)))
density_per_square = aggregate(abLarge_unscaled ~ quadratID + time + site, 
                               data = droso_natural, 
                               function(x) unique(x))
yearly_density_per_square = aggregate(abLarge_unscaled ~ time + site, 
                                      data = density_per_square, 
                                      function(x) sum(x))
yearly_density_per_square$abLarge_unscaled = yearly_density_per_square$abLarge_unscaled/nbSquares$quadratID


droso_natural$time = factor(droso_natural$time)
droso_natural$site = factor(droso_natural$site)


## 2.2. Subsetting data for number of flowers only ----
# ------------------------------------------------

droso_natural$nbFlow = droso_natural$fs * droso_natural$fps

droso_nbFlow = droso_natural[- which(is.na(droso_natural$nbFlow)), ] # Reproductive individuals




###########################################################################
#
# 3. Modelling number of flowers ----
#
###########################################################################

# Null model 
nbFlow_null = gam(nbFlow ~ s(time, bs = "re") + s(site, bs = "re"),
                  data = droso_nbFlow, family = "nb", method = "REML",
                  gamma = 1.4, select = T)


## 3.1. Best model with climate only ----
# ----------------------------------

# Find the best lagged window:

# Previous winter temperature (mean max temperature from January to April of the census year, i.e., prior to the census)
# January = Census date - 4
# April = Census date - 1

nbFlow_prevwinterT = update(nbFlow_null,    ~ . + s(prevwinterT, bs = "cr", k = 3), 
                            data = droso_nbFlow)

# Previous fall rainfall (cumulative rainfall from September to November of the year prior to the census year)
# September = Census date - 8
# November = Census date - 6

nbFlow_prevfallR = update(nbFlow_null,    ~ . + s(prevfallR, bs = "cr", k = 3), 
                          data = droso_nbFlow)

# Previous winter rainfall (cumulative rainfall from January to April of the census year, i.e., prior to the census)
# January = Census date - 4
# April = Census date - 1

nbFlow_prevwinterR = update(nbFlow_null,    ~ . + s(prevwinterR, bs = "cr", k = 3), 
                            data = droso_nbFlow)


# Put models in two lists (temperature and rainfall)
mod_listT = list(nbFlow_null, nbFlow_prevwinterT)
mod_listR = list(nbFlow_null, nbFlow_prevfallR, nbFlow_prevwinterT)


# Model names
names(mod_listT) = c("null", "previous_winterT")
names(mod_listR) = c("null", "previous_fallR", "previous_winterR")

model.sel(mod_listR)

model.sel(mod_listT)

AICctab(mod_listR[["previous_winterR"]], mod_listT[["previous_winterT"]]) 

# Rainfall in previous fall (September-November)
# No effect of temperature


# Checking correlation
cor(droso_nbFlow$prevwinterT, droso_nbFlow$prevwinterR) # Yes
cor(droso_nbFlow$prevwinterT, droso_nbFlow$abLarge) # No
cor(droso_nbFlow$prevwinterT, droso_nbFlow$size) # No
cor(droso_nbFlow$prevwinterT, droso_nbFlow$TSFcont) # No

cor(droso_nbFlow$prevwinterR, droso_nbFlow$abLarge) # No
cor(droso_nbFlow$prevwinterR, droso_nbFlow$size) # No
cor(droso_nbFlow$prevwinterR, droso_nbFlow$TSFcont) # No

cor(droso_nbFlow$abLarge, droso_nbFlow$size) # No
cor(droso_nbFlow$abLarge, droso_nbFlow$TSFcont) # No

cor(droso_nbFlow$size, droso_nbFlow$TSFcont) # No


# Model selection 

RT1a = gam(nbFlow ~ prevwinterR + 
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT1b = gam(nbFlow ~ s(prevwinterR, k = 3, bs = "cr") + 
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT1c = gam(nbFlow ~ prevwinterT + 
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT1d = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

AICtab(nbFlow_null, RT1a, RT1b, RT1c, RT1d, base = T)


## 3.2. Best model including size and abundance ----
# ---------------------------------------------

RT2a = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    size +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT2b = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT2c = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    abLarge +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT2d = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(abLarge, k = 3, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT2e = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    TSFcont +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT2f = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(TSFcont, k = 3, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT1d, RT2a, RT2b, RT2c, RT2d, RT2e, RT2f, base = T)


RT3a = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    abLarge +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT3b = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    s(abLarge, k = 3, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT3c = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT3d = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    s(TSFcont, k = 3, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT2b, RT3a, RT3b, RT3c, RT3d, base = T)


RT4a = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    abLarge +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT4b = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    s(abLarge, k = 3, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT3c, RT4a, RT4b, base = T)


RT5a = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    s(prevwinterT, site, bs = "re") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT5b = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    s(size, site, bs = "re") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT5c = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    s(TSFcont, site, bs = "re") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT3c, RT5a, RT5b, RT5c, base = T)


RT5a = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    prevwinterT:size +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT5b = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    ti(prevwinterT, size, k = 3, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT5c = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    prevwinterT:TSFcont +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT5d = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    ti(prevwinterT, TSFcont, k = 3, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT5e = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    prevwinterT:abLarge +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT5f = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    ti(prevwinterT, abLarge, k = 3, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT5g = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    size:TSFcont +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT5h = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    ti(size, TSFcont, k = 3, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT5i = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    size:abLarge +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT5j = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    ti(size, abLarge, k = 3, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT5k = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    TSFcont:abLarge +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT5l = gam(nbFlow ~ s(prevwinterT, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    TSFcont +
                    ti(TSFcont, abLarge, k = 3, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT3c, RT5a, RT5b, RT5c, RT5d, RT5e, RT5f, RT5g, RT5h, RT5i, RT5j,
       RT5k, RT5l, base = T)

nbFlow_natural = RT3c

save(nbFlow_natural, file = "Output/Models/NbFlowers_GAM_Natural.RData")