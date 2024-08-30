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
                      function(x) length(unique(x)))
density_per_square = aggregate(abLarge_unscaled ~ quadratID + time + site, 
                               data = droso_anthropogenic, 
                               function(x) unique(x))
yearly_density_per_square = aggregate(abLarge_unscaled ~ time + site, 
                                      data = density_per_square, 
                                      function(x) sum(x))
yearly_density_per_square$abLarge_unscaled = yearly_density_per_square$abLarge_unscaled/nbSquares$quadratID


droso_anthropogenic$time = factor(droso_anthropogenic$time)
droso_anthropogenic$site = factor(droso_anthropogenic$site)


## 2.2. Subsetting data for number of flowers only ----
# ------------------------------------------------

droso_anthropogenic$nbFlow = droso_anthropogenic$fs * droso_anthropogenic$fps

droso_nbFlow = droso_anthropogenic[- which(is.na(droso_anthropogenic$nbFlow)), ] # Reproductive individuals




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
mod_listR = list(nbFlow_null, nbFlow_prevfallR, nbFlow_prevwinterR)


# Model names
names(mod_listT) = c("null", "previous_winterT")
names(mod_listR) = c("null", "previous_fallR", "previous_winterR")

model.sel(mod_listR)

model.sel(mod_listT)

AICctab(mod_listR[["previous_fallR"]], mod_listT[["null"]]) 

# Rainfall in previous fall (September-November)
# No effect of temperature

# Checking correlation
cor(droso_nbFlow$prevfallR, droso_nbFlow$abLarge) # No
cor(droso_nbFlow$prevfallR, droso_nbFlow$size) # No

cor(droso_nbFlow$abLarge, droso_nbFlow$size) # No


# Model selection 

RT1a = gam(nbFlow ~ prevfallR + 
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT1b = gam(nbFlow ~ s(prevfallR, k = 3, bs = "cr") + 
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

AICtab(nbFlow_null, RT1a, RT1b, base = T)


## 3.2. Best model including size and abundance ----
# ---------------------------------------------

RT2a = gam(nbFlow ~ s(prevfallR, k = 3, bs = "cr") + 
                    size +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT2b = gam(nbFlow ~ s(prevfallR, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT2c = gam(nbFlow ~ s(prevfallR, k = 3, bs = "cr") + 
                    abLarge +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT2d = gam(nbFlow ~ s(prevfallR, k = 3, bs = "cr") + 
                    s(abLarge, k = 3, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT1b, RT2a, RT2b, RT2c, RT2d, base = T)


RT3a = gam(nbFlow ~ s(prevfallR, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    abLarge +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT3b = gam(nbFlow ~ s(prevfallR, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    s(abLarge, k = 3, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT2b, RT3a, RT3b, base = T)


RT4a = gam(nbFlow ~ s(prevfallR, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    s(prevfallR, site, bs = "re") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT4b = gam(nbFlow ~ s(prevfallR, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    s(size, site, bs = "re") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT2b, RT4a, RT4b, base = T)


RT5 = gam(nbFlow ~ s(prevfallR, k = 3, bs = "cr") + 
                   s(size, k = 4, bs = "cr") +
                   s(prevfallR, site, bs = "re") +
                   s(size, site, bs = "re") +
                   s(time, bs = "re") + s(site, bs = "re"), 
          data = droso_nbFlow, family = "nb", method = "REML", 
          gamma = 1.4, select = T)

AICtab(RT4a, RT5, base = T)


RT6a = gam(nbFlow ~ s(prevfallR, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    s(prevfallR, site, bs = "re") +
                    s(size, site, bs = "re") +
                    prevfallR:size +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT6b = gam(nbFlow ~ s(prevfallR, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    s(prevfallR, site, bs = "re") +
                    s(size, site, bs = "re") +
                    ti(prevfallR, size, k = 3, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT6c = gam(nbFlow ~ s(prevfallR, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    s(prevfallR, site, bs = "re") +
                    s(size, site, bs = "re") +
                    prevfallR:abLarge +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT6d = gam(nbFlow ~ s(prevfallR, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    s(prevfallR, site, bs = "re") +
                    s(size, site, bs = "re") +
                    ti(prevfallR, abLarge, k = 3, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT6e = gam(nbFlow ~ s(prevfallR, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    s(prevfallR, site, bs = "re") +
                    s(size, site, bs = "re") +
                    size:abLarge +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

RT6f = gam(nbFlow ~ s(prevfallR, k = 3, bs = "cr") + 
                    s(size, k = 4, bs = "cr") +
                    s(prevfallR, site, bs = "re") +
                    s(size, site, bs = "re") +
                    ti(size, abLarge, k = 3, bs = "cr") +
                    s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_nbFlow, family = "nb", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT5, RT6a, RT6b, RT6c, RT6d, RT6e, RT6f, base = T)

nbFlow_anthropogenic = RT5

save(nbFlow_anthropogenic, file = "Output/Models/NbFlowers_GAM_Anthropogenic.RData")