############################################################################
#
# The aim of this script is to fit GAMs to estimate dewy-pine flowering
# probability as a function of size, population density, temperature, 
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


## 2.2. Subsetting data for flowering only ----
# ----------------------------------------

droso_flowering = droso_anthropogenic[- which(is.na(droso_anthropogenic$fl)), ] 




###########################################################################
#
# 3. Modelling flowering probability ----
#
###########################################################################

# Null model
flowering_null = gam(fl ~ s(time, bs = "re") + s(site, bs = "re"),
                     data = droso_flowering, family = "binomial", method = "REML",
                     gamma = 1.4, select = T)


## 3.1. Best model with climate only ----
# ----------------------------------

# Find the best window:

# Previous winter temperature (mean max temperature from January to April of the census year, i.e., prior to the census)
# January = Census date - 4
# April = Census date - 1

flowering_prevwinterT = update(flowering_null,    ~ . + s(prevwinterT, bs = "cr", k = 3), 
                               data = droso_flowering)

# Previous fall rainfall (cumulative rainfall from September to November of the year prior to the census year)
# September = Census date - 8
# November = Census date - 6

flowering_prevfallR = update(flowering_null,    ~ . + s(prevfallR, bs = "cr", k = 3), 
                             data = droso_flowering)

# Previous winter rainfall (cumulative rainfall from January to April of the census year, i.e., prior to the census)
# January = Census date - 4
# April = Census date - 1

flowering_prevwinterR = update(flowering_null,    ~ . + s(prevwinterR, bs = "cr", k = 3), 
                               data = droso_flowering)


# Put models in two lists (temperature and rainfall)
mod_listT = list(flowering_null, flowering_prevwinterT)
mod_listR = list(flowering_null, flowering_prevfallR, flowering_prevwinterR)


# Model names
names(mod_listT) = c("null", "previous_winterT")
names(mod_listR) = c("null", "previous_fallR", "previous_winterR")

model.sel(mod_listR)

model.sel(mod_listT)

AICctab(mod_listR[["previous_winterR"]], mod_listT[["null"]]) 

# Cumulative rainfall in previous winter (January to April)
# No effect of mean max. temperature


# Checking correlation
cor(droso_flowering$prevwinterT, droso_flowering$prevfallR) # No
cor(droso_flowering$prevwinterT, droso_flowering$abLarge) # No
cor(droso_flowering$prevwinterT, droso_flowering$size) # No
cor(droso_flowering$prevwinterT, droso_flowering$TSFcont) # No

cor(droso_flowering$prevfallR, droso_flowering$abLarge) # No
cor(droso_flowering$prevfallR, droso_flowering$size) # No
cor(droso_flowering$prevfallR, droso_flowering$TSFcont) # No

cor(droso_flowering$abLarge, droso_flowering$size) # No
cor(droso_flowering$abLarge, droso_flowering$TSFcont) # No

cor(droso_flowering$size, droso_flowering$TSFcont) # No


# Model selection

RT1a = gam(fl ~ prevwinterR + 
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT1b = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") + 
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(flowering_null, RT1a, RT1b, base = T)


## 3.2. Best model including size and abundance ----
# ---------------------------------------------

RT2a = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                size +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT2b = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                s(size, k = 3, bs = "cr") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT2c = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                abLarge +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT2d = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                s(abLarge, k = 3, bs = "cr") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT1b, RT2a, RT2b, RT2c, RT2d, base = T)


RT3a = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                s(size, k = 3, bs = "cr") +
                abLarge +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT3b = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                s(size, k = 3, bs = "cr") +
                s(abLarge, k = 3, bs = "cr") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT2b, RT3a, RT3b, base = T)


RT4a = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                s(size, k = 3, bs = "cr") +
                s(abLarge, k = 3, bs = "cr") +
                s(prevwinterR, site, bs = "re") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT4b = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                s(size, k = 3, bs = "cr") +
                s(abLarge, k = 3, bs = "cr") +
                s(size, site, bs = "re") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT4c = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                s(size, k = 3, bs = "cr") +
                s(abLarge, k = 3, bs = "cr") +
                s(abLarge, site, bs = "re") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT3b, RT4a, RT4b, RT4c, base = T)


RT5a = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                s(size, k = 3, bs = "cr") +
                s(abLarge, k = 3, bs = "cr") +
                s(size, site, bs = "re") +
                s(prevwinterR, site, bs = "re") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT5b = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                s(size, k = 3, bs = "cr") +
                s(abLarge, k = 3, bs = "cr") +
                s(size, site, bs = "re") +
                s(abLarge, site, bs = "re") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT4b, RT5a, RT5b, base = T)


RT6 = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
               s(size, k = 3, bs = "cr") +
               s(abLarge, k = 3, bs = "cr") +
               s(size, site, bs = "re") +
               s(prevwinterR, site, bs = "re") +
               s(abLarge, site, bs = "re") +
               s(time, bs = "re") + s(site, bs = "re"),
          data = droso_flowering, family = "binomial", method = "REML",
          gamma = 1.4, select = T)

AICtab(RT5a, RT6, base = T)


RT7a = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                s(size, k = 3, bs = "cr") +
                s(abLarge, k = 3, bs = "cr") +
                s(size, site, bs = "re") +
                s(prevwinterR, site, bs = "re") +
                s(abLarge, site, bs = "re") +
                prevwinterR:size +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT7b = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                s(size, k = 3, bs = "cr") +
                s(abLarge, k = 3, bs = "cr") +
                s(size, site, bs = "re") +
                s(prevwinterR, site, bs = "re") +
                s(abLarge, site, bs = "re") +
                ti(prevwinterR, size, k = 3, bs = "cr") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT7c = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                s(size, k = 3, bs = "cr") +
                s(abLarge, k = 3, bs = "cr") +
                s(size, site, bs = "re") +
                s(prevwinterR, site, bs = "re") +
                s(abLarge, site, bs = "re") +
                prevwinterR:abLarge +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT7d = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                s(size, k = 3, bs = "cr") +
                s(abLarge, k = 3, bs = "cr") +
                s(size, site, bs = "re") +
                s(prevwinterR, site, bs = "re") +
                s(abLarge, site, bs = "re") +
                ti(prevwinterR, abLarge, k = 3, bs = "cr") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT7e = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                s(size, k = 3, bs = "cr") +
                s(abLarge, k = 3, bs = "cr") +
                s(size, site, bs = "re") +
                s(prevwinterR, site, bs = "re") +
                s(abLarge, site, bs = "re") +
                size:abLarge +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT7f = gam(fl ~ s(prevwinterR, k = 3, bs = "cr") +
                s(size, k = 3, bs = "cr") +
                s(abLarge, k = 3, bs = "cr") +
                s(size, site, bs = "re") +
                s(prevwinterR, site, bs = "re") +
                s(abLarge, site, bs = "re") +
                ti(size, abLarge, k = 3, bs = "cr") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT6, RT7a, RT7b, RT7c, RT7d, RT7e, RT7f, base = T)


flowering_anthropogenic = RT6

save(flowering_anthropogenic, file = "Output/Models/FloweringProb_GAM_Anthropogenic.RData")