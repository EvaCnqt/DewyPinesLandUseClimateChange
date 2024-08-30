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


## 2.2. Subsetting data for flowering only ----
# ----------------------------------------

droso_flowering = droso_natural[- which(is.na(droso_natural$fl)), ] 




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

# Previous winter temperature (mean max temperature from January to Apr of the census year, i.e., prior to the census)
# January = Census date - 4
# Apr = Census date - 1

flowering_prevwinterT = update(flowering_null,    ~ . + s(prevwinterT, bs = "cr", k = 3), 
                               data = droso_flowering)

# Previous fall rainfall (cumulative rainfall from September to November of the year prior to the census year)
# September = Census date - 8
# November = Census date - 6

flowering_prevfallR = update(flowering_null,    ~ . + s(prevfallR, bs = "cr", k = 3), 
                             data = droso_flowering)

# Previous winter rainfall (cumulative rainfall from January to Apr of the census year, i.e., prior to the census)
# January = Census date - 4
# Apr = Census date - 1

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

AICctab(mod_listR[["previous_fallR"]], mod_listT[["previous_winterT"]]) 

# Cumulative rainfall in previous fall (September to November)
# Mean max. temperature in previous winter (January to Apr)


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

RT1a = gam(fl ~ prevwinterT + 
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT1b = gam(fl ~ s(prevwinterT, k = 3, bs = "cr") + 
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT2a = gam(fl ~ prevfallR + 
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT2b = gam(fl ~ s(prevfallR, k = 3, bs = "cr") + 
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(flowering_null, RT1a, RT1b, RT2a, RT2b, base = T)


RT3a = gam(fl ~ prevfallR + 
                prevwinterT +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT3b = gam(fl ~ prevfallR + 
                s(prevwinterT, k = 3, bs = "cr") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT2a, RT3a, RT3b, base = T)


RT4a = gam(fl ~ prevfallR + 
                prevfallR:prevwinterT +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT4b = gam(fl ~ prevfallR + 
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT2a, RT4a, RT4b, base = T)


## 3.2. Best model including size and abundance ----
# ---------------------------------------------

RT5a = gam(fl ~ prevfallR + 
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                size +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT5b = gam(fl ~ prevfallR + 
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                s(size, k = 3, bs = "cr") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT5c = gam(fl ~ prevfallR + 
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                abLarge +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT5d = gam(fl ~ prevfallR + 
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                s(abLarge, k = 3, bs = "cr") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT5e = gam(fl ~ prevfallR + 
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                TSFcont +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT5f = gam(fl ~ prevfallR + 
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                s(TSFcont, k = 3, bs = "cr") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT4b, RT5a, RT5b, RT5c, RT5d, RT5e, RT5f, base = T)


RT6a = gam(fl ~ prevfallR + 
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                size +
                abLarge +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT6b = gam(fl ~ prevfallR + 
                s(prevwinterT, k = 3, bs = "cr") +
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                size +
                s(abLarge, k = 3, bs = "cr") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT6c = gam(fl ~ prevfallR + 
                s(prevwinterT, k = 3, bs = "cr") +
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                size +
                TSFcont +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT6d = gam(fl ~ prevfallR + 
                s(prevwinterT, k = 3, bs = "cr") +
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                size +
                s(TSFcont, k = 3, bs = "cr") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT5a, RT6a, RT6b, RT6c, RT6d, base = T)


RT7a = gam(fl ~ prevfallR + 
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                size +
                s(TSFcont, k = 3, bs = "cr") +
                abLarge + 
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT7b = gam(fl ~ prevfallR + 
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                size +
                s(TSFcont, k = 3, bs = "cr") +
                s(abLarge, k = 3, bs = "cr") + 
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT6d, RT7a, RT7b, base = T)


RT8a = gam(fl ~ prevfallR + 
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                size +
                s(TSFcont, k = 3, bs = "cr") +
                s(prevfallR, site, bs = "re") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8b = gam(fl ~ prevfallR + 
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                size +
                s(TSFcont, k = 3, bs = "cr") +
                s(size, site, bs = "re") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8c = gam(fl ~ prevfallR + 
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                size +
                s(TSFcont, k = 3, bs = "cr") +
                s(TSFcont, site, bs = "re") +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT6d, RT8a, RT8b, RT8c, base = T)


RT8a = gam(fl ~ prevfallR + 
                ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                size +
                s(TSFcont, k = 3, bs = "cr") +
                prevfallR:size +
                s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8b = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(prevfallR, size, k = 3, bs = "cr") +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8c = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             prevfallR:TSFcont +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8d = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(prevfallR, TSFcont, k = 3, bs = "cr") +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8e = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             prevfallR:abLarge +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8f = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(prevfallR, abLarge, k = 3, bs = "cr") +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8g = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             size:prevwinterT +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8h = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(size, prevwinterT, k = 3, bs = "cr") +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8i = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             size:TSFcont +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8j = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(size, TSFcont, k = 3, bs = "cr") +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8k = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             size:abLarge +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8l = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(size, abLarge, k = 3, bs = "cr") +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8m = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             TSFcont:abLarge +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8n = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(TSFcont, abLarge, k = 3, bs = "cr") +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8o = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             TSFcont:prevwinterT +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8p = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT6d, RT8a, RT8b, RT8c, RT8d, RT8e, RT8f, RT8g, RT8h, RT8i, RT8j, 
       RT8k, RT8l, RT8m, RT8n, RT8o, RT8p, base = T)


RT9a = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
             prevfallR:size +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT9b = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
             ti(prevfallR, size, k = 3, bs = "cr") +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT9c = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
             prevfallR:TSFcont +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT9d = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
             ti(prevfallR, TSFcont, k = 3, bs = "cr") +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT9e = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
             prevfallR:abLarge +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT9f = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
             ti(prevfallR, abLarge, k = 3, bs = "cr") +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT9g = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
             size:TSFcont +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT9h = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
             ti(size, TSFcont, k = 3, bs = "cr") +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT9i = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
             size:abLarge +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT9j = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
             ti(size, abLarge, k = 3, bs = "cr") +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT9k = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
             size:prevwinterT +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT9l = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
             ti(size, prevwinterT, k = 3, bs = "cr") +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT9m = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
             TSFcont:abLarge +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT9n = gam(fl ~ prevfallR + 
             ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
             size +
             s(TSFcont, k = 3, bs = "cr") +
             ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
             ti(TSFcont, abLarge, k = 3, bs = "cr") +
             s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_flowering, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT8p, RT9a, RT9b, RT9c, RT9d, RT9e, RT9f, RT9g, RT9h, RT9i, RT9j, 
       RT9k, RT9l, RT9m, RT9n, base = T)


RT10a = gam(fl ~ prevfallR + 
              ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
              ti(TSFcont, abLarge, k = 3, bs = "cr") +
              prevfallR:size +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT10b = gam(fl ~ prevfallR + 
              ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
              ti(TSFcont, abLarge, k = 3, bs = "cr") +
              ti(prevfallR, size, k = 3, bs = "cr") +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT10c = gam(fl ~ prevfallR + 
              ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
              ti(TSFcont, abLarge, k = 3, bs = "cr") +
              prevfallR:TSFcont +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT10d = gam(fl ~ prevfallR + 
              ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
              ti(TSFcont, abLarge, k = 3, bs = "cr") +
              ti(prevfallR, TSFcont, k = 3, bs = "cr") +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT10e = gam(fl ~ prevfallR + 
              ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
              ti(TSFcont, abLarge, k = 3, bs = "cr") +
              prevfallR:abLarge +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT10f = gam(fl ~ prevfallR + 
              ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
              ti(TSFcont, abLarge, k = 3, bs = "cr") +
              ti(prevfallR, abLarge, k = 3, bs = "cr") +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT10g = gam(fl ~ prevfallR + 
              ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
              ti(TSFcont, abLarge, k = 3, bs = "cr") +
              size:TSFcont +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT10h = gam(fl ~ prevfallR + 
              ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
              ti(TSFcont, abLarge, k = 3, bs = "cr") +
              ti(size, TSFcont, k = 3, bs = "cr") +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT10i = gam(fl ~ prevfallR + 
              ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
              ti(TSFcont, abLarge, k = 3, bs = "cr") +
              size:abLarge +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT10j = gam(fl ~ prevfallR + 
              ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
              ti(TSFcont, abLarge, k = 3, bs = "cr") +
              ti(size, abLarge, k = 3, bs = "cr") +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT10k = gam(fl ~ prevfallR + 
              ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
              ti(TSFcont, abLarge, k = 3, bs = "cr") +
              size:prevwinterT +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT10l = gam(fl ~ prevfallR + 
              ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
              ti(TSFcont, abLarge, k = 3, bs = "cr") +
              ti(size, prevwinterT, k = 3, bs = "cr") +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT9n, RT10a, RT10b, RT10c, RT10d, RT10e, RT10f, RT10g, RT10h, RT10i, RT10j, 
       RT10k, RT10l, base = T)


RT11a = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 prevfallR:size +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11b = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(prevfallR, size, k = 3, bs = "cr") +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11c = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 prevfallR:abLarge +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11d = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(prevfallR, abLarge, k = 3, bs = "cr") +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11e = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 size:TSFcont +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11f = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(size, TSFcont, k = 3, bs = "cr") +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11g = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 size:abLarge +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11h = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(size, abLarge, k = 3, bs = "cr") +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11i = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 size:prevwinterT +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11j = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(size, prevwinterT, k = 3, bs = "cr") +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT10d, RT11a, RT11b, RT11c, RT11d, RT11e, RT11f, RT11g, RT11h, RT11i, RT11j, base = T)


RT12a = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(prevfallR, abLarge, k = 3, bs = "cr") +
                 prevfallR:size +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12b = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(prevfallR, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, size, k = 3, bs = "cr") +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12c = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(prevfallR, abLarge, k = 3, bs = "cr") +
                 size:TSFcont +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12d = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(prevfallR, abLarge, k = 3, bs = "cr") +
                 ti(size, TSFcont, k = 3, bs = "cr") +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12e = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(prevfallR, abLarge, k = 3, bs = "cr") +
                 size:abLarge +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12f = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(prevfallR, abLarge, k = 3, bs = "cr") +
                 ti(size, abLarge, k = 3, bs = "cr") +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12g = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(prevfallR, abLarge, k = 3, bs = "cr") +
                 size:prevwinterT +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12h = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(prevfallR, abLarge, k = 3, bs = "cr") +
                 ti(size, prevwinterT, k = 3, bs = "cr") +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT11d, RT12a, RT12b, RT12c, RT12d, RT12e, RT12f, RT12g, RT12h, base = T)


RT13a = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(prevfallR, abLarge, k = 3, bs = "cr") +
                 ti(size, abLarge, k = 3, bs = "cr") +
                 prevfallR:size +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13b = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(prevfallR, abLarge, k = 3, bs = "cr") +
                 ti(size, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, size, k = 3, bs = "cr") +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13c = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(prevfallR, abLarge, k = 3, bs = "cr") +
                 ti(size, abLarge, k = 3, bs = "cr") +
                 size:TSFcont +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13d = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(prevfallR, abLarge, k = 3, bs = "cr") +
                 ti(size, abLarge, k = 3, bs = "cr") +
                 ti(size, TSFcont, k = 3, bs = "cr") +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13e = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(prevfallR, abLarge, k = 3, bs = "cr") +
                 ti(size, abLarge, k = 3, bs = "cr") +
                 size:prevwinterT +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13f = gam(fl ~ prevfallR + 
                 ti(prevfallR, prevwinterT, k = 3, bs = "cr") +
                 size +
                 s(TSFcont, k = 3, bs = "cr") +
                 ti(TSFcont, prevwinterT, k = 3, bs = "cr") +
                 ti(TSFcont, abLarge, k = 3, bs = "cr") +
                 ti(prevfallR, TSFcont, k = 3, bs = "cr") +
                 ti(prevfallR, abLarge, k = 3, bs = "cr") +
                 ti(size, abLarge, k = 3, bs = "cr") +
                 ti(size, prevwinterT, k = 3, bs = "cr") +
                 s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_flowering, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT12f, RT13a, RT13b, RT13c, RT13d, RT13e, RT13f, base = T)

flowering_natural = RT12f

save(flowering_natural, file = "Output/Models/FloweringProb_GAM_Natural.RData")