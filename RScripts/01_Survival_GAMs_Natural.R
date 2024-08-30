############################################################################
#
# The aim of this script is to fit GAMs to estimate dewy-pine survival
# as a function of size, population density, temperature, and rainfall
# in various populations.
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


## 2.2. Subsetting data for survival only ----
# ---------------------------------------

droso_surv = droso_natural[!is.na(droso_natural$surv), ]



###########################################################################
#
# 3. Modelling survival ----
#
###########################################################################

# Null model 
surv_null = gam(surv ~ s(time, bs = "re") + s(site, bs = "re"),
                data = droso_surv, family = "binomial", method = "REML",
                gamma = 1.4, select = T)


## 3.1. Best model with climate only ----
# ----------------------------------

# Find the best window:

# Summer temperature (mean max temperature from May to September of the census year)
# May = Census date
# September = Census date + 4

surv_summerT = update(surv_null,    ~ . + s(summerT, bs = "cr", k = 3), 
                      data = droso_surv)

# Fall rainfall (cumulative rainfall from September to November of the census year)
# September = Census date + 4
# November = Census date + 6

surv_fallR = update(surv_null,    ~ . + s(fallR, bs = "cr", k = 3), 
                    data = droso_surv)

# Winter rainfall (cumulative rainfall from January to April of the next year)
# January = Census date + 8
# April = Census date + 11

surv_winterR = update(surv_null,    ~ . + s(winterR, bs = "cr", k = 3), 
                      data = droso_surv)

# Fall and winter rainfall (cumulative rainfall from September of the census year to April of the next year)
# September = Census date + 4
# April = Census date + 11

surv_fallwinterR = update(surv_null,    ~ . + s(fallwinterR, bs = "cr", k = 3), 
                          data = droso_surv)


# Put models in two lists (temperature and rainfall)
mod_listT = list(surv_null, surv_summerT)
mod_listR = list(surv_null, surv_fallR, surv_winterR, surv_fallwinterR)


# Model names
names(mod_listT) = c("null", "summerT")
names(mod_listR) = c("null", "fallR", "winterR", "fall_winterR")

model.sel(mod_listR)

model.sel(mod_listT)

AICctab(mod_listR[["fallR"]], mod_listT[["summerT"]]) 

# Cumulative rainfall in next fall (September to November)
# Mean max. temperature in next summer (May to September)


# Checking correlation
cor(droso_surv$summerT, droso_surv$fallR) # No
cor(droso_surv$summerT, droso_surv$abLarge) # No
cor(droso_surv$summerT, droso_surv$size) # No
cor(droso_surv$summerT, droso_surv$TSFcont) # No

cor(droso_surv$fallR, droso_surv$abLarge) # No
cor(droso_surv$fallR, droso_surv$size) # No
cor(droso_surv$fallR, droso_surv$TSFcont) # No

cor(droso_surv$abLarge, droso_surv$size) # No
cor(droso_surv$abLarge, droso_surv$TSFcont) # No

cor(droso_surv$size, droso_surv$TSFcont) # No


# Model selection with both rainfall and temperature

RT1a = gam(surv ~ fallR + 
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT1b = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT2a = gam(surv ~ summerT + 
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT2b = gam(surv ~ s(summerT, k = 3, bs = "cr") + 
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(surv_null, RT1a, RT1b, RT2a, RT2b, base = T)


RT3a = gam(surv ~ summerT + 
                  fallR +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT3b = gam(surv ~ summerT + 
                  s(fallR, k = 3, bs = "cr") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT2a, RT3a, RT3b, base = T)


RT4a = gam(surv ~ summerT + 
                  fallR +
                  fallR:summerT +
                  s(time, bs = "re") + s(site, bs = "re"), 
          data = droso_surv, family = "binomial", method = "REML", 
          gamma = 1.4, select = T)

RT4b = gam(surv ~ summerT + 
                  fallR +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT3a, RT4a, RT4b, base = T)


## 3.2. Best model including size and abundance ----
# ---------------------------------------------

RT5a = gam(surv ~ summerT + 
                  fallR +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  size +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT5b = gam(surv ~ summerT + 
                  fallR +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(size, k = 3, bs = "cr") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT6a = gam(surv ~ summerT + 
                  fallR +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  abLarge +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT6b = gam(surv ~ summerT + 
                  fallR +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(abLarge, k = 3, bs = "cr") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT7a = gam(surv ~ summerT + 
                  fallR +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  TSFcont +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT7b = gam(surv ~ summerT + 
                  fallR +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(TSFcont, k = 3, bs = "cr") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)


AICtab(RT4b, RT5a, RT5b, RT6a, RT6b, RT7a, RT7b, base = T)


RT8a = gam(surv ~ summerT + 
                  fallR +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(size, k = 3, bs = "cr") +
                  abLarge +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8b = gam(surv ~ summerT + 
                  fallR +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(size, k = 3, bs = "cr") +
                  s(abLarge, k = 3, bs = "cr") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8c = gam(surv ~ summerT + 
                  fallR +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(size, k = 3, bs = "cr") +
                  TSFcont +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8d = gam(surv ~ summerT + 
                  fallR +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(size, k = 3, bs = "cr") +
                  s(TSFcont, k = 3, bs = "cr") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT5b, RT8a, RT8b, RT8c, RT8d, base = T)


RT9a = gam(surv ~ summerT + 
                  fallR +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(size, k = 3, bs = "cr") +
                  s(TSFcont, k = 3, bs = "cr") +
                  abLarge +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT9b = gam(surv ~ summerT + 
                  fallR +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(size, k = 3, bs = "cr") +
                  s(TSFcont, k = 3, bs = "cr") +
                  s(abLarge, k = 3, bs = "cr") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT8d, RT9a, RT9b, base = T)


RT10a = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT10b = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(fallR, site, bs = "re") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT10c = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT10d = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(TSFcont, site, bs = "re") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT8d, RT10a, RT10b, RT10c, RT10d, base = T)


RT11a = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11b = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   s(size, site, bs = "re") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11c = gam(surv ~ summerT + 
              fallR +
              ti(fallR, summerT, k = 3, bs = "cr") +
              s(size, k = 3, bs = "cr") +
              s(TSFcont, k = 3, bs = "cr") +
              s(summerT, site, bs = "re") +
              s(TSFcont, site, bs = "re") +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT10a, RT11a, RT11b, RT11c, base = T)


RT12a = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   summerT:size +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12b = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(summerT, size, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12c = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   summerT:TSFcont +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12d = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(summerT, TSFcont, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12e = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   summerT:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12f = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(summerT, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12g = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   fallR:size +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12h = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12i = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   fallR:TSFcont +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12j = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(fallR, TSFcont, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12k = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   fallR:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12l = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12m = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   size:TSFcont +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12n = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12o = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   size:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12p = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12q = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   TSFcont:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12r = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(TSFcont, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT10a, RT12a, RT12b, RT12c, RT12d, RT12e, RT12f, RT12g, RT12h,
       RT12i, RT12j, RT12k, RT12l, RT12m, RT12n, RT12o, RT12p, RT12q, RT12r, base = T)


RT13a = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   summerT:size +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13b = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(summerT, size, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13c = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   summerT:TSFcont +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13d = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(summerT, TSFcont, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13e = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   summerT:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13f = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(summerT, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13g = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   fallR:size +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13h = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13i = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13j = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, TSFcont, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13k = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   fallR:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13l = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13m = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   size:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13n = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(size, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13o = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   TSFcont:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13p = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(TSFcont, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT12j, RT13a, RT13b, RT13c, RT13d, RT13e, RT13f, RT13g, RT13h,
       RT13i, RT13j, RT13k, RT13l, RT13m, RT13n, RT13o, RT13p, base = T)


RT14a = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   summerT:size +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT14b = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   ti(summerT, size, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT14c = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   summerT:TSFcont +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT14d = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   ti(summerT, TSFcont, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT14e = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   summerT:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT14f = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   ti(summerT, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT14g = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:size +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT14h = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:size +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT14i = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT14j = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   ti(fallR, TSFcont, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT14k = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   size:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT14l = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   ti(size, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT14m = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   TSFcont:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT14n = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   ti(TSFcont, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT13l, RT14a, RT14b, RT14c, RT14d, RT14e, RT14f, RT14g,
       RT14h, RT14i, RT14j, RT14k, RT14l, RT14m, RT14n, base = T)


RT15a = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   summerT:size +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT15b = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT15c = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   summerT:TSFcont +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT15d = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, TSFcont, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT15e = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   summerT:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT15f = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT15g = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   fallR:size +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT15h = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(fallR, size, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT15i = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   size:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT15j = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(size, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT15k = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   TSFcont:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT15l = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(TSFcont, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT14d, RT15a, RT15b, RT15c, RT15d, RT15e, RT15f, RT15g, RT15h, 
       RT15i, RT15j, RT15k, RT15l, base = T)


RT16a = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:TSFcont +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT16b = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   ti(summerT, TSFcont, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT16c = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT16d = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   ti(summerT, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT16e = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   fallR:size +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT16f = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT16g = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   size:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT16h = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   ti(size, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT16i = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   TSFcont:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT16j = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   ti(TSFcont, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT15d, RT16a, RT16b, RT16c, RT16d, RT16e, RT16f, RT16g, 
       RT16h, RT16i, RT16j, base = T)


RT17a = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   summerT:TSFcont +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT17b = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   ti(summerT, TSFcont, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT17c = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   fallR:size +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT17d = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   ti(fallR, size, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT17e = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   size:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT17f = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   ti(size, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT17g = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   TSFcont:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT17h = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   ti(TSFcont, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT16c, RT17a, RT17b, RT17c, RT17d, RT17e, RT17f, RT17g, RT17h, base = T)


RT18a = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   ti(summerT, TSFcont, k = 3, bs = "cr") +
                   fallR:size +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT18b = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   ti(summerT, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT18c = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   ti(summerT, TSFcont, k = 3, bs = "cr") +
                   size:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT18d = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   ti(summerT, TSFcont, k = 3, bs = "cr") +
                   ti(size, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT18e = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   ti(summerT, TSFcont, k = 3, bs = "cr") +
                   TSFcont:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT18f = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   ti(summerT, TSFcont, k = 3, bs = "cr") +
                   ti(TSFcont, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT17b, RT18a, RT18b, RT18c, RT18d, RT18e, RT18f, base = T)


RT19a = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   ti(summerT, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   size:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT19b = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   ti(summerT, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   ti(size, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT19c = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   ti(summerT, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   TSFcont:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT19d = gam(surv ~ summerT + 
                   fallR +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(TSFcont, k = 3, bs = "cr") +
                   s(summerT, site, bs = "re") +
                   ti(size, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   fallR:TSFcont +
                   ti(summerT, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   ti(summerT, TSFcont, k = 3, bs = "cr") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   ti(TSFcont, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT18b, RT19a, RT19b, RT19c, RT19d, base = T)


surv_natural = RT18b

save(surv_natural, file = "Output/Models/Survival_GAM_Natural.RData")