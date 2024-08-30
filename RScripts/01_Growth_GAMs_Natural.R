############################################################################
#
# The aim of this script is to fit GAMs to estimate dewy-pine growth
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


## 2.2. Subsetting data for growth only ----
# -------------------------------------

droso_growth = droso_natural[!is.na(droso_natural$sizeNext), ]




###########################################################################
#
# 3. Modelling growth ----
#
###########################################################################

# Null model 
growth_null = gam(sizeNext ~ s(time, bs = "re") + s(site, bs = "re"),
                  data = droso_growth, family = "scat", method = "REML",
                  gamma = 1.4, select = T)


## 3.1. Best model with climate only ----
# ----------------------------------

# Find the best lagged window:

# Summer temperature (mean max temperature from May to September of the census year)
# May = Census date
# September = Census date + 4

growth_summerT = update(growth_null,    ~ . + s(summerT, bs = "cr", k = 3), 
                        data = droso_growth)

# Fall rainfall (cumulative rainfall from September to November of the census year)
# September = Census date + 4
# November = Census date + 6

growth_fallR = update(growth_null,    ~ . + s(fallR, bs = "cr", k = 3), 
                      data = droso_growth)

# Winter rainfall (cumulative rainfall from January to April of the next year)
# January = Census date + 8
# April = Census date + 11

growth_winterR = update(growth_null,    ~ . + s(winterR, bs = "cr", k = 3), 
                        data = droso_growth)

# Fall and winter rainfall (cumulative rainfall from September of the census year to April of the next year)
# September = Census date + 4
# April = Census date + 11

growth_fallwinterR = update(growth_null,    ~ . + s(fallwinterR, bs = "cr", k = 3), 
                            data = droso_growth)


# Put models in two lists (temperature and rainfall)
mod_listT = list(growth_null, growth_summerT)
mod_listR = list(growth_null, growth_fallR, growth_winterR, growth_fallwinterR)


# Model names
names(mod_listT) = c("null", "summerT")
names(mod_listR) = c("null", "fallR", "winterR", "fall_winterR")

model.sel(mod_listR)

model.sel(mod_listT)

AICctab(mod_listR[["fallR"]], mod_listT[["null"]]) 

# Cumulative rainfall in next fall (September to November)
# No effect of temperature


# Checking correlation
cor(droso_growth$summerT, droso_growth$fallR) # No
cor(droso_growth$summerT, droso_growth$abLarge) # No
cor(droso_growth$summerT, droso_growth$size) # No
cor(droso_growth$summerT, droso_growth$TSFcont) # No

cor(droso_growth$fallR, droso_growth$abLarge) # No
cor(droso_growth$fallR, droso_growth$size) # No
cor(droso_growth$fallR, droso_growth$TSFcont) # No

cor(droso_growth$abLarge, droso_growth$size) # No
cor(droso_growth$abLarge, droso_growth$TSFcont) # No

cor(droso_growth$size, droso_growth$TSFcont) # No


# Model selection 

RT1a = gam(sizeNext ~ fallR + 
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT1b = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") + 
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(growth_null, RT1a, RT1b, base = T)


## 3.2. Best model including size and abundance ----
# ---------------------------------------------

RT2a = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT2b = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      s(size, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT3a = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      abLarge +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT3b = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      s(abLarge, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT4a = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      TSFcont +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT4b = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      s(TSFcont, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT1b, RT2a, RT2b, RT3a, RT3b, RT4a, RT4b, base = T)


RT5a = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      abLarge +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT5b = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(abLarge, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT5c = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      TSFcont +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT5d = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT2a, RT5a, RT5b, RT5c, RT5d, base = T)


RT6a = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT6b = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      s(abLarge, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT5d, RT6a, RT6b, base = T)


RT7a = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(fallR, site, bs = "re") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT7b = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT7c = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(TSFcont, site, bs = "re") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT7d = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(abLarge, site, bs = "re") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT6a, RT7a, RT7b, RT7c, RT7d, base = T)


RT8a = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      s(fallR, site, bs = "re") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT8b = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      s(TSFcont, site, bs = "re") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT8c = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      s(abLarge, site, bs = "re") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT7b, RT8a, RT8b, RT8c, base = T)


RT9a = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      fallR:size +
                      s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT9b = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      ti(fallR, size, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT9c = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      fallR:TSFcont +
                      s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT9d = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      ti(fallR, TSFcont, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT9e = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      fallR:abLarge +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT9f = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      ti(fallR, abLarge, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)


RT9g = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      size:TSFcont +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT9h = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      ti(size, TSFcont, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT9i = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      size:abLarge +
                      s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT9j = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      ti(size, abLarge, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT9k = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      TSFcont:abLarge +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT9l = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                      size +
                      s(TSFcont, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      ti(TSFcont, abLarge, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT7b, RT9a, RT9b, RT9c, RT9d, RT9e, RT9f, RT9g, 
       RT9h, RT9i, RT9j, RT9k, RT9l, base = T)


RT10a = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                       size +
                       s(TSFcont, k = 3, bs = "cr") +
                       abLarge +
                       s(size, site, bs = "re") +
                       ti(fallR, TSFcont, k = 3, bs = "cr") +
                       fallR:size +
                       s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT10b = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                       size +
                       s(TSFcont, k = 3, bs = "cr") +
                       abLarge +
                       s(size, site, bs = "re") +
                       ti(fallR, TSFcont, k = 3, bs = "cr") +
                       ti(fallR, size, k = 3, bs = "cr") +
                       s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT10c = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                       size +
                       s(TSFcont, k = 3, bs = "cr") +
                       abLarge +
                       s(size, site, bs = "re") +
                       ti(fallR, TSFcont, k = 3, bs = "cr") +
                       fallR:abLarge +
                       s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT10d = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                       size +
                       s(TSFcont, k = 3, bs = "cr") +
                       abLarge +
                       s(size, site, bs = "re") +
                       ti(fallR, TSFcont, k = 3, bs = "cr") +
                       ti(fallR, abLarge, k = 3, bs = "cr") +
                       s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT10e = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                       size +
                       s(TSFcont, k = 3, bs = "cr") +
                       abLarge +
                       s(size, site, bs = "re") +
                       ti(fallR, TSFcont, k = 3, bs = "cr") +
                       size:TSFcont +
                       s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT10f = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                       size +
                       s(TSFcont, k = 3, bs = "cr") +
                       abLarge +
                       s(size, site, bs = "re") +
                       ti(fallR, TSFcont, k = 3, bs = "cr") +
                       ti(size, TSFcont, k = 3, bs = "cr") +
                       s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT10g = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              abLarge +
              s(size, site, bs = "re") +
              ti(fallR, TSFcont, k = 3, bs = "cr") +
              size:abLarge +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT10h = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              abLarge +
              s(size, site, bs = "re") +
              ti(fallR, TSFcont, k = 3, bs = "cr") +
              ti(size, abLarge, k = 3, bs = "cr") +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT10i = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                       size +
                       s(TSFcont, k = 3, bs = "cr") +
                       abLarge +
                       s(size, site, bs = "re") +
                       ti(fallR, TSFcont, k = 3, bs = "cr") +
                       TSFcont:abLarge +
                       s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT10j = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                       size +
                       s(TSFcont, k = 3, bs = "cr") +
                       abLarge +
                       s(size, site, bs = "re") +
                       ti(fallR, TSFcont, k = 3, bs = "cr") +
                       ti(TSFcont, abLarge, k = 3, bs = "cr") +
                       s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT9d, RT10a, RT10b, RT10c, RT10d, RT10e, RT10f, RT10g,
       RT10h, RT10i, RT10j, base = T)


RT11a = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                       size +
                       s(TSFcont, k = 3, bs = "cr") +
                       abLarge +
                       s(size, site, bs = "re") +
                       ti(fallR, TSFcont, k = 3, bs = "cr") +
                       ti(fallR, abLarge, k = 3, bs = "cr") +
                       fallR:size +
                       s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT11b = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                       size +
                       s(TSFcont, k = 3, bs = "cr") +
                       abLarge +
                       s(size, site, bs = "re") +
                       ti(fallR, TSFcont, k = 3, bs = "cr") +
                       ti(fallR, abLarge, k = 3, bs = "cr") +
                       ti(fallR, size, k = 3, bs = "cr") +
                       s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT11c = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              abLarge +
              s(size, site, bs = "re") +
              ti(fallR, TSFcont, k = 3, bs = "cr") +
              ti(fallR, abLarge, k = 3, bs = "cr") +
              size:TSFcont +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT11d = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              abLarge +
              s(size, site, bs = "re") +
              ti(fallR, TSFcont, k = 3, bs = "cr") +
              ti(fallR, abLarge, k = 3, bs = "cr") +
              ti(size, TSFcont, k = 3, bs = "cr") +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT11e = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                       size +
                       s(TSFcont, k = 3, bs = "cr") +
                       abLarge +
                       s(size, site, bs = "re") +
                       ti(fallR, TSFcont, k = 3, bs = "cr") +
                       ti(fallR, abLarge, k = 3, bs = "cr") +
                       size:abLarge +
                       s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT11f = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                       size +
                       s(TSFcont, k = 3, bs = "cr") +
                       abLarge +
                       s(size, site, bs = "re") +
                       ti(fallR, TSFcont, k = 3, bs = "cr") +
                       ti(fallR, abLarge, k = 3, bs = "cr") +
                       ti(size, abLarge, k = 3, bs = "cr") +
                       s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT11g = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              abLarge +
              s(size, site, bs = "re") +
              ti(fallR, TSFcont, k = 3, bs = "cr") +
              ti(fallR, abLarge, k = 3, bs = "cr") +
              TSFcont:abLarge +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT11h = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              abLarge +
              s(size, site, bs = "re") +
              ti(fallR, TSFcont, k = 3, bs = "cr") +
              ti(fallR, abLarge, k = 3, bs = "cr") +
              ti(TSFcont, abLarge, k = 3, bs = "cr") +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT10d, RT11a, RT11b, RT11c, RT11d, RT11e, RT11f,
       RT11g, RT11h, base = T)


RT12a = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
                       size +
                       s(TSFcont, k = 3, bs = "cr") +
                       abLarge +
                       s(size, site, bs = "re") +
                       ti(fallR, TSFcont, k = 3, bs = "cr") +
                       ti(fallR, abLarge, k = 3, bs = "cr") +
                       ti(size, abLarge, k = 3, bs = "cr") +
                       fallR:size +
                       s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT12b = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              abLarge +
              s(size, site, bs = "re") +
              ti(fallR, TSFcont, k = 3, bs = "cr") +
              ti(fallR, abLarge, k = 3, bs = "cr") +
              ti(size, abLarge, k = 3, bs = "cr") +
              ti(fallR, size, k = 3, bs = "cr") +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT12c = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              abLarge +
              s(size, site, bs = "re") +
              ti(fallR, TSFcont, k = 3, bs = "cr") +
              ti(fallR, abLarge, k = 3, bs = "cr") +
              ti(size, abLarge, k = 3, bs = "cr") +
              size:TSFcont +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT12d = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              abLarge +
              s(size, site, bs = "re") +
              ti(fallR, TSFcont, k = 3, bs = "cr") +
              ti(fallR, abLarge, k = 3, bs = "cr") +
              ti(size, abLarge, k = 3, bs = "cr") +
              ti(size, TSFcont, k = 3, bs = "cr") +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT12e = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              abLarge +
              s(size, site, bs = "re") +
              ti(fallR, TSFcont, k = 3, bs = "cr") +
              ti(fallR, abLarge, k = 3, bs = "cr") +
              ti(size, abLarge, k = 3, bs = "cr") +
              TSFcont:abLarge +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

RT12f = gam(sizeNext ~ s(fallR, k = 3, bs = "cr") +
              size +
              s(TSFcont, k = 3, bs = "cr") +
              abLarge +
              s(size, site, bs = "re") +
              ti(fallR, TSFcont, k = 3, bs = "cr") +
              ti(fallR, abLarge, k = 3, bs = "cr") +
              ti(size, abLarge, k = 3, bs = "cr") +
              ti(TSFcont, abLarge, k = 3, bs = "cr") +
              s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_growth, family = "scat", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT11f, RT12a, RT12b, RT12c, RT12d, RT12e, RT12f, base = T)

growth_natural = RT11f

save(growth_natural, file = "Output/Models/Growth_GAM_Natural.RData")