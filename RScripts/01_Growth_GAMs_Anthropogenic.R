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


## 2.2. Subsetting data for growth only ----
# -------------------------------------

droso_growth = droso_anthropogenic[!is.na(droso_anthropogenic$sizeNext), ]




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

AICctab(mod_listR[["null"]], mod_listT[["summerT"]]) 

# No effect of cumulative rainfall
# Mean max. temperature in next summer (May to September)


# Checking correlation
cor(droso_growth$summerT, droso_growth$abLarge) # No
cor(droso_growth$summerT, droso_growth$size) # No

cor(droso_growth$fallR, droso_growth$abLarge) # No
cor(droso_growth$fallR, droso_growth$size) # No

cor(droso_growth$abLarge, droso_growth$size) # No


# Model selection 

RT1a = gam(sizeNext ~ summerT + 
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT1b = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") + 
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(growth_null, RT1a, RT1b, base = T)


## 3.2. Best model including size and abundance ----
# ---------------------------------------------

RT2a = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      size +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT2b = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      s(size, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT3a = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      abLarge +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT3b = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      s(abLarge, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT1b, RT2a, RT2b, RT3a, RT3b, base = T)


RT4a = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      s(size, k = 3, bs = "cr") +
                      abLarge +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT4b = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      s(size, k = 3, bs = "cr") +
                      s(abLarge, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT2b, RT4a, RT4b, base = T)


RT5a = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      s(size, k = 3, bs = "cr") +
                      abLarge +
                      s(summerT, site, bs = "re") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT5b = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      s(size, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT5c = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      s(size, k = 3, bs = "cr") +
                      abLarge +
                      s(abLarge, site, bs = "re") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT4a, RT5a, RT5b, RT5c, base = T)


RT6a = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      s(size, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      s(summerT, site, bs = "re") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT6b = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      s(size, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      s(abLarge, site, bs = "re") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT5b, RT6a, RT6b, base = T)


RT7 = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                     s(size, k = 3, bs = "cr") +
                     abLarge +
                     s(size, site, bs = "re") +
                     s(summerT, site, bs = "re") +
                     s(abLarge, site, bs = "re") +
                     s(time, bs = "re") + s(site, bs = "re"), 
          data = droso_growth, family = "scat", method = "REML", 
          gamma = 1.4, select = T)

AICtab(RT6a, RT7, base = T)


RT8a = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      s(size, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      s(summerT, site, bs = "re") +
                      summerT:size +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT8b = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      s(size, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      s(summerT, site, bs = "re") +
                      ti(summerT, size, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT8c = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      s(size, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      s(summerT, site, bs = "re") +
                      summerT:abLarge +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT8d = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      s(size, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      s(summerT, site, bs = "re") +
                      ti(summerT, abLarge, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT8e = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      s(size, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      s(summerT, site, bs = "re") +
                      size:abLarge +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

RT8f = gam(sizeNext ~ s(summerT, k = 3, bs = "cr") +
                      s(size, k = 3, bs = "cr") +
                      abLarge +
                      s(size, site, bs = "re") +
                      s(summerT, site, bs = "re") +
                      ti(size, abLarge, k = 3, bs = "cr") +
                      s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_growth, family = "scat", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT6a, RT8a, RT8b, RT8c, RT8d, RT8e, RT8f, base = T)

growth_anthropogenic = RT6a

save(growth_anthropogenic, file = "Output/Models/Growth_GAM_Anthropogenic.RData")