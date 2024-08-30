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


## 2.2. Subsetting data for survival only ----
# ---------------------------------------

droso_surv = droso_anthropogenic[!is.na(droso_anthropogenic$surv), ]




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

cor(droso_surv$fallR, droso_surv$abLarge) # No
cor(droso_surv$fallR, droso_surv$size) # No

cor(droso_surv$abLarge, droso_surv$size) # No



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


RT3a = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  summerT +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT3b = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  s(summerT, k = 3, bs = "cr") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT1b, RT3a, RT3b, base = T)


RT4a = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  summerT +
                  fallR:summerT +
                 s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT4b = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  summerT +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT3a, RT4a, RT4b, base = T)


## 3.2. Best model including size and abundance ----
# ---------------------------------------------

RT5a = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  summerT +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  size +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT5b = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  summerT +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(size, k = 3, bs = "cr") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT6a = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  summerT +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  abLarge +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT6b = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  summerT +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(abLarge, k = 3, bs = "cr") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT4b, RT5a, RT5b, RT6a, RT6b, base = T)


RT7a = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  summerT +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(size, k = 3, bs = "cr") +
                  abLarge +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT7b = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  summerT +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(size, k = 3, bs = "cr") +
                  s(abLarge, k = 3, bs = "cr") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT5b, RT7a, RT7b, base = T)


RT8a = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  summerT +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(size, k = 3, bs = "cr") +
                  s(fallR, site, bs = "re") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8b = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  summerT +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(size, k = 3, bs = "cr") +
                  s(summerT, site, bs = "re") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT8c = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  summerT +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(size, k = 3, bs = "cr") +
                  s(size, site, bs = "re") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT5b, RT8a, RT8b, RT8c, base = T)


RT9a = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  summerT +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(size, k = 3, bs = "cr") +
                  s(size, site, bs = "re") +
                  s(fallR, site, bs = "re") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

RT9b = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  summerT +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(size, k = 3, bs = "cr") +
                  s(size, site, bs = "re") +
                  s(summerT, site, bs = "re") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT8c, RT9a, RT9b, base = T)


RT10 = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                  summerT +
                  ti(fallR, summerT, k = 3, bs = "cr") +
                  s(size, k = 3, bs = "cr") +
                  s(size, site, bs = "re") +
                  s(fallR, site, bs = "re") +
                  s(summerT, site, bs = "re") +
                  s(time, bs = "re") + s(site, bs = "re"), 
           data = droso_surv, family = "binomial", method = "REML", 
           gamma = 1.4, select = T)

AICtab(RT9a, RT10, base = T)


RT11a = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   summerT:size +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11b = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(summerT, size, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11c = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   summerT:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11d = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(summerT, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11e = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   fallR:size +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11f = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11g = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   fallR:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11h = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11i = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   size:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT11j = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(size, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT9a, RT11a, RT11b, RT11c, RT11d, RT11e, RT11f, RT11g, RT11h, RT11i, RT11j, base = T)


RT12a = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   summerT:size +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12b = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   ti(summerT, size, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12c = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   summerT:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12d = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   ti(summerT, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12e = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   fallR:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12f = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12g = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   size:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT12h = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   ti(size, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT11f, RT12a, RT12b, RT12c, RT12d, RT12e, RT12f, RT12g, RT12h, base = T)


RT13a = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   size:abLarge +
                   summerT:size +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13b = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   size:abLarge +
                   ti(summerT, size, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13c = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   size:abLarge +
                   summerT:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13d = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   size:abLarge +
                   ti(summerT, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13e = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   size:abLarge +
                   fallR:abLarge +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

RT13f = gam(surv ~ s(fallR, k = 3, bs = "cr") + 
                   summerT +
                   ti(fallR, summerT, k = 3, bs = "cr") +
                   s(size, k = 3, bs = "cr") +
                   s(size, site, bs = "re") +
                   s(fallR, site, bs = "re") +
                   ti(fallR, size, k = 3, bs = "cr") +
                   size:abLarge +
                   ti(fallR, abLarge, k = 3, bs = "cr") +
                   s(time, bs = "re") + s(site, bs = "re"), 
            data = droso_surv, family = "binomial", method = "REML", 
            gamma = 1.4, select = T)

AICtab(RT12g, RT13a, RT13b, RT13c, RT13d, RT13e, RT13f, base = T)


surv_anthropogenic = RT12g

save(surv_anthropogenic, file = "Output/Models/Survival_GAM_Anthropogenic.RData")