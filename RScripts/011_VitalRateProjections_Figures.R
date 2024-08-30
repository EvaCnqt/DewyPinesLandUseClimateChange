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

library(ggplot2)
library(viridis)
library(ggbeeswarm)
library(plyr)
library(cowplot)
library(patchwork)


## 1.3. Loading data ----
# ------------------

# Dewy-pine data
droso_natural = read.csv("Data/droso_natural.csv")
droso_anthropogenic = read.csv("Data/droso_anthropogenic.csv")


# Vital-rate models
load("Output/Models/Survival_GAM_Anthropogenic.RData")
load("Output/Models/Survival_GAM_Natural.RData")

load("Output/Models/Growth_GAM_Anthropogenic.RData")
load("Output/Models/Growth_GAM_Natural.RData")

load("Output/Models/FloweringProb_GAM_Anthropogenic.RData")
load("Output/Models/FloweringProb_GAM_Natural.RData")

load("Output/Models/NbFlowers_GAM_Anthropogenic.RData")
load("Output/Models/NbFlowers_GAM_Natural.RData")

load("Output/Models/SeedlingSize_GAM_Anthropogenic.RData")
load("Output/Models/SeedlingSize_GAM_Natural.RData")




###########################################################################
#
# 2. Preparing data ----
#
###########################################################################

## 2.1. Calculating mean density per quadrat for density standardization ----
# ----------------------------------------------------------------------

nbSquares_natural = aggregate(quadratID ~ time + site, 
                              data = droso_natural, 
                              function(x) length(unique(x)))
density_per_square_natural = aggregate(abLarge_unscaled ~ quadratID + time + site, 
                                       data = droso_natural, 
                                       function(x) unique(x))
yearly_density_per_square_natural = aggregate(abLarge_unscaled ~ time + site, 
                                              data = density_per_square_natural, 
                                              function(x) sum(x))
yearly_density_per_square_natural$abLarge_unscaled = yearly_density_per_square_natural$abLarge_unscaled/nbSquares_natural$quadratID

nbSquares_anthropogenic = aggregate(quadratID ~ time + site, 
                                data = droso_anthropogenic, 
                                function(x) length(unique(x)))
density_per_square_anthropogenic = aggregate(abLarge_unscaled ~ quadratID + time + site, 
                                         data = droso_anthropogenic, 
                                         function(x) unique(x))
yearly_density_per_square_anthropogenic = aggregate(abLarge_unscaled ~ time + site, 
                                                data = density_per_square_anthropogenic, 
                                                function(x) sum(x))
yearly_density_per_square_anthropogenic$abLarge_unscaled = yearly_density_per_square_anthropogenic$abLarge_unscaled/nbSquares_anthropogenic$quadratID


## 2.2. Formatting time and site as factors ----
# -----------------------------------------

droso_natural$time = factor(droso_natural$time)
droso_natural$site = factor(droso_natural$site)

droso_anthropogenic$time = factor(droso_anthropogenic$time)
droso_anthropogenic$site = factor(droso_anthropogenic$site)


## 2.3. Covariate values ----
# ----------------------

fallR_values_anthropogenic = round(seq(quantile(droso_anthropogenic$fallR_unscaled, probs = 0.20), 
                                   quantile(droso_anthropogenic$fallR_unscaled, probs = 0.80), 
                                   length.out = 50), 2)
prevfallR_values_anthropogenic = round(seq(quantile(droso_anthropogenic$prevfallR_unscaled, probs = 0.20), 
                                       quantile(droso_anthropogenic$prevfallR_unscaled, probs = 0.80), 
                                       length.out = 50), 2)
prevwinterR_values_anthropogenic = round(seq(quantile(droso_anthropogenic$prevwinterR_unscaled, probs = 0.20), 
                                         quantile(droso_anthropogenic$prevwinterR_unscaled, probs = 0.80), 
                                         length.out = 50), 2)
summerT_values_anthropogenic = round(seq(quantile(droso_anthropogenic$summerT_unscaled, probs = 0.20), 
                                     quantile(droso_anthropogenic$summerT_unscaled, probs = 0.80), 
                                     length.out = 50), 2)
prevwinterT_values_anthropogenic = round(seq(quantile(droso_anthropogenic$prevwinterT_unscaled, probs = 0.20), 
                                         quantile(droso_anthropogenic$prevwinterT_unscaled, probs = 0.80), 
                                         length.out = 50), 2)
size_values_anthropogenic = round(seq(min(droso_anthropogenic$size_unscaled), 
                                  max(droso_anthropogenic$size_unscaled), 
                                  length.out = 50), 2)
abLarge_values_anthropogenic = unique(round(seq(min(droso_anthropogenic$abLarge_unscaled), 
                                            max(droso_anthropogenic$abLarge_unscaled), 
                                            length.out = 50)))

fallR_timeSeries_anthropogenic = aggregate(fallR_unscaled ~ time + site, 
                                       droso_anthropogenic, mean)
prevfallR_timeSeries_anthropogenic = aggregate(prevfallR_unscaled ~ time + site, 
                                           droso_anthropogenic, mean)
prevwinterR_timeSeries_anthropogenic = aggregate(prevwinterR_unscaled ~ time + site, 
                                             droso_anthropogenic, mean)
summerT_timeSeries_anthropogenic = aggregate(summerT_unscaled ~ time + site, 
                                         droso_anthropogenic, mean)
prevwinterT_timeSeries_anthropogenic = aggregate(prevwinterT_unscaled ~ time + site, 
                                             droso_anthropogenic, mean)


fallR_values_natural = round(seq(quantile(droso_natural$fallR_unscaled, probs = 0.20), 
                                 quantile(droso_natural$fallR_unscaled, probs = 0.80), 
                                 length.out = 50), 2)
prevfallR_values_natural = round(seq(quantile(droso_natural$prevfallR_unscaled, probs = 0.20), 
                                     quantile(droso_natural$prevfallR_unscaled, probs = 0.80), 
                                     length.out = 50), 2)
prevwinterR_values_natural = round(seq(quantile(droso_natural$prevwinterR_unscaled, probs = 0.20), 
                                       quantile(droso_natural$prevwinterR_unscaled, probs = 0.80), 
                                       length.out = 50), 2)
summerT_values_natural = round(seq(quantile(droso_natural$summerT_unscaled, probs = 0.20),
                                   quantile(droso_natural$summerT_unscaled, probs = 0.80), 
                                   length.out = 50), 2)
prevwinterT_values_natural = round(seq(quantile(droso_natural$prevwinterT_unscaled, probs = 0.20),
                                       quantile(droso_natural$prevwinterT_unscaled, probs = 0.80), 
                                       length.out = 50), 2)
size_values_natural = round(seq(min(droso_natural$size_unscaled), 
                                max(droso_natural$size_unscaled), 
                                length.out = 50), 2)
abLarge_values_natural = unique(round(seq(min(droso_natural$abLarge_unscaled), 
                                          max(droso_natural$abLarge_unscaled), 
                                          length.out = 50)))
TSFcont_values_natural = unique(round(seq(min(droso_natural$TSFcont_unscaled), 
                                          max(droso_natural$TSFcont_unscaled), 
                                          length.out = 50)))

fallR_timeSeries_natural = aggregate(fallR_unscaled ~ time + site, 
                                     droso_natural, mean)
prevfallR_timeSeries_natural = aggregate(prevfallR_unscaled ~ time + site, 
                                         droso_natural, mean)
prevwinterR_timeSeries_natural = aggregate(prevwinterR_unscaled ~ time + site, 
                                           droso_natural, mean)
summerT_timeSeries_natural = aggregate(summerT_unscaled ~ time + site, 
                                       droso_natural, mean)
prevwinterT_timeSeries_natural = aggregate(prevwinterT_unscaled ~ time + site, 
                                           droso_natural, mean)




###########################################################################
#
# 3. Plotting predictions - Main text ----
#
###########################################################################

## Plot theme function

colBG = "white" # Plot background color
colPlot = "black"     # Plot color
font = "Helvetica"    # Plot font
fontSize = 10         # Plot font size


# Define ggplot theme
theme_general = function(){ 
  
  theme_minimal() %+replace%    # Replace elements we want to change
    
    theme(axis.text.x = element_text(colour = colPlot, size = 7, family = font, 
                                     margin = margin(t = 2, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(colour = colPlot, size = 7, family = font, 
                                     margin = margin(t = 0, r = 2, b = 0, l = 0)),
          axis.title.x = element_text(colour = colPlot, size = 8, family = font, 
                                      margin = margin(t = 5, r = 0, b = 0, l = 0)), 
          axis.title.y = element_text(colour = colPlot, size = 8, family = font, 
                                      margin = margin(t = 0, r = 5, b = 0, l = 0), angle = 90),
          strip.text = element_text(colour = colPlot, size = 8, family = font, 
                                    margin = margin(t = 0, r = 0, b = 5, l = 0)),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = colPlot),
          axis.ticks = element_line(colour = colPlot),
          plot.background = element_rect(fill = colBG, color = colBG, linewidth = 0),
          plot.title = element_blank(),
          legend.position = "right", 
          legend.justification = "right", 
          legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.text = element_text(colour = colPlot, size = 7, family = font),
          legend.title = element_text(colour = colPlot, size = 8, family = font),
          legend.key.width = unit(35, "pt")) 
}


## 3.1. Mean vital rates per habitat ----
# ----------------------------------

mean_vr = rbind(expand.grid(size = 0,
                            fallR = 0,
                            summerT = 0,
                            prevfallR = 0,
                            prevwinterR = 0,
                            prevwinterT = 0,
                            abLarge = 0,
                            TSFcont = 0,
                            time = unique(droso_anthropogenic$time[which(droso_anthropogenic$time != 2022)]),
                            site = unique(droso_anthropogenic$site),
                            habitat = "Anthropogenic",
                            vr = c("Survival", "Growth", "Flowering",
                                   "Number of flowers", "Seedling size")),
                expand.grid(size = 0,
                            fallR = 0,
                            summerT = 0,
                            prevfallR = 0,
                            prevwinterR = 0,
                            prevwinterT = 0,
                            abLarge = 0,
                            TSFcont = 0,
                            time = unique(droso_natural$time[which(droso_natural$time != 2022)]),
                            site = unique(droso_natural$site),
                            habitat = "Natural",
                            vr = c("Survival", "Growth", "Flowering",
                                   "Number of flowers", "Seedling size")))

mean_vr$vr_mean_pred = NA
mean_vr$vr_mean_obs = NA


# Survival

# Natural
mean_vr$vr_mean_pred[which(mean_vr$habitat == "Natural" & 
                             mean_vr$vr == "Survival")] = plogis(predict(surv_natural, 
                                                                            newdata = mean_vr[which(mean_vr$habitat == "Natural" & 
                                                                                                      mean_vr$vr == "Survival"), ], 
                                                                            se.fit = T)$fit)

mean_vr$vr_mean_obs[which(mean_vr$habitat == "Natural" & 
                            mean_vr$vr == "Survival")] = apply(mean_vr[which(mean_vr$habitat == "Natural" & 
                                                                               mean_vr$vr == "Survival"), ], 1,
                                                               function(x) mean(droso_natural$surv[which(!is.na(droso_natural$surv) &
                                                                                                           droso_natural$time == as.character(x[9]) &
                                                                                                           droso_natural$site == as.character(x[10]))], na.rm = T))

# Anthropogenic
mean_vr$vr_mean_pred[which(mean_vr$habitat == "Anthropogenic" & 
                             mean_vr$vr == "Survival")] = plogis(predict(surv_anthropogenic, 
                                                                            newdata = mean_vr[which(mean_vr$habitat == "Anthropogenic" & 
                                                                                                      mean_vr$vr == "Survival"), ], 
                                                                            se.fit = T)$fit)

mean_vr$vr_mean_obs[which(mean_vr$habitat == "Anthropogenic" & 
                            mean_vr$vr == "Survival")] = apply(mean_vr[which(mean_vr$habitat == "Anthropogenic" & 
                                                                               mean_vr$vr == "Survival"), ], 1,
                                                               function(x) mean(droso_anthropogenic$surv[which(!is.na(droso_anthropogenic$surv) &
                                                                                                             droso_anthropogenic$time == as.character(x[9]) &
                                                                                                             droso_anthropogenic$site == as.character(x[10]))], na.rm = T))


# Growth

# Natural
mean_vr$vr_mean_pred[which(mean_vr$habitat == "Natural" & 
                             mean_vr$vr == "Growth")] = predict(growth_natural, 
                                                                newdata = mean_vr[which(mean_vr$habitat == "Natural" &
                                                                                          mean_vr$vr == "Growth"), ], 
                                                                se.fit = T)$fit

mean_vr$vr_mean_obs[which(mean_vr$habitat == "Natural" & 
                            mean_vr$vr == "Growth")] = apply(mean_vr[which(mean_vr$habitat == "Natural" & 
                                                                             mean_vr$vr == "Growth"), ], 1,
                                                             function(x) mean(droso_natural$sizeNext[which(!is.na(droso_natural$sizeNext) &
                                                                                                             droso_natural$time == as.character(x[9]) &
                                                                                                             droso_natural$site == as.character(x[10]))], na.rm = T))

# Anthropogenic
mean_vr$vr_mean_pred[which(mean_vr$habitat == "Anthropogenic" & 
                             mean_vr$vr == "Growth")] = predict(growth_anthropogenic, 
                                                                newdata = mean_vr[which(mean_vr$habitat == "Anthropogenic" &
                                                                                          mean_vr$vr == "Growth"), ], 
                                                                se.fit = T)$fit

mean_vr$vr_mean_obs[which(mean_vr$habitat == "Anthropogenic" & 
                            mean_vr$vr == "Growth")] = apply(mean_vr[which(mean_vr$habitat == "Anthropogenic" & 
                                                                             mean_vr$vr == "Growth"), ], 1,
                                                             function(x) mean(droso_anthropogenic$sizeNext[which(!is.na(droso_anthropogenic$sizeNext) &
                                                                                                               droso_anthropogenic$time == as.character(x[9]) &
                                                                                                               droso_anthropogenic$site == as.character(x[10]))], na.rm = T))


# Seedling size

# Natural
mean_vr$vr_mean_pred[which(mean_vr$habitat == "Natural" & 
                             mean_vr$vr == "Seedling size")] = predict(seedlingSize_natural, 
                                                                       newdata = mean_vr[which(mean_vr$habitat == "Natural" &
                                                                                                 mean_vr$vr == "Seedling size"), ], 
                                                                       se.fit = T)$fit

mean_vr$vr_mean_obs[which(mean_vr$habitat == "Natural" & 
                            mean_vr$vr == "Seedling size")] = apply(mean_vr[which(mean_vr$habitat == "Natural" & 
                                                                                    mean_vr$vr == "Seedling size"), ], 1,
                                                                    function(x) mean(droso_natural$size_unscaled[which(!is.na(droso_natural$size_unscaled) &
                                                                                                                         droso_natural$stage == "SD" &
                                                                                                                         droso_natural$time == as.character(x[9]) &
                                                                                                                         droso_natural$site == as.character(x[10]))], na.rm = T))

# Anthropogenic
mean_vr$vr_mean_pred[which(mean_vr$habitat == "Anthropogenic" & 
                             mean_vr$vr == "Seedling size")] = predict(seedlingSize_anthropogenic, 
                                                                       newdata = mean_vr[which(mean_vr$habitat == "Anthropogenic" &
                                                                                                 mean_vr$vr == "Seedling size"), ],
                                                                       se.fit = T)$fit

mean_vr$vr_mean_obs[which(mean_vr$habitat == "Anthropogenic" & 
                            mean_vr$vr == "Seedling size")] = apply(mean_vr[which(mean_vr$habitat == "Anthropogenic" & 
                                                                                    mean_vr$vr == "Seedling size"), ], 1,
                                                                    function(x) mean(droso_anthropogenic$size_unscaled[which(!is.na(droso_anthropogenic$size_unscaled) &
                                                                                                                           droso_anthropogenic$stage == "SD" &
                                                                                                                           droso_anthropogenic$time == as.character(x[9]) &
                                                                                                                           droso_anthropogenic$site == as.character(x[10]))], na.rm = T))




# Flowering

# Natural
mean_vr$vr_mean_pred[which(mean_vr$habitat == "Natural" & 
                             mean_vr$vr == "Flowering")] = plogis(predict(flowering_natural, 
                                                                             newdata = mean_vr[which(mean_vr$habitat == "Natural" & 
                                                                                                       mean_vr$vr == "Flowering"), ],
                                                                             se.fit = T)$fit)

mean_vr$vr_mean_obs[which(mean_vr$habitat == "Natural" & 
                            mean_vr$vr == "Flowering")] = apply(mean_vr[which(mean_vr$habitat == "Natural" &
                                                                                mean_vr$vr == "Flowering"), ], 1,
                                                                function(x) mean(droso_natural$fl[which(!is.na(droso_natural$fl) &
                                                                                                          droso_natural$time == as.character(x[9]) &
                                                                                                          droso_natural$site == as.character(x[10]))], na.rm = T))

# Anthropogenic
mean_vr$vr_mean_pred[which(mean_vr$habitat == "Anthropogenic" & 
                             mean_vr$vr == "Flowering")] = plogis(predict(flowering_anthropogenic, 
                                                                             newdata = mean_vr[which(mean_vr$habitat == "Anthropogenic" &
                                                                                                       mean_vr$vr == "Flowering"), ],
                                                                             se.fit = T)$fit)

mean_vr$vr_mean_obs[which(mean_vr$habitat == "Anthropogenic" & 
                            mean_vr$vr == "Flowering")] = apply(mean_vr[which(mean_vr$habitat == "Anthropogenic" &
                                                                                mean_vr$vr == "Flowering"), ], 1,
                                                                function(x) mean(droso_anthropogenic$fl[which(!is.na(droso_anthropogenic$fl) &
                                                                                                            droso_anthropogenic$time == as.character(x[9]) &
                                                                                                            droso_anthropogenic$site == as.character(x[10]))], na.rm = T))


# Number of flowers

droso_natural$nbFlow = droso_natural$fs * droso_natural$fps
droso_anthropogenic$nbFlow = droso_anthropogenic$fs * droso_anthropogenic$fps

mean_vr$size[which(mean_vr$habitat == "Anthropogenic")] = (mean(droso_anthropogenic$size_unscaled[which(!is.na(droso_anthropogenic$nbFlow))]) - mean(droso_anthropogenic$size_unscaled)) / (2 * sd(droso_anthropogenic$size_unscaled))
mean_vr$size[which(mean_vr$habitat == "Natural")] = (mean(droso_natural$size_unscaled[which(!is.na(droso_natural$nbFlow))]) - mean(droso_natural$size_unscaled)) / (2 * sd(droso_natural$size_unscaled))


# Natural
mean_vr$vr_mean_pred[which(mean_vr$habitat == "Natural" & 
                             mean_vr$vr == "Number of flowers")] = exp(predict(nbFlow_natural, 
                                                                               newdata = mean_vr[which(mean_vr$habitat == "Natural" & 
                                                                                                         mean_vr$vr == "Number of flowers"), ],
                                                                               se.fit = T)$fit)

mean_vr$vr_mean_obs[which(mean_vr$habitat == "Natural" & 
                            mean_vr$vr == "Number of flowers")] = apply(mean_vr[which(mean_vr$habitat == "Natural" &
                                                                                        mean_vr$vr == "Number of flowers"), ], 1,
                                                                        function(x) mean(droso_natural$nbFlow[which(droso_natural$fl == 1 &
                                                                                                                      !is.na(droso_natural$nbFlow) &
                                                                                                                      droso_natural$time == as.character(x[9]) &
                                                                                                                      droso_natural$site == as.character(x[10]))], na.rm = T))

# Anthropogenic
mean_vr$vr_mean_pred[which(mean_vr$habitat == "Anthropogenic" & 
                             mean_vr$vr == "Number of flowers")] = exp(predict(nbFlow_anthropogenic, 
                                                                               newdata = mean_vr[which(mean_vr$habitat == "Anthropogenic" & 
                                                                                                         mean_vr$vr == "Number of flowers"), ],
                                                                               se.fit = T)$fit)

mean_vr$vr_mean_obs[which(mean_vr$habitat == "Anthropogenic" & 
                            mean_vr$vr == "Number of flowers")] = apply(mean_vr[which(mean_vr$habitat == "Anthropogenic" &
                                                                                        mean_vr$vr == "Number of flowers"), ], 1,
                                                                        function(x) mean(droso_anthropogenic$nbFlow[which(droso_anthropogenic$fl == 1 &
                                                                                                                        !is.na(droso_anthropogenic$nbFlow) &
                                                                                                                        droso_anthropogenic$time == as.character(x[9]) &
                                                                                                                        droso_anthropogenic$site == as.character(x[10]))], na.rm = T))


mean_vr$vr_mean_obs[which(is.nan(mean_vr$vr_mean_obs))] = NA


mean_vr$habitat = factor(mean_vr$habitat, 
                         levels = c("Natural", "Anthropogenic"))

mean_vr_CI = aggregate(vr_mean_pred ~ habitat + vr, 
                       data = mean_vr,
                       FUN = function(x) quantile(x, probs = c(0.025, 0.975), na.rm = T))

png(filename = "Output/Plots/Figure2.png", 
    width = 8,
    height = 12,
    units = "cm",
    bg = "white",
    res = 600)

ggplot(data = mean_vr, mapping = aes(x = habitat, y = vr_mean_pred)) +
  facet_wrap(~ vr, scales = "free", nrow = 3) +
  geom_boxplot(ymin = mean_vr_CI$vr_mean[, 1],
               ymax = mean_vr_CI$vr_mean[, 2],
               outlier.size = 0.1, size = 0.1,
               linewidth = 0.4,
               fill = alpha("grey", 0.5)) + 
  geom_point(aes(y = mean_vr$vr_mean_obs),
             position = position_quasirandom(width = 0.1),
             size = 0.8, alpha = 0.5, shape = 16,
             colour = "lightsalmon") +
  stat_summary(fun = mean, aes(shape = "Mean"), geom = "point", shape = 17, size = 1, position = position_dodge(width = 0.75)) +
  xlab("Habitat type") +
  ylab("Vital-rate value") +
  theme_general()

dev.off()


## 3.2. Survival as a function of temperature per rainfall (other covariates at average value) ----
# --------------------------------------------------------------------------------------------

surv_pred_noRE = rbind(expand.grid(size = 0,
                                   fallR = (c(100, 150, 200) - mean(fallR_timeSeries_anthropogenic$fallR_unscaled)) / (2 * sd(fallR_timeSeries_anthropogenic$fallR_unscaled)),
                                   summerT = (summerT_values_anthropogenic - mean(summerT_timeSeries_anthropogenic$summerT_unscaled)) / (2 * sd(summerT_timeSeries_anthropogenic$summerT_unscaled)),
                                   abLarge = 0,
                                   TSFcont = 0,
                                   dormancy = "Anthropogenic"),
                       expand.grid(size = 0,
                                   fallR = (c(100, 150, 200) - mean(fallR_timeSeries_natural$fallR_unscaled)) / (2 * sd(fallR_timeSeries_natural$fallR_unscaled)),
                                   summerT = (summerT_values_natural - mean(summerT_timeSeries_natural$summerT_unscaled)) / (2 * sd(summerT_timeSeries_natural$summerT_unscaled)),
                                   abLarge = 0,
                                   TSFcont = 0,
                                   dormancy = "Natural"))


surv_pred_noRE$surv = NA
surv_pred_noRE$upr = NA
surv_pred_noRE$lwr = NA


surv_pred_noRE$surv[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                     newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                     se.fit = T,
                                                                                     exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                     newdata.guaranteed = T)$fit)
surv_pred_noRE$upr[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                    newdata.guaranteed = T)$fit + 
                                                                              1.96 * predict(surv_natural, 
                                                                                             newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                             newdata.guaranteed = T)$se.fit)
surv_pred_noRE$lwr[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                    newdata.guaranteed = T)$fit - 
                                                                              1.96 * predict(surv_natural, 
                                                                                             newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                             newdata.guaranteed = T)$se.fit)

surv_pred_noRE$surv[which(surv_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                       newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                       se.fit = T,
                                                                                       exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                       newdata.guaranteed = T)$fit)
surv_pred_noRE$upr[which(surv_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                      newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                      se.fit = T,
                                                                                      exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                      newdata.guaranteed = T)$fit + 
                                                                                1.96 * predict(surv_anthropogenic, 
                                                                                               newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                               se.fit = T,
                                                                                               exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                               newdata.guaranteed = T)$se.fit)
surv_pred_noRE$lwr[which(surv_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                      newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                      se.fit = T,
                                                                                      exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                      newdata.guaranteed = T)$fit - 
                                                                                1.96 * predict(surv_anthropogenic, 
                                                                                               newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                               se.fit = T,
                                                                                               exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                               newdata.guaranteed = T)$se.fit)


surv_pred_noRE$summerT_unscaled[which(surv_pred_noRE$dormancy == "Natural")] = surv_pred_noRE$summerT[which(surv_pred_noRE$dormancy == "Natural")] * (2 * sd(summerT_timeSeries_natural$summerT_unscaled)) + mean(summerT_timeSeries_natural$summerT_unscaled)
surv_pred_noRE$summerT_unscaled[which(surv_pred_noRE$dormancy == "Anthropogenic")] = surv_pred_noRE$summerT[which(surv_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(summerT_timeSeries_anthropogenic$summerT_unscaled)) + mean(summerT_timeSeries_anthropogenic$summerT_unscaled)

surv_pred_noRE$fallR_unscaled[which(surv_pred_noRE$dormancy == "Natural")] = surv_pred_noRE$fallR[which(surv_pred_noRE$dormancy == "Natural")] * (2 * sd(fallR_timeSeries_natural$fallR_unscaled)) + mean(fallR_timeSeries_natural$fallR_unscaled)
surv_pred_noRE$fallR_unscaled[which(surv_pred_noRE$dormancy == "Anthropogenic")] = surv_pred_noRE$fallR[which(surv_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(fallR_timeSeries_anthropogenic$fallR_unscaled)) + mean(fallR_timeSeries_anthropogenic$fallR_unscaled)

surv_pred_noRE$fallR_unscaled = as.factor(surv_pred_noRE$fallR_unscaled)

droso_natural$summerT_unscaled_rounded = round(droso_natural$summerT_unscaled, 2)
droso_natural$fallR_unscaled_rounded = round_any(droso_natural$fallR_unscaled, 10)

surv_pred_noRE$dormancy = factor(surv_pred_noRE$dormancy, 
                                 levels = c("Natural", "Anthropogenic"))

surv_plot_temp_rain = ggplot(surv_pred_noRE, aes(summerT_unscaled, surv,
                                                 col = fallR_unscaled, group = fallR_unscaled)) +
  facet_wrap(~ dormancy) + 
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = fallR_unscaled), alpha = 0.1, col = NA) +
  scale_colour_manual(name = "Next fall\ncumulative rainfall\n(Sep-Nov) (in mm)",
                      values = palette.colors(3, palette = "Dark 2")) +
  scale_fill_manual(name = "Next fall\ncumulative rainfall\n(Sep-Nov) (in mm)",
                    values = palette.colors(3, palette = "Dark 2")) +
  xlab("Next summer average max. daily\ntemperature (May-Sep) (in °C)") +
  ylab("Survival\nprobability") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"))




## 3.3. Survival as a function of rainfall per density (other covariates at average value) ----
# --------------------------------------------------------------------------------------------

surv_pred_noRE = expand.grid(size = 0,
                             fallR = (fallR_values_natural - mean(fallR_timeSeries_natural$fallR_unscaled)) / (2 * sd(fallR_timeSeries_natural$fallR_unscaled)),
                             summerT = 0,
                             abLarge = (c(2, 6, 10) - mean(yearly_density_per_square_natural$abLarge_unscaled)) / (2 * sd(yearly_density_per_square_natural$abLarge_unscaled)),
                             TSFcont = 0,dormancy = "Natural")


surv_pred_noRE$surv = NA
surv_pred_noRE$upr = NA
surv_pred_noRE$lwr = NA


surv_pred_noRE$surv[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                     newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                     se.fit = T,
                                                                                     exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                     newdata.guaranteed = T)$fit)
surv_pred_noRE$upr[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                    newdata.guaranteed = T)$fit + 
                                                                              1.96 * predict(surv_natural, 
                                                                                             newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                             newdata.guaranteed = T)$se.fit)
surv_pred_noRE$lwr[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                    newdata.guaranteed = T)$fit - 
                                                                              1.96 * predict(surv_natural, 
                                                                                             newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                             newdata.guaranteed = T)$se.fit)

surv_pred_noRE$abLarge_unscaled[which(surv_pred_noRE$dormancy == "Natural")] = surv_pred_noRE$abLarge[which(surv_pred_noRE$dormancy == "Natural")] * (2 * sd(yearly_density_per_square_natural$abLarge_unscaled)) + mean(yearly_density_per_square_natural$abLarge_unscaled)

surv_pred_noRE$fallR_unscaled[which(surv_pred_noRE$dormancy == "Natural")] = surv_pred_noRE$fallR[which(surv_pred_noRE$dormancy == "Natural")] * (2 * sd(fallR_timeSeries_natural$fallR_unscaled)) + mean(fallR_timeSeries_natural$fallR_unscaled)

surv_pred_noRE$abLarge_unscaled = as.factor(surv_pred_noRE$abLarge_unscaled)

surv_plot_rain_density = ggplot(surv_pred_noRE, aes(fallR_unscaled, surv,
                                                 col = abLarge_unscaled, group = abLarge_unscaled)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = abLarge_unscaled), alpha = 0.1, col = NA) +
  scale_colour_manual(name = "Aboveground\ndensity of large\nindividuals\n(in ind./m2)",
                      values = palette.colors(3, palette = "Dark 2")) +
  scale_fill_manual(name = "Aboveground\ndensity of large\nindividuals\n(in ind./m2)",
                    values = palette.colors(3, palette = "Dark 2")) +
  xlab("Next fall\ncumulative rainfall\n(Sep-Nov) (in mm)") +
  ylab("Survival\nprobability") +
  ggtitle("Natural") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"),
        axis.line.y = element_line(colour = "white"),
        axis.ticks.y = element_line(colour = "white"),
        axis.title.y = element_text(colour = "white"),
        axis.text.y = element_text(colour = "white"),
        plot.title = element_text(colour = colPlot, size = 8, family = font,
                                  margin = margin(t = 0, r = 0, b = 5, l = 0)))




## 3.4. Survival as a function of temperature per density (other covariates at average value) ----
# --------------------------------------------------------------------------------------------

surv_pred_noRE = expand.grid(size = 0,
                             fallR = 0,
                             summerT = (summerT_values_natural - mean(summerT_timeSeries_natural$summerT_unscaled)) / (2 * sd(summerT_timeSeries_natural$summerT_unscaled)),
                             abLarge = (c(2, 6, 10) - mean(yearly_density_per_square_natural$abLarge_unscaled)) / (2 * sd(yearly_density_per_square_natural$abLarge_unscaled)),
                             TSFcont = 0,dormancy = "Natural")


surv_pred_noRE$surv = NA
surv_pred_noRE$upr = NA
surv_pred_noRE$lwr = NA


surv_pred_noRE$surv[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                     newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                     se.fit = T,
                                                                                     exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                     newdata.guaranteed = T)$fit)
surv_pred_noRE$upr[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                    newdata.guaranteed = T)$fit + 
                                                                              1.96 * predict(surv_natural, 
                                                                                             newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                             newdata.guaranteed = T)$se.fit)
surv_pred_noRE$lwr[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                    newdata.guaranteed = T)$fit - 
                                                                              1.96 * predict(surv_natural, 
                                                                                             newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                             newdata.guaranteed = T)$se.fit)

surv_pred_noRE$abLarge_unscaled[which(surv_pred_noRE$dormancy == "Natural")] = surv_pred_noRE$abLarge[which(surv_pred_noRE$dormancy == "Natural")] * (2 * sd(yearly_density_per_square_natural$abLarge_unscaled)) + mean(yearly_density_per_square_natural$abLarge_unscaled)

surv_pred_noRE$summerT_unscaled[which(surv_pred_noRE$dormancy == "Natural")] = surv_pred_noRE$summerT[which(surv_pred_noRE$dormancy == "Natural")] * (2 * sd(summerT_timeSeries_natural$summerT_unscaled)) + mean(summerT_timeSeries_natural$summerT_unscaled)

surv_pred_noRE$abLarge_unscaled = as.factor(surv_pred_noRE$abLarge_unscaled)

surv_plot_temp_density = ggplot(surv_pred_noRE, aes(summerT_unscaled, surv,
                                                    col = abLarge_unscaled, group = abLarge_unscaled)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = abLarge_unscaled), alpha = 0.1, col = NA) +
  scale_colour_manual(name = "Aboveground\ndensity of large\nindividuals\n(in ind./m2)",
                      values = palette.colors(3, palette = "Dark 2")) +
  scale_fill_manual(name = "Aboveground\ndensity of large\nindividuals\n(in ind./m2)",
                    values = palette.colors(3, palette = "Dark 2")) +
  xlab("Next summer average\nmax. daily temperature\n(May-Sep) (in ºC)") +
  ylab("Survival\nprobability") +
  ggtitle("Natural") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(colour = colPlot, size = 8, family = font,
                                  margin = margin(t = 0, r = 0, b = 5, l = 0)))


## 3.5. Growth as a function of rainfall (other covariates at average value) ----
# ----------------------------------------------------------------------------

growth_pred_noRE = expand.grid(size = 0,
                               fallR = (fallR_values_natural - mean(fallR_timeSeries_natural$fallR_unscaled)) / (2 * sd(fallR_timeSeries_natural$fallR_unscaled)),
                               summerT = 0,
                               abLarge = 0,
                               TSFcont = 0,
                               dormancy = "Natural")


growth_pred_noRE$sizeNext = NA
growth_pred_noRE$upr = NA
growth_pred_noRE$lwr = NA


growth_pred_noRE$sizeNext[which(growth_pred_noRE$dormancy == "Natural")] = (predict(growth_natural, 
                                                                                    newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Natural"), ], 
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)", "s(site)", "s(size,site)"), 
                                                                                    newdata.guaranteed = T)$fit)
growth_pred_noRE$upr[which(growth_pred_noRE$dormancy == "Natural")] = (predict(growth_natural, 
                                                                               newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Natural"), ], 
                                                                               se.fit = T,
                                                                               exclude = c("s(time)", "s(site)", "s(size,site)"), 
                                                                               newdata.guaranteed = T)$fit + 
                                                                         1.96 * predict(growth_natural, 
                                                                                        newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Natural"), ], 
                                                                                        se.fit = T,
                                                                                        exclude = c("s(time)", "s(site)", "s(size,site)"), 
                                                                                        newdata.guaranteed = T)$se.fit)
growth_pred_noRE$lwr[which(growth_pred_noRE$dormancy == "Natural")] = (predict(growth_natural, 
                                                                               newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Natural"), ], 
                                                                               se.fit = T,
                                                                               exclude = c("s(time)", "s(site)", "s(size,site)"), 
                                                                               newdata.guaranteed = T)$fit - 
                                                                         1.96 * predict(growth_natural, 
                                                                                        newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Natural"), ], 
                                                                                        se.fit = T,
                                                                                        exclude = c("s(time)", "s(site)", "s(size,site)"), 
                                                                                        newdata.guaranteed = T)$se.fit)

growth_pred_noRE$fallR_unscaled[which(growth_pred_noRE$dormancy == "Natural")] = growth_pred_noRE$fallR[which(growth_pred_noRE$dormancy == "Natural")] * (2 * sd(fallR_timeSeries_natural$fallR_unscaled)) + mean(fallR_timeSeries_natural$fallR_unscaled)


growth_plot_rainfall = ggplot(growth_pred_noRE, aes(fallR_unscaled, sizeNext)) +
  facet_wrap(~ dormancy, scales = "free_x") +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, col = NA) +
  xlab("Next fall cumulative\nrainfall (Sep-Nov) (in mm)") +
  ylab("\nSize next year") +
  theme_general()


## 3.6. Growth as a function of density (other covariates at average value) ----
# ------------------------------------------------------------------------

growth_pred_noRE = rbind(expand.grid(size = 0,
                                     fallR = 0,
                                     summerT = 0,
                                     abLarge = (abLarge_values_anthropogenic - mean(yearly_density_per_square_anthropogenic$abLarge_unscaled)) / (2 * sd(yearly_density_per_square_anthropogenic$abLarge_unscaled)),
                                     TSFcont = 0,
                                     dormancy = "Anthropogenic"),
                         expand.grid(size = 0,
                                     fallR = 0,
                                     summerT = 0,
                                     abLarge = (abLarge_values_natural - mean(yearly_density_per_square_natural$abLarge_unscaled)) / (2 * sd(yearly_density_per_square_natural$abLarge_unscaled)),
                                     TSFcont = 0,
                                     dormancy = "Natural"))


growth_pred_noRE$sizeNext = NA
growth_pred_noRE$upr = NA
growth_pred_noRE$lwr = NA


growth_pred_noRE$sizeNext[which(growth_pred_noRE$dormancy == "Natural")] = (predict(growth_natural, 
                                                                                    newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Natural"), ], 
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)", "s(site)", "s(size,site)"), 
                                                                                    newdata.guaranteed = T)$fit)
growth_pred_noRE$upr[which(growth_pred_noRE$dormancy == "Natural")] = (predict(growth_natural, 
                                                                               newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Natural"), ], 
                                                                               se.fit = T,
                                                                               exclude = c("s(time)", "s(site)", "s(size,site)"), 
                                                                               newdata.guaranteed = T)$fit + 
                                                                         1.96 * predict(growth_natural, 
                                                                                        newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Natural"), ], 
                                                                                        se.fit = T,
                                                                                        exclude = c("s(time)", "s(site)", "s(size,site)"), 
                                                                                        newdata.guaranteed = T)$se.fit)
growth_pred_noRE$lwr[which(growth_pred_noRE$dormancy == "Natural")] = (predict(growth_natural, 
                                                                               newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Natural"), ], 
                                                                               se.fit = T,
                                                                               exclude = c("s(time)", "s(site)", "s(size,site)"), 
                                                                               newdata.guaranteed = T)$fit - 
                                                                         1.96 * predict(growth_natural, 
                                                                                        newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Natural"), ], 
                                                                                        se.fit = T,
                                                                                        exclude = c("s(time)", "s(site)", "s(size,site)"), 
                                                                                        newdata.guaranteed = T)$se.fit)

growth_pred_noRE$sizeNext[which(growth_pred_noRE$dormancy == "Anthropogenic")] = (predict(growth_anthropogenic, 
                                                                                      newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                      se.fit = T,
                                                                                      exclude = c("s(time)", "s(site)", "s(size,site)", "s(summerT,site)"), 
                                                                                      newdata.guaranteed = T)$fit)
growth_pred_noRE$upr[which(growth_pred_noRE$dormancy == "Anthropogenic")] = (predict(growth_anthropogenic, 
                                                                                 newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                 se.fit = T,
                                                                                 exclude = c("s(time)", "s(site)", "s(size,site)", "s(summerT,site)"), 
                                                                                 newdata.guaranteed = T)$fit + 
                                                                           1.96 * predict(growth_anthropogenic, 
                                                                                          newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                          se.fit = T,
                                                                                          exclude = c("s(time)", "s(site)", "s(size,site)", "s(summerT,site)"), 
                                                                                          newdata.guaranteed = T)$se.fit)
growth_pred_noRE$lwr[which(growth_pred_noRE$dormancy == "Anthropogenic")] = (predict(growth_anthropogenic, 
                                                                                 newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                 se.fit = T,
                                                                                 exclude = c("s(time)", "s(site)", "s(size,site)", "s(summerT,site)"), 
                                                                                 newdata.guaranteed = T)$fit - 
                                                                           1.96 * predict(growth_anthropogenic, 
                                                                                          newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                          se.fit = T,
                                                                                          exclude = c("s(time)", "s(site)", "s(size,site)", "s(summerT,site)"), 
                                                                                          newdata.guaranteed = T)$se.fit)

growth_pred_noRE$abLarge_unscaled[which(growth_pred_noRE$dormancy == "Natural")] = growth_pred_noRE$abLarge[which(growth_pred_noRE$dormancy == "Natural")] * (2 * sd(yearly_density_per_square_natural$abLarge_unscaled)) + mean(yearly_density_per_square_natural$abLarge_unscaled)
growth_pred_noRE$abLarge_unscaled[which(growth_pred_noRE$dormancy == "Anthropogenic")] = growth_pred_noRE$abLarge[which(growth_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(yearly_density_per_square_anthropogenic$abLarge_unscaled)) + mean(yearly_density_per_square_anthropogenic$abLarge_unscaled)

growth_pred_noRE$dormancy = factor(growth_pred_noRE$dormancy, 
                                   levels = c("Natural", "Anthropogenic"))


growth_plot_density = ggplot(growth_pred_noRE, aes(abLarge_unscaled, sizeNext)) +
  facet_wrap(~ dormancy, scales = "free_x") +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, col = NA) +
  xlab("Aboveground density of\nlarge individuals (in ind./m2)") +
  ylab("\nSize next year") +
  theme_general()




## 3.7. Flowering as a function of density (other covariates at average value) ----
# --------------------------------------------------------------------------------------------

flowering_pred_noRE = rbind(expand.grid(size = 0,
                                        prevfallR = 0,
                                        prevwinterR = 0,
                                        prevwinterT = 0,
                                        abLarge = (abLarge_values_anthropogenic - mean(yearly_density_per_square_anthropogenic$abLarge_unscaled)) / (2 * sd(yearly_density_per_square_anthropogenic$abLarge_unscaled)),
                                        TSFcont = 0,
                                        dormancy = "Anthropogenic"),
                            expand.grid(size = 0,
                                        prevfallR = 0,
                                        prevwinterR = 0,
                                        prevwinterT = 0,
                                        abLarge = (abLarge_values_natural - mean(yearly_density_per_square_natural$abLarge_unscaled)) / (2 * sd(yearly_density_per_square_natural$abLarge_unscaled)),
                                        TSFcont = 0,
                                        dormancy = "Natural"))


flowering_pred_noRE$flowering = NA
flowering_pred_noRE$upr = NA
flowering_pred_noRE$lwr = NA


flowering_pred_noRE$flowering[which(flowering_pred_noRE$dormancy == "Natural")] = plogis(predict(flowering_natural, 
                                                                                                    newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                                    se.fit = T,
                                                                                                    exclude = c("s(time)", "s(site)"), 
                                                                                                    newdata.guaranteed = T)$fit)
flowering_pred_noRE$upr[which(flowering_pred_noRE$dormancy == "Natural")] = plogis(predict(flowering_natural, 
                                                                                              newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                              se.fit = T,
                                                                                              exclude = c("s(time)", "s(site)"), 
                                                                                              newdata.guaranteed = T)$fit + 
                                                                                        1.96 * predict(flowering_natural, 
                                                                                                       newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                                       se.fit = T,
                                                                                                       exclude = c("s(time)", "s(site)"), 
                                                                                                       newdata.guaranteed = T)$se.fit)
flowering_pred_noRE$lwr[which(flowering_pred_noRE$dormancy == "Natural")] = plogis(predict(flowering_natural, 
                                                                                              newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                              se.fit = T,
                                                                                              exclude = c("s(time)", "s(site)"), 
                                                                                              newdata.guaranteed = T)$fit - 
                                                                                        1.96 * predict(flowering_natural, 
                                                                                                       newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                                       se.fit = T,
                                                                                                       exclude = c("s(time)", "s(site)"), 
                                                                                                       newdata.guaranteed = T)$se.fit)

flowering_pred_noRE$flowering[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(flowering_anthropogenic, 
                                                                                                      newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                      se.fit = T,
                                                                                                      exclude = c("s(time)", "s(site)", "s(size,site)", "s(prevwinterR,site)", "s(abLarge,site)"), 
                                                                                                      newdata.guaranteed = T)$fit)
flowering_pred_noRE$upr[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(flowering_anthropogenic, 
                                                                                                newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                se.fit = T,
                                                                                                exclude = c("s(time)", "s(site)", "s(size,site)", "s(prevwinterR,site)", "s(abLarge,site)"), 
                                                                                                newdata.guaranteed = T)$fit + 
                                                                                          1.96 * predict(flowering_anthropogenic, 
                                                                                                         newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                         se.fit = T,
                                                                                                         exclude = c("s(time)", "s(site)", "s(size,site)", "s(prevwinterR,site)", "s(abLarge,site)"), 
                                                                                                         newdata.guaranteed = T)$se.fit)
flowering_pred_noRE$lwr[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(flowering_anthropogenic, 
                                                                                                newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                se.fit = T,
                                                                                                exclude = c("s(time)", "s(site)", "s(size,site)", "s(prevwinterR,site)", "s(abLarge,site)"), 
                                                                                                newdata.guaranteed = T)$fit - 
                                                                                          1.96 * predict(flowering_anthropogenic, 
                                                                                                         newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                         se.fit = T,
                                                                                                         exclude = c("s(time)", "s(site)", "s(size,site)", "s(prevwinterR,site)", "s(abLarge,site)"), 
                                                                                                         newdata.guaranteed = T)$se.fit)


flowering_pred_noRE$abLarge_unscaled[which(flowering_pred_noRE$dormancy == "Natural")] = flowering_pred_noRE$abLarge[which(flowering_pred_noRE$dormancy == "Natural")] * (2 * sd(yearly_density_per_square_natural$abLarge_unscaled)) + mean(yearly_density_per_square_natural$abLarge_unscaled)
flowering_pred_noRE$abLarge_unscaled[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = flowering_pred_noRE$abLarge[which(flowering_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(yearly_density_per_square_anthropogenic$abLarge_unscaled)) + mean(yearly_density_per_square_anthropogenic$abLarge_unscaled)

flowering_pred_noRE$dormancy = factor(flowering_pred_noRE$dormancy, 
                                      levels = c("Natural", "Anthropogenic"))

flowering_plot_density = ggplot(flowering_pred_noRE, aes(abLarge_unscaled, flowering)) +
  facet_wrap(~ dormancy, scales = "free_x") +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, col = NA) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  xlab("Aboveground density of\nlarge individuals (in ind./m2)") +
  ylab("Flowering\nprobability") +
  theme_general()


## 3.8. Flowering as a function of temperature per rainfall (other covariates at average value) ----
# --------------------------------------------------------------------------------------------

flowering_pred_noRE = expand.grid(size = 0,
                                  prevfallR = (c(100, 150, 180) - mean(prevfallR_timeSeries_natural$prevfallR_unscaled)) / (2 * sd(prevfallR_timeSeries_natural$prevfallR_unscaled)),
                                  prevwinterT = (prevwinterT_values_natural - mean(prevwinterT_timeSeries_natural$prevwinterT_unscaled)) / (2 * sd(prevwinterT_timeSeries_natural$prevwinterT_unscaled)),
                                  abLarge = 0,
                                  TSFcont = 0,
                                  dormancy = "Natural")


flowering_pred_noRE$flowering = NA
flowering_pred_noRE$upr = NA
flowering_pred_noRE$lwr = NA


flowering_pred_noRE$flowering[which(flowering_pred_noRE$dormancy == "Natural")] = plogis(predict(flowering_natural,
                                                                                     newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ],
                                                                                     se.fit = T,
                                                                                     exclude = c("s(time)", "s(site)"),
                                                                                     newdata.guaranteed = T)$fit)
flowering_pred_noRE$upr[which(flowering_pred_noRE$dormancy == "Natural")] = plogis(predict(flowering_natural,
                                                                                    newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ],
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)", "s(site)"),
                                                                                    newdata.guaranteed = T)$fit +
                                                                              1.96 * predict(flowering_natural,
                                                                                             newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ],
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)", "s(site)"),
                                                                                             newdata.guaranteed = T)$se.fit)
flowering_pred_noRE$lwr[which(flowering_pred_noRE$dormancy == "Natural")] = plogis(predict(flowering_natural,
                                                                                    newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ],
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)", "s(site)"),
                                                                                    newdata.guaranteed = T)$fit -
                                                                              1.96 * predict(flowering_natural,
                                                                                             newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ],
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)", "s(site)"),
                                                                                             newdata.guaranteed = T)$se.fit)

flowering_pred_noRE$prevwinterT_unscaled[which(flowering_pred_noRE$dormancy == "Natural")] = flowering_pred_noRE$prevwinterT[which(flowering_pred_noRE$dormancy == "Natural")] * (2 * sd(prevwinterT_timeSeries_natural$prevwinterT_unscaled)) + mean(prevwinterT_timeSeries_natural$prevwinterT_unscaled)

flowering_pred_noRE$prevfallR_unscaled[which(flowering_pred_noRE$dormancy == "Natural")] = flowering_pred_noRE$prevfallR[which(flowering_pred_noRE$dormancy == "Natural")] * (2 * sd(prevfallR_timeSeries_natural$prevfallR_unscaled)) + mean(prevfallR_timeSeries_natural$prevfallR_unscaled)

flowering_pred_noRE$prevfallR_unscaled = as.factor(flowering_pred_noRE$prevfallR_unscaled)


flowering_plot_temp_rain = ggplot(flowering_pred_noRE, aes(prevwinterT_unscaled, flowering,
                                                 col = prevfallR_unscaled, 
                                                 group = prevfallR_unscaled)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = prevfallR_unscaled), alpha = 0.1, col = NA) +
  scale_colour_manual(name = "Previous fall\ncumulative rainfall\n(Sep-Nov) (in mm)",
                      values = palette.colors(3, palette = "Dark 2")) +
  scale_fill_manual(name = "Previous fall\ncumulative rainfall\n(Sep-Nov) (in mm)",
                    values = palette.colors(3, palette = "Dark 2")) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  xlab("Previous winter average max. daily\ntemperature (Jan-Apr) (in °C)") +
  ylab("Flowering\nprobability") +
  ggtitle("Natural") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(colour = colPlot, size = 8, family = font,
                                  margin = margin(t = 0, r = 0, b = 5, l = 0)))



## 3.9. Plot final figure ----
# -----------------------

surv_plot_temp_rain_legend = get_plot_component(surv_plot_temp_rain, 'guide-box-right')
surv_plot_rain_density_legend = get_plot_component(surv_plot_rain_density, 'guide-box-right')
flowering_plot_temp_rain_legend = get_plot_component(flowering_plot_temp_rain, 'guide-box-right')

png(filename = "Output/Plots/Figure3.png", 
    width = 16, 
    height = 15, 
    units = "cm", 
    bg = "white", 
    res = 600, 
    type = "cairo")

ggdraw() +
  draw_plot(surv_plot_temp_rain + theme(legend.position = "none"), x = 0, y = .68, width = .4, height = .33) +
  draw_plot(growth_plot_rainfall, x = 0.56, y = .68, width = .28, height = .33) +
  draw_plot(surv_plot_temp_rain_legend, x = 0.455, y = .71, width = .1, height = .33) +
  draw_plot(growth_plot_density, x = 0, y = 0.35, width = .4, height = .33) +
  draw_plot(flowering_plot_density, x = .46, y = 0.35, width = .4, height = .33) +
  draw_plot(flowering_plot_temp_rain + theme(legend.position = "none"), x = 0, y = 0.015, width = .28, height = .33) +
  draw_plot(flowering_plot_temp_rain_legend, x = 0.338, y = 0.045, width = .1, height = .33) +
  draw_plot(surv_plot_rain_density + theme(legend.position = "none"), x = 0.625, y = 0, width = .245, height = .348) +
  draw_plot(surv_plot_temp_density + theme(legend.position = "none"), x = 0.455, y = 0, width = .245, height = .348) +
  draw_plot(surv_plot_rain_density_legend, x = 0.9, y = .06, width = .1, height = .33) +
  
  draw_plot_label(label = c("(a)", "(b)", "(c)",
                            "(d)", "(e)", "(f)"),
                  size = 10,
                  fontface = "plain",
                  x = c(0.005, 0.565, 0.005, 0.465, 0.005, 0.465),
                  y = c(1, 1, 0.67, 0.67, 0.335, 0.335))
  
dev.off() 




###########################################################################
#
# 4. Plotting predictions - Appendix - Among-site differences in mean ----
#
###########################################################################

mean_size_repro_anthropogenic = mean(droso_anthropogenic$size_unscaled[which(!is.na(droso_anthropogenic$nbFlow))])
mean_size_repro_natural = mean(droso_natural$size_unscaled[which(!is.na(droso_natural$nbFlow))])


## 4.1. Average site-specific vital rates ----
# ---------------------------------------

# Survival
surv_pred = expand.grid(size = 0,
                        fallR = 0,
                        summerT = 0,
                        abLarge = 0,
                        TSFcont = 0,
                        time = levels(droplevels(droso_anthropogenic$time[which(droso_anthropogenic$time != 2022)])),
                        site = c(levels(droso_anthropogenic$site), levels(droso_natural$site)))

surv_pred$vr = "Survival probability"
surv_pred$dormancy = "Anthropogenic"
surv_pred$dormancy[which(surv_pred$site %in% levels(droso_natural$site))] = "Natural"

surv_pred$mean_pred[which(surv_pred$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                newdata = surv_pred[which(surv_pred$dormancy == "Natural"), ], 
                                                                                se.fit = T)$fit)
surv_pred$upr[which(surv_pred$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                               newdata = surv_pred[which(surv_pred$dormancy == "Natural"), ], 
                                                                               se.fit = T)$fit + 
                                                                         1.96 * predict(surv_natural, 
                                                                                        newdata = surv_pred[which(surv_pred$dormancy == "Natural"), ], 
                                                                                        se.fit = T)$se.fit)
surv_pred$lwr[which(surv_pred$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                               newdata = surv_pred[which(surv_pred$dormancy == "Natural"), ], 
                                                                               se.fit = T)$fit -
                                                                         1.96 * predict(surv_natural, 
                                                                                        newdata = surv_pred[which(surv_pred$dormancy == "Natural"), ], 
                                                                                        se.fit = T)$se.fit)

surv_pred$mean_obs[which(surv_pred$dormancy == "Natural" & 
                         surv_pred$vr == "Survival probability")] = apply(surv_pred[which(surv_pred$dormancy == "Natural" & 
                                                                                surv_pred$vr == "Survival probability"), ], 1,
                                                               function(x) mean(droso_natural$surv[which(!is.na(droso_natural$surv) &
                                                                                                           droso_natural$time == as.character(x[6]) &
                                                                                                           droso_natural$site == as.character(x[7]))], na.rm = T))


surv_pred$mean_pred[which(surv_pred$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                  newdata = surv_pred[which(surv_pred$dormancy == "Anthropogenic"), ], 
                                                                                  se.fit = T)$fit)
surv_pred$upr[which(surv_pred$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                 newdata = surv_pred[which(surv_pred$dormancy == "Anthropogenic"), ], 
                                                                                 se.fit = T)$fit + 
                                                                           1.96 * predict(surv_anthropogenic, 
                                                                                          newdata = surv_pred[which(surv_pred$dormancy == "Anthropogenic"), ], 
                                                                                          se.fit = T)$se.fit)
surv_pred$lwr[which(surv_pred$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                 newdata = surv_pred[which(surv_pred$dormancy == "Anthropogenic"), ], 
                                                                                 se.fit = T)$fit -
                                                                           1.96 * predict(surv_anthropogenic, 
                                                                                          newdata = surv_pred[which(surv_pred$dormancy == "Anthropogenic"), ], 
                                                                                          se.fit = T)$se.fit)

surv_pred$mean_obs[which(surv_pred$dormancy == "Anthropogenic" & 
                           surv_pred$vr == "Survival probability")] = apply(surv_pred[which(surv_pred$dormancy == "Anthropogenic" & 
                                                                                              surv_pred$vr == "Survival probability"), ], 1,
                                                                            function(x) mean(droso_anthropogenic$surv[which(!is.na(droso_anthropogenic$surv) &
                                                                                                                          droso_anthropogenic$time == as.character(x[6]) &
                                                                                                                          droso_anthropogenic$site == as.character(x[7]))], na.rm = T))


# Growth
growth_pred = expand.grid(size = 0,
                          fallR = 0,
                          summerT = 0,
                          abLarge = 0,
                          TSFcont = 0,
                          time = levels(droplevels(droso_anthropogenic$time[which(droso_anthropogenic$time != 2022)])),
                          site = c(levels(droso_anthropogenic$site), levels(droso_natural$site)))

growth_pred$vr = "Growth"
growth_pred$dormancy = "Anthropogenic"
growth_pred$dormancy[which(growth_pred$site %in% levels(droso_natural$site))] = "Natural"

growth_pred$mean_pred[which(growth_pred$dormancy == "Natural")] = (predict(growth_natural, 
                                                                           newdata = growth_pred[which(growth_pred$dormancy == "Natural"), ], 
                                                                           se.fit = T)$fit)
growth_pred$upr[which(growth_pred$dormancy == "Natural")] = (predict(growth_natural,
                                                                          newdata = growth_pred[which(growth_pred$dormancy == "Natural"), ], 
                                                                          se.fit = T)$fit + 
                                                                         1.96 * predict(growth_natural, 
                                                                                        newdata = growth_pred[which(growth_pred$dormancy == "Natural"), ], 
                                                                                        se.fit = T)$se.fit)
growth_pred$lwr[which(growth_pred$dormancy == "Natural")] = (predict(growth_natural, 
                                                                               newdata = growth_pred[which(growth_pred$dormancy == "Natural"), ], 
                                                                               se.fit = T)$fit -
                                                                         1.96 * predict(growth_natural, 
                                                                                        newdata = growth_pred[which(growth_pred$dormancy == "Natural"), ], 
                                                                                        se.fit = T)$se.fit)

growth_pred$mean_obs[which(growth_pred$dormancy == "Natural" & 
                             growth_pred$vr == "Growth")] = apply(growth_pred[which(growth_pred$dormancy == "Natural" & 
                                                                                      growth_pred$vr == "Growth"), ], 1,
                                                             function(x) mean(droso_natural$sizeNext[which(!is.na(droso_natural$sizeNext) &
                                                                                                             droso_natural$time == as.character(x[6]) &
                                                                                                             droso_natural$site == as.character(x[7]))], na.rm = T))


growth_pred$mean_pred[which(growth_pred$dormancy == "Anthropogenic")] = (predict(growth_anthropogenic, 
                                                                                  newdata = growth_pred[which(growth_pred$dormancy == "Anthropogenic"), ], 
                                                                                  se.fit = T)$fit)
growth_pred$upr[which(growth_pred$dormancy == "Anthropogenic")] = (predict(growth_anthropogenic, 
                                                                                 newdata = growth_pred[which(growth_pred$dormancy == "Anthropogenic"), ], 
                                                                                 se.fit = T)$fit + 
                                                                           1.96 * predict(growth_anthropogenic, 
                                                                                          newdata = growth_pred[which(growth_pred$dormancy == "Anthropogenic"), ], 
                                                                                          se.fit = T)$se.fit)
growth_pred$lwr[which(growth_pred$dormancy == "Anthropogenic")] = (predict(growth_anthropogenic, 
                                                                                 newdata = growth_pred[which(growth_pred$dormancy == "Anthropogenic"), ], 
                                                                                 se.fit = T)$fit -
                                                                           1.96 * predict(growth_anthropogenic, 
                                                                                          newdata = growth_pred[which(growth_pred$dormancy == "Anthropogenic"), ], 
                                                                                          se.fit = T)$se.fit)

growth_pred$mean_obs[which(growth_pred$dormancy == "Anthropogenic" & 
                             growth_pred$vr == "Growth")] = apply(growth_pred[which(growth_pred$dormancy == "Anthropogenic" & 
                                                                                      growth_pred$vr == "Growth"), ], 1,
                                                                  function(x) mean(droso_anthropogenic$sizeNext[which(!is.na(droso_anthropogenic$sizeNext) &
                                                                                                                    droso_anthropogenic$time == as.character(x[6]) &
                                                                                                                    droso_anthropogenic$site == as.character(x[7]))], na.rm = T))


# Flowering
flowering_pred = expand.grid(size = 0,
                             prevwinterR = 0,
                             prevfallR = 0,
                             prevwinterT = 0,
                             abLarge = 0,
                             TSFcont = 0,
                             time = levels(droplevels(droso_anthropogenic$time[which(droso_anthropogenic$time != 2022)])),
                             site = c(levels(droso_anthropogenic$site), levels(droso_natural$site)))

flowering_pred$vr = "Flowering probability"
flowering_pred$dormancy = "Anthropogenic"
flowering_pred$dormancy[which(flowering_pred$site %in% levels(droso_natural$site))] = "Natural"

flowering_pred$mean_pred[which(flowering_pred$dormancy == "Natural")] = plogis(predict(flowering_natural, 
                                                                       newdata = flowering_pred[which(flowering_pred$dormancy == "Natural"), ], 
                                                                       se.fit = T)$fit)
flowering_pred$upr[which(flowering_pred$dormancy == "Natural")] = plogis(predict(flowering_natural,
                                                                     newdata = flowering_pred[which(flowering_pred$dormancy == "Natural"), ], 
                                                                     se.fit = T)$fit + 
                                                               1.96 * predict(flowering_natural, 
                                                                              newdata = flowering_pred[which(flowering_pred$dormancy == "Natural"), ], 
                                                                              se.fit = T)$se.fit)
flowering_pred$lwr[which(flowering_pred$dormancy == "Natural")] = plogis(predict(flowering_natural, 
                                                                     newdata = flowering_pred[which(flowering_pred$dormancy == "Natural"), ], 
                                                                     se.fit = T)$fit -
                                                               1.96 * predict(flowering_natural, 
                                                                              newdata = flowering_pred[which(flowering_pred$dormancy == "Natural"), ], 
                                                                              se.fit = T)$se.fit)

flowering_pred$mean_obs[which(flowering_pred$dormancy == "Natural" & 
                                flowering_pred$vr == "Flowering probability")] = apply(flowering_pred[which(flowering_pred$dormancy == "Natural" &
                                                                                                              flowering_pred$vr == "Flowering probability"), ], 1,
                                                                function(x) mean(droso_natural$fl[which(!is.na(droso_natural$fl) &
                                                                                                          droso_natural$size_unscaled >= min(droso_natural$size_unscaled[which(droso_natural$fl == 1)]) &
                                                                                                          droso_natural$time == as.character(x[7]) &
                                                                                                          droso_natural$site == as.character(x[8]))], na.rm = T))


flowering_pred$mean_pred[which(flowering_pred$dormancy == "Anthropogenic")] = plogis(predict(flowering_anthropogenic, 
                                                                         newdata = flowering_pred[which(flowering_pred$dormancy == "Anthropogenic"), ], 
                                                                         se.fit = T)$fit)
flowering_pred$upr[which(flowering_pred$dormancy == "Anthropogenic")] = plogis(predict(flowering_anthropogenic, 
                                                                       newdata = flowering_pred[which(flowering_pred$dormancy == "Anthropogenic"), ], 
                                                                       se.fit = T)$fit + 
                                                                 1.96 * predict(flowering_anthropogenic, 
                                                                                newdata = flowering_pred[which(flowering_pred$dormancy == "Anthropogenic"), ], 
                                                                                se.fit = T)$se.fit)
flowering_pred$lwr[which(flowering_pred$dormancy == "Anthropogenic")] = plogis(predict(flowering_anthropogenic, 
                                                                       newdata = flowering_pred[which(flowering_pred$dormancy == "Anthropogenic"), ], 
                                                                       se.fit = T)$fit -
                                                                 1.96 * predict(flowering_anthropogenic, 
                                                                                newdata = flowering_pred[which(flowering_pred$dormancy == "Anthropogenic"), ], 
                                                                                se.fit = T)$se.fit)

flowering_pred$mean_obs[which(flowering_pred$dormancy == "Anthropogenic" & 
                                flowering_pred$vr == "Flowering probability")] = apply(flowering_pred[which(flowering_pred$dormancy == "Anthropogenic" &
                                                                                                              flowering_pred$vr == "Flowering probability"), ], 1,
                                                                                       function(x) mean(droso_anthropogenic$fl[which(!is.na(droso_anthropogenic$fl) &
                                                                                                                                   droso_anthropogenic$size_unscaled >= min(droso_anthropogenic$size_unscaled[which(droso_anthropogenic$fl == 1)]) &
                                                                                                                                   droso_anthropogenic$time == as.character(x[7]) &
                                                                                                                                   droso_anthropogenic$site == as.character(x[8]))], na.rm = T))


# Number of flowers
nbFlowers_pred = expand.grid(size = 0,
                             prevfallR = 0,
                             prevwinterT = 0,
                             TSFcont = 0,
                             time = levels(droplevels(droso_anthropogenic$time[which(droso_anthropogenic$time != 2022)])),
                             site = c(levels(droso_anthropogenic$site), levels(droso_natural$site)))

nbFlowers_pred$vr = "Number of flowers"
nbFlowers_pred$dormancy = "Anthropogenic"
nbFlowers_pred$dormancy[which(nbFlowers_pred$site %in% levels(droso_natural$site))] = "Natural"

nbFlowers_pred$size[which(nbFlowers_pred$dormancy == "Natural")] = (mean_size_repro_natural - mean(droso_natural$size_unscaled, na.rm = T)) / (2 * sd(droso_natural$size_unscaled, na.rm = T))
nbFlowers_pred$size[which(nbFlowers_pred$dormancy == "Anthropogenic")] = (mean_size_repro_anthropogenic - mean(droso_anthropogenic$size_unscaled, na.rm = T)) / (2 * sd(droso_anthropogenic$size_unscaled, na.rm = T))

nbFlowers_pred$mean_pred[which(nbFlowers_pred$dormancy == "Natural")] = exp(predict(nbFlow_natural, 
                                                                                      newdata = nbFlowers_pred[which(nbFlowers_pred$dormancy == "Natural"), ], 
                                                                                      se.fit = T)$fit)
nbFlowers_pred$upr[which(nbFlowers_pred$dormancy == "Natural")] = exp(predict(nbFlow_natural,
                                                                                    newdata = nbFlowers_pred[which(nbFlowers_pred$dormancy == "Natural"), ], 
                                                                                    se.fit = T)$fit + 
                                                                              1.96 * predict(nbFlow_natural, 
                                                                                             newdata = nbFlowers_pred[which(nbFlowers_pred$dormancy == "Natural"), ], 
                                                                                             se.fit = T)$se.fit)
nbFlowers_pred$lwr[which(nbFlowers_pred$dormancy == "Natural")] = exp(predict(nbFlow_natural, 
                                                                                    newdata = nbFlowers_pred[which(nbFlowers_pred$dormancy == "Natural"), ], 
                                                                                    se.fit = T)$fit -
                                                                              1.96 * predict(nbFlow_natural, 
                                                                                             newdata = nbFlowers_pred[which(nbFlowers_pred$dormancy == "Natural"), ], 
                                                                                             se.fit = T)$se.fit)

nbFlowers_pred$mean_obs[which(nbFlowers_pred$dormancy == "Natural" & 
                                nbFlowers_pred$vr == "Number of flowers")] = apply(nbFlowers_pred[which(nbFlowers_pred$dormancy == "Natural" &
                                                                                                          nbFlowers_pred$vr == "Number of flowers"), ], 1,
                                                                        function(x) mean(droso_natural$nbFlow[which(!is.na(droso_natural$fl) &
                                                                                                                      !is.na(droso_natural$nbFlow) &
                                                                                                                      droso_natural$time == as.character(x[5]) &
                                                                                                                      droso_natural$site == as.character(x[6]))], na.rm = T))


nbFlowers_pred$mean_pred[which(nbFlowers_pred$dormancy == "Anthropogenic")] = exp(predict(nbFlow_anthropogenic, 
                                                                                        newdata = nbFlowers_pred[which(nbFlowers_pred$dormancy == "Anthropogenic"), ], 
                                                                                        se.fit = T)$fit)
nbFlowers_pred$upr[which(nbFlowers_pred$dormancy == "Anthropogenic")] = exp(predict(nbFlow_anthropogenic, 
                                                                                      newdata = nbFlowers_pred[which(nbFlowers_pred$dormancy == "Anthropogenic"), ], 
                                                                                      se.fit = T)$fit + 
                                                                                1.96 * predict(nbFlow_anthropogenic, 
                                                                                               newdata = nbFlowers_pred[which(nbFlowers_pred$dormancy == "Anthropogenic"), ], 
                                                                                               se.fit = T)$se.fit)
nbFlowers_pred$lwr[which(nbFlowers_pred$dormancy == "Anthropogenic")] = exp(predict(nbFlow_anthropogenic, 
                                                                                      newdata = nbFlowers_pred[which(nbFlowers_pred$dormancy == "Anthropogenic"), ], 
                                                                                      se.fit = T)$fit -
                                                                                1.96 * predict(nbFlow_anthropogenic, 
                                                                                               newdata = nbFlowers_pred[which(nbFlowers_pred$dormancy == "Anthropogenic"), ], 
                                                                                               se.fit = T)$se.fit)

nbFlowers_pred$mean_obs[which(nbFlowers_pred$dormancy == "Anthropogenic" & 
                                nbFlowers_pred$vr == "Number of flowers")] = apply(nbFlowers_pred[which(nbFlowers_pred$dormancy == "Anthropogenic" &
                                                                                                          nbFlowers_pred$vr == "Number of flowers"), ], 1,
                                                                                   function(x) mean(droso_anthropogenic$nbFlow[which(!is.na(droso_anthropogenic$fl) &
                                                                                                                                 !is.na(droso_anthropogenic$nbFlow) &
                                                                                                                                   droso_anthropogenic$time == as.character(x[5]) &
                                                                                                                                   droso_anthropogenic$site == as.character(x[6]))], na.rm = T))


# Seedling size
seedlingSize_pred = expand.grid(prevfallR = 0,
                                prevwinterT = 0,
                                abLarge = 0,
                                TSFcont = 0,
                                time = levels(droplevels(droso_anthropogenic$time[which(droso_anthropogenic$time != 2022)])),
                                site = c(levels(droso_anthropogenic$site), levels(droso_natural$site)))

seedlingSize_pred$vr = "Seedling size"
seedlingSize_pred$dormancy = "Anthropogenic"
seedlingSize_pred$dormancy[which(seedlingSize_pred$site %in% levels(droso_natural$site))] = "Natural"

seedlingSize_pred$mean_pred[which(seedlingSize_pred$dormancy == "Natural")] = (predict(seedlingSize_natural, 
                                                                                newdata = seedlingSize_pred[which(seedlingSize_pred$dormancy == "Natural"), ], 
                                                                                se.fit = T)$fit)
seedlingSize_pred$upr[which(seedlingSize_pred$dormancy == "Natural")] = (predict(seedlingSize_natural,
                                                                              newdata = seedlingSize_pred[which(seedlingSize_pred$dormancy == "Natural"), ], 
                                                                              se.fit = T)$fit + 
                                                                        1.96 * predict(seedlingSize_natural, 
                                                                                       newdata = seedlingSize_pred[which(seedlingSize_pred$dormancy == "Natural"), ], 
                                                                                       se.fit = T)$se.fit)
seedlingSize_pred$lwr[which(seedlingSize_pred$dormancy == "Natural")] = (predict(seedlingSize_natural, 
                                                                              newdata = seedlingSize_pred[which(seedlingSize_pred$dormancy == "Natural"), ], 
                                                                              se.fit = T)$fit -
                                                                        1.96 * predict(seedlingSize_natural, 
                                                                                       newdata = seedlingSize_pred[which(seedlingSize_pred$dormancy == "Natural"), ], 
                                                                                       se.fit = T)$se.fit)

seedlingSize_pred$mean_obs[which(seedlingSize_pred$dormancy == "Natural" & 
                                   seedlingSize_pred$vr == "Seedling size")] = apply(seedlingSize_pred[which(seedlingSize_pred$dormancy == "Natural" & 
                                                                                                               seedlingSize_pred$vr == "Seedling size"), ], 1,
                                                                    function(x) mean(droso_natural$size_unscaled[which(!is.na(droso_natural$size_unscaled) &
                                                                                                                         droso_natural$stage == "SD" &
                                                                                                                         droso_natural$time == as.character(x[5]) &
                                                                                                                         droso_natural$site == as.character(x[6]))], na.rm = T))


seedlingSize_pred$mean_pred[which(seedlingSize_pred$dormancy == "Anthropogenic")] = (predict(seedlingSize_anthropogenic, 
                                                                                  newdata = seedlingSize_pred[which(seedlingSize_pred$dormancy == "Anthropogenic"), ], 
                                                                                  se.fit = T)$fit)
seedlingSize_pred$upr[which(seedlingSize_pred$dormancy == "Anthropogenic")] = (predict(seedlingSize_anthropogenic, 
                                                                                newdata = seedlingSize_pred[which(seedlingSize_pred$dormancy == "Anthropogenic"), ], 
                                                                                se.fit = T)$fit + 
                                                                          1.96 * predict(seedlingSize_anthropogenic, 
                                                                                         newdata = seedlingSize_pred[which(seedlingSize_pred$dormancy == "Anthropogenic"), ], 
                                                                                         se.fit = T)$se.fit)
seedlingSize_pred$lwr[which(seedlingSize_pred$dormancy == "Anthropogenic")] = (predict(seedlingSize_anthropogenic, 
                                                                                newdata = seedlingSize_pred[which(seedlingSize_pred$dormancy == "Anthropogenic"), ], 
                                                                                se.fit = T)$fit -
                                                                          1.96 * predict(seedlingSize_anthropogenic, 
                                                                                         newdata = seedlingSize_pred[which(seedlingSize_pred$dormancy == "Anthropogenic"), ], 
                                                                                         se.fit = T)$se.fit)

seedlingSize_pred$mean_obs[which(seedlingSize_pred$dormancy == "Anthropogenic" & 
                                   seedlingSize_pred$vr == "Seedling size")] = apply(seedlingSize_pred[which(seedlingSize_pred$dormancy == "Anthropogenic" & 
                                                                                                               seedlingSize_pred$vr == "Seedling size"), ], 1,
                                                                                     function(x) mean(droso_anthropogenic$size_unscaled[which(!is.na(droso_anthropogenic$size_unscaled) &
                                                                                                                                            droso_anthropogenic$stage == "SD" &
                                                                                                                                            droso_anthropogenic$time == as.character(x[5]) &
                                                                                                                                            droso_anthropogenic$site == as.character(x[6]))], na.rm = T))


avg_vr_sites = rbind(surv_pred[, c("dormancy", "time", "site", "vr", "mean_pred", "lwr", "upr", "mean_obs")],
                     growth_pred[, c("dormancy", "time", "site", "vr", "mean_pred", "lwr", "upr", "mean_obs")],
                     flowering_pred[, c("dormancy", "time", "site", "vr", "mean_pred", "lwr", "upr", "mean_obs")],
                     nbFlowers_pred[, c("dormancy", "time", "site", "vr", "mean_pred", "lwr", "upr", "mean_obs")],
                     seedlingSize_pred[, c("dormancy", "time", "site", "vr", "mean_pred", "lwr", "upr", "mean_obs")])

avg_vr_sites$site = factor(avg_vr_sites$site,
                           labels = c("Bujeo", "Montera\ndel\nTorero", "Prisoneros",
                                      "Sierra del\nRetín\nDisturbed", "Sierra\nCarbonera\nDisturbed",
                                      "Sierra\nCarbonera\nYoung", "Sierra del\nRetín\nYoung",
                                      "Vertedero"))

avg_vr_sites$dormancy = factor(avg_vr_sites$dormancy, 
                               levels = c("Natural", "Anthropogenic"))

avg_vr_sites_CI = aggregate(mean_pred ~ site + dormancy + vr, 
                            data = avg_vr_sites,
                            FUN = function(x) quantile(x, probs = c(0.025, 0.975), na.rm = T))


png("Output/Plots/AppendixS1_Figure1.png", 
    width = 15, height = 18, units = "cm",
    res = 600,
    pointsize = 12, bg = "white")

ggplot(avg_vr_sites, aes(y = mean_obs, x = site)) +
  facet_grid(vr ~ dormancy, scales = "free") +
  geom_boxplot(aes(y = mean_pred),
               ymin = avg_vr_sites_CI$mean_pred[, 1],
               ymax = avg_vr_sites_CI$mean_pred[, 2],
               outlier.size = 0.1, size = 0.1,
               linewidth = 0.4,
               fill = alpha("grey", 0.5)) + 
  geom_point(aes(y = mean_obs),
             position = position_quasirandom(width = 0.1),
             size = 0.8, alpha = 0.5, shape = 16,
             colour = "lightsalmon") +
  stat_summary(fun = mean, aes(shape = "Mean"), geom = "point", shape = 17, size = 1, position = position_dodge(width = 0.75)) +
  xlab("Population") +
  ylab("Vital-rate value") +
  theme_general()

dev.off()




###########################################################################
#
# 5. Plotting predictions - Appendix - Among-site differences in effects ----
#
###########################################################################

## 5.1. Survival as a function of rainfall (other covariates at average value) ----
# ----------------------------------------------------------------------------

surv_pred_noRE = expand.grid(size = 0,
                             fallR = (fallR_values_anthropogenic - mean(fallR_timeSeries_anthropogenic$fallR_unscaled)) / (2 * sd(fallR_timeSeries_anthropogenic$fallR_unscaled)),
                             summerT = 0,
                             abLarge = 0,
                             dormancy = "Anthropogenic",
                             site = unique(droso_anthropogenic$site))


surv_pred_noRE$surv = NA
surv_pred_noRE$upr = NA
surv_pred_noRE$lwr = NA


surv_pred_noRE$surv[which(surv_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                     newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                     se.fit = T,
                                                                                     exclude = c("s(time)"), 
                                                                                     newdata.guaranteed = T)$fit)
surv_pred_noRE$upr[which(surv_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)"), 
                                                                                    newdata.guaranteed = T)$fit + 
                                                                              1.96 * predict(surv_anthropogenic, 
                                                                                             newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)"), 
                                                                                             newdata.guaranteed = T)$se.fit)
surv_pred_noRE$lwr[which(surv_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)"), 
                                                                                    newdata.guaranteed = T)$fit - 
                                                                              1.96 * predict(surv_anthropogenic, 
                                                                                             newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)"), 
                                                                                             newdata.guaranteed = T)$se.fit)

surv_pred_noRE$fallR_unscaled[which(surv_pred_noRE$dormancy == "Anthropogenic")] = surv_pred_noRE$fallR[which(surv_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(fallR_timeSeries_anthropogenic$fallR_unscaled)) + mean(fallR_timeSeries_anthropogenic$fallR_unscaled)

surv_plot_rain_sites = ggplot(surv_pred_noRE, aes(fallR_unscaled, surv)) +
  facet_wrap(~ site, 
             labeller = labeller(site = c("Bujeo" = "Bujeo",
                                          "MonteraTorero" = "Montera del\nTorero",
                                          "Prisoneros" = "Prisoneros",
                                          "Retin" = "Sierra del\nRetín\nDisturbed",
                                          "SCarbDist" = "Sierra\nCarbonera\nDisturbed"))) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, col = NA) +
  xlab("Next fall cumulative rainfall\n(Sep-Nov) (in mm)") +
  ylab("Survival probability") +
  theme_general()




## 5.2. Growth as a function of temperature (other covariates at average value) ----
# -----------------------------------------------------------------------------

growth_pred_noRE = expand.grid(size = 0,
                               summerT = (summerT_values_anthropogenic - mean(summerT_timeSeries_anthropogenic$summerT_unscaled)) / (2 * sd(summerT_timeSeries_anthropogenic$summerT_unscaled)),
                               abLarge = 0,
                               dormancy = "Anthropogenic",
                               site = unique(droso_anthropogenic$site))


growth_pred_noRE$sizeNext = NA
growth_pred_noRE$upr = NA
growth_pred_noRE$lwr = NA


growth_pred_noRE$sizeNext[which(growth_pred_noRE$dormancy == "Anthropogenic")] = (predict(growth_anthropogenic, 
                                                                                       newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                       se.fit = T,
                                                                                       exclude = c("s(time)"), 
                                                                                       newdata.guaranteed = T)$fit)
growth_pred_noRE$upr[which(growth_pred_noRE$dormancy == "Anthropogenic")] = (predict(growth_anthropogenic, 
                                                                                      newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                      se.fit = T,
                                                                                      exclude = c("s(time)"), 
                                                                                      newdata.guaranteed = T)$fit + 
                                                                                1.96 * predict(growth_anthropogenic, 
                                                                                               newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                               se.fit = T,
                                                                                               exclude = c("s(time)"), 
                                                                                               newdata.guaranteed = T)$se.fit)
growth_pred_noRE$lwr[which(growth_pred_noRE$dormancy == "Anthropogenic")] = (predict(growth_anthropogenic, 
                                                                                      newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                      se.fit = T,
                                                                                      exclude = c("s(time)"), 
                                                                                      newdata.guaranteed = T)$fit - 
                                                                                1.96 * predict(growth_anthropogenic, 
                                                                                               newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                               se.fit = T,
                                                                                               exclude = c("s(time)"), 
                                                                                               newdata.guaranteed = T)$se.fit)

growth_pred_noRE$summerT_unscaled[which(growth_pred_noRE$dormancy == "Anthropogenic")] = growth_pred_noRE$summerT[which(growth_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(summerT_timeSeries_anthropogenic$summerT_unscaled)) + mean(summerT_timeSeries_anthropogenic$summerT_unscaled)

growth_plot_temp_sites = ggplot(growth_pred_noRE, aes(summerT_unscaled, sizeNext)) +
  facet_wrap(~ site, 
             labeller = labeller(site = c("Bujeo" = "Bujeo",
                                          "MonteraTorero" = "Montera del\nTorero",
                                          "Prisoneros" = "Prisoneros",
                                          "Retin" = "Sierra del\nRetín\nDisturbed",
                                          "SCarbDist" = "Sierra\nCarbonera\nDisturbed"))) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, col = NA) +
  xlab("Next summer mean max. daily temperature\n(May-Sep) (in ºC)") +
  ylab("Size next year") +
  theme_general()


## 5.3. Flowering as a function of rainfall (other covariates at average value) ----
# -----------------------------------------------------------------------------

flowering_pred_noRE = expand.grid(size = 0,
                                  prevwinterR = (prevwinterR_values_anthropogenic - mean(prevwinterR_timeSeries_anthropogenic$prevwinterR_unscaled)) / (2 * sd(prevwinterR_timeSeries_anthropogenic$prevwinterR_unscaled)),
                                  abLarge = 0,
                                  dormancy = "Anthropogenic",
                                  site = unique(droso_anthropogenic$site))


flowering_pred_noRE$flowering = NA
flowering_pred_noRE$upr = NA
flowering_pred_noRE$lwr = NA


flowering_pred_noRE$flowering[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(flowering_anthropogenic, 
                                                                                      newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                      se.fit = T,
                                                                                      exclude = c("s(time)"), 
                                                                                      newdata.guaranteed = T)$fit)
flowering_pred_noRE$upr[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(flowering_anthropogenic, 
                                                                                 newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                 se.fit = T,
                                                                                 exclude = c("s(time)"), 
                                                                                 newdata.guaranteed = T)$fit + 
                                                                           1.96 * predict(flowering_anthropogenic, 
                                                                                          newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                          se.fit = T,
                                                                                          exclude = c("s(time)"), 
                                                                                          newdata.guaranteed = T)$se.fit)
flowering_pred_noRE$lwr[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(flowering_anthropogenic, 
                                                                                 newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                 se.fit = T,
                                                                                 exclude = c("s(time)"), 
                                                                                 newdata.guaranteed = T)$fit - 
                                                                           1.96 * predict(flowering_anthropogenic, 
                                                                                          newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                          se.fit = T,
                                                                                          exclude = c("s(time)"), 
                                                                                          newdata.guaranteed = T)$se.fit)

flowering_pred_noRE$prevwinterR_unscaled[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = flowering_pred_noRE$prevwinterR[which(flowering_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(prevwinterR_timeSeries_anthropogenic$prevwinterR_unscaled)) + mean(prevwinterR_timeSeries_anthropogenic$prevwinterR_unscaled)

flowering_plot_rain_sites = ggplot(flowering_pred_noRE, aes(prevwinterR_unscaled, flowering)) +
  facet_wrap(~ site, 
             labeller = labeller(site = c("Bujeo" = "Bujeo",
                                          "MonteraTorero" = "Montera del\nTorero",
                                          "Prisoneros" = "Prisoneros",
                                          "Retin" = "Sierra del\nRetín\nDisturbed",
                                          "SCarbDist" = "Sierra\nCarbonera\nDisturbed"))) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, col = NA) +
  xlab("Previous winter cumulative rainfall\n(Jan-Apr) (in mm)") +
  ylab("Flowering probability") +
  theme_general()


## 5.4. Number of flowers as a function of rainfall (other covariates at average value) ----
# -------------------------------------------------------------------------------------

nbFlowers_pred_noRE = expand.grid(size = (mean_size_repro_anthropogenic - mean(droso_anthropogenic$size_unscaled, na.rm = T)) / (2 * sd(droso_anthropogenic$size_unscaled, na.rm = T)),
                                  prevfallR = (prevfallR_values_anthropogenic - mean(prevfallR_timeSeries_anthropogenic$prevfallR_unscaled)) / (2 * sd(prevfallR_timeSeries_anthropogenic$prevfallR_unscaled)),
                                  dormancy = "Anthropogenic",
                                  site = unique(droso_anthropogenic$site))


nbFlowers_pred_noRE$nbFlowers = NA
nbFlowers_pred_noRE$upr = NA
nbFlowers_pred_noRE$lwr = NA


nbFlowers_pred_noRE$nbFlowers[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic")] = exp(predict(nbFlow_anthropogenic, 
                                                                                                      newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                      se.fit = T,
                                                                                                      exclude = c("s(time)"), 
                                                                                                      newdata.guaranteed = T)$fit)
nbFlowers_pred_noRE$upr[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic")] = exp(predict(nbFlow_anthropogenic, 
                                                                                                newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                se.fit = T,
                                                                                                exclude = c("s(time)"), 
                                                                                                newdata.guaranteed = T)$fit + 
                                                                                          1.96 * predict(nbFlow_anthropogenic, 
                                                                                                         newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                         se.fit = T,
                                                                                                         exclude = c("s(time)"), 
                                                                                                         newdata.guaranteed = T)$se.fit)
nbFlowers_pred_noRE$lwr[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic")] = exp(predict(nbFlow_anthropogenic, 
                                                                                                newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                se.fit = T,
                                                                                                exclude = c("s(time)"), 
                                                                                                newdata.guaranteed = T)$fit - 
                                                                                          1.96 * predict(nbFlow_anthropogenic, 
                                                                                                         newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                         se.fit = T,
                                                                                                         exclude = c("s(time)"), 
                                                                                                         newdata.guaranteed = T)$se.fit)

nbFlowers_pred_noRE$prevfallR_unscaled[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic")] = nbFlowers_pred_noRE$prevfallR[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(prevfallR_timeSeries_anthropogenic$prevfallR_unscaled)) + mean(prevfallR_timeSeries_anthropogenic$prevfallR_unscaled)

nbFlowers_plot_rain_sites = ggplot(nbFlowers_pred_noRE, aes(prevfallR_unscaled, nbFlowers)) +
  facet_wrap(~ site, 
             labeller = labeller(site = c("Bujeo" = "Bujeo",
                                          "MonteraTorero" = "Montera del\nTorero",
                                          "Prisoneros" = "Prisoneros",
                                          "Retin" = "Sierra del\nRetín\nDisturbed",
                                          "SCarbDist" = "Sierra\nCarbonera\nDisturbed"))) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, col = NA) +
  xlab("Previous fall cumulative rainfall\n(Sep-Nov) (in mm)") +
  ylab("Number of flowers per individual") +
  theme_general()


## 5.5. Seedling size as a function of temperature (other covariates at average value) ----
# ------------------------------------------------------------------------------------

seedlingSize_pred_noRE = expand.grid(prevwinterT = (prevwinterT_values_anthropogenic - mean(prevwinterT_timeSeries_anthropogenic$prevwinterT_unscaled)) / (2 * sd(prevwinterT_timeSeries_anthropogenic$prevwinterT_unscaled)),
                                     abLarge = 0,
                                     dormancy = "Anthropogenic",
                                     site = unique(droso_anthropogenic$site))


seedlingSize_pred_noRE$seedlingSize = NA
seedlingSize_pred_noRE$upr = NA
seedlingSize_pred_noRE$lwr = NA


seedlingSize_pred_noRE$seedlingSize[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic")] = (predict(seedlingSize_anthropogenic, 
                                                                                                newdata = seedlingSize_pred_noRE[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                se.fit = T,
                                                                                                exclude = c("s(time)"), 
                                                                                                newdata.guaranteed = T)$fit)
seedlingSize_pred_noRE$upr[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic")] = (predict(seedlingSize_anthropogenic, 
                                                                                          newdata = seedlingSize_pred_noRE[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                          se.fit = T,
                                                                                          exclude = c("s(time)"), 
                                                                                          newdata.guaranteed = T)$fit + 
                                                                                    1.96 * predict(seedlingSize_anthropogenic, 
                                                                                                   newdata = seedlingSize_pred_noRE[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                   se.fit = T,
                                                                                                   exclude = c("s(time)"), 
                                                                                                   newdata.guaranteed = T)$se.fit)
seedlingSize_pred_noRE$lwr[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic")] = (predict(seedlingSize_anthropogenic, 
                                                                                          newdata = seedlingSize_pred_noRE[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                          se.fit = T,
                                                                                          exclude = c("s(time)"), 
                                                                                          newdata.guaranteed = T)$fit - 
                                                                                    1.96 * predict(seedlingSize_anthropogenic, 
                                                                                                   newdata = seedlingSize_pred_noRE[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                   se.fit = T,
                                                                                                   exclude = c("s(time)"), 
                                                                                                   newdata.guaranteed = T)$se.fit)

seedlingSize_pred_noRE$prevwinterT_unscaled[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic")] = seedlingSize_pred_noRE$prevwinterT[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(prevwinterT_timeSeries_anthropogenic$prevwinterT_unscaled)) + mean(prevwinterT_timeSeries_anthropogenic$prevwinterT_unscaled)

seedlingSize_plot_temp_sites = ggplot(seedlingSize_pred_noRE, aes(prevwinterT_unscaled, seedlingSize)) +
  facet_wrap(~ site, 
             labeller = labeller(site = c("Bujeo" = "Bujeo",
                                          "MonteraTorero" = "Montera del\nTorero",
                                          "Prisoneros" = "Prisoneros",
                                          "Retin" = "Sierra del\nRetín\nDisturbed",
                                          "SCarbDist" = "Sierra\nCarbonera\nDisturbed"))) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, col = NA) +
  xlab("Previous winter mean max. daily temperature\n(Jan-Apr) (in mm)") +
  ylab("Seedling size") +
  theme_general()


## 5.6. Final plot ----
# ----------------

png(filename = "Output/Plots/AppendixS1_Figure2.png", 
    width = 16, 
    height = 16, 
    units = "cm", 
    bg = "white", 
    res = 600, 
    type = "cairo")

surv_plot_rain_sites + growth_plot_temp_sites + nbFlowers_plot_rain_sites + seedlingSize_plot_temp_sites +
  plot_layout(ncol = 2) + 
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") &
  theme(plot.tag = element_text(size = 10))

dev.off()




###########################################################################
#
# 6. Plotting predictions - Appendix - Density effects ----
#
###########################################################################

## 6.1. Seedling size as a function of rainfall per density (other covariates at average value) ----
# ---------------------------------------------------------------------------------------------

seedlingSize_pred_noRE = rbind(expand.grid(prevwinterT = (prevwinterT_values_anthropogenic - mean(prevwinterT_timeSeries_anthropogenic$prevwinterT_unscaled)) / (2 * sd(prevwinterT_timeSeries_anthropogenic$prevwinterT_unscaled)),
                                           abLarge = (c(2, 6, 10) - mean(yearly_density_per_square_anthropogenic$abLarge_unscaled)) / (2 * sd(yearly_density_per_square_anthropogenic$abLarge_unscaled)),
                                           TSFcont = 0,
                                           dormancy = "Anthropogenic"),
                               expand.grid(prevwinterT = (prevwinterT_values_natural - mean(prevwinterT_timeSeries_natural$prevwinterT_unscaled)) / (2 * sd(prevwinterT_timeSeries_natural$prevwinterT_unscaled)),
                                           abLarge = (c(2, 6, 10) - mean(yearly_density_per_square_natural$abLarge_unscaled)) / (2 * sd(yearly_density_per_square_natural$abLarge_unscaled)),
                                           TSFcont = 0,
                                           dormancy = "Natural"))


seedlingSize_pred_noRE$seedlingSize = NA
seedlingSize_pred_noRE$upr = NA
seedlingSize_pred_noRE$lwr = NA


seedlingSize_pred_noRE$seedlingSize[which(seedlingSize_pred_noRE$dormancy == "Natural")] = (predict(seedlingSize_natural, 
                                                                                                    newdata = seedlingSize_pred_noRE[which(seedlingSize_pred_noRE$dormancy == "Natural"), ], 
                                                                                                    se.fit = T,
                                                                                                    exclude = c("s(time)", "s(site)", "s(abLarge,site)", "s(TSFcont,site)"), 
                                                                                                    newdata.guaranteed = T)$fit)
seedlingSize_pred_noRE$upr[which(seedlingSize_pred_noRE$dormancy == "Natural")] = (predict(seedlingSize_natural, 
                                                                                              newdata = seedlingSize_pred_noRE[which(seedlingSize_pred_noRE$dormancy == "Natural"), ], 
                                                                                              se.fit = T,
                                                                                              exclude = c("s(time)", "s(site)", "s(abLarge,site)", "s(TSFcont,site)"), 
                                                                                              newdata.guaranteed = T)$fit + 
                                                                                        1.96 * predict(seedlingSize_natural, 
                                                                                                       newdata = seedlingSize_pred_noRE[which(seedlingSize_pred_noRE$dormancy == "Natural"), ], 
                                                                                                       se.fit = T,
                                                                                                       exclude = c("s(time)", "s(site)", "s(abLarge,site)", "s(TSFcont,site)"), 
                                                                                                       newdata.guaranteed = T)$se.fit)
seedlingSize_pred_noRE$lwr[which(seedlingSize_pred_noRE$dormancy == "Natural")] = (predict(seedlingSize_natural, 
                                                                                              newdata = seedlingSize_pred_noRE[which(seedlingSize_pred_noRE$dormancy == "Natural"), ], 
                                                                                              se.fit = T,
                                                                                              exclude = c("s(time)", "s(site)", "s(abLarge,site)", "s(TSFcont,site)"), 
                                                                                              newdata.guaranteed = T)$fit - 
                                                                                        1.96 * predict(seedlingSize_natural, 
                                                                                                       newdata = seedlingSize_pred_noRE[which(seedlingSize_pred_noRE$dormancy == "Natural"), ], 
                                                                                                       se.fit = T,
                                                                                                       exclude = c("s(time)", "s(site)", "s(abLarge,site)", "s(TSFcont,site)"), 
                                                                                                       newdata.guaranteed = T)$se.fit)

seedlingSize_pred_noRE$seedlingSize[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic")] = (predict(seedlingSize_anthropogenic, 
                                                                                                      newdata = seedlingSize_pred_noRE[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                      se.fit = T,
                                                                                                      exclude = c("s(time)", "s(site)", "s(prevwinterT,site)", "s(abLarge,site)"), 
                                                                                                      newdata.guaranteed = T)$fit)
seedlingSize_pred_noRE$upr[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic")] = (predict(seedlingSize_anthropogenic, 
                                                                                                newdata = seedlingSize_pred_noRE[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                se.fit = T,
                                                                                                exclude = c("s(time)", "s(site)", "s(prevwinterT,site)", "s(abLarge,site)"), 
                                                                                                newdata.guaranteed = T)$fit + 
                                                                                          1.96 * predict(seedlingSize_anthropogenic, 
                                                                                                         newdata = seedlingSize_pred_noRE[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                         se.fit = T,
                                                                                                         exclude = c("s(time)", "s(site)", "s(prevwinterT,site)", "s(abLarge,site)"), 
                                                                                                         newdata.guaranteed = T)$se.fit)
seedlingSize_pred_noRE$lwr[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic")] = (predict(seedlingSize_anthropogenic, 
                                                                                                newdata = seedlingSize_pred_noRE[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                se.fit = T,
                                                                                                exclude = c("s(time)", "s(site)", "s(prevwinterT,site)", "s(abLarge,site)"), 
                                                                                                newdata.guaranteed = T)$fit - 
                                                                                          1.96 * predict(seedlingSize_anthropogenic, 
                                                                                                         newdata = seedlingSize_pred_noRE[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                         se.fit = T,
                                                                                                         exclude = c("s(time)", "s(site)", "s(prevwinterT,site)", "s(abLarge,site)"), 
                                                                                                         newdata.guaranteed = T)$se.fit)


seedlingSize_pred_noRE$prevwinterT_unscaled[which(seedlingSize_pred_noRE$dormancy == "Natural")] = seedlingSize_pred_noRE$prevwinterT[which(seedlingSize_pred_noRE$dormancy == "Natural")] * (2 * sd(prevwinterT_timeSeries_natural$prevwinterT_unscaled)) + mean(prevwinterT_timeSeries_natural$prevwinterT_unscaled)
seedlingSize_pred_noRE$prevwinterT_unscaled[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic")] = seedlingSize_pred_noRE$prevwinterT[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(prevwinterT_timeSeries_anthropogenic$prevwinterT_unscaled)) + mean(prevwinterT_timeSeries_anthropogenic$prevwinterT_unscaled)

seedlingSize_pred_noRE$abLarge_unscaled[which(seedlingSize_pred_noRE$dormancy == "Natural")] = seedlingSize_pred_noRE$abLarge[which(seedlingSize_pred_noRE$dormancy == "Natural")] * (2 * sd(yearly_density_per_square_natural$abLarge_unscaled)) + mean(yearly_density_per_square_natural$abLarge_unscaled)
seedlingSize_pred_noRE$abLarge_unscaled[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic")] = seedlingSize_pred_noRE$abLarge[which(seedlingSize_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(yearly_density_per_square_anthropogenic$abLarge_unscaled)) + mean(yearly_density_per_square_anthropogenic$abLarge_unscaled)

seedlingSize_pred_noRE$abLarge_unscaled = as.factor(seedlingSize_pred_noRE$abLarge_unscaled)
seedlingSize_pred_noRE$dormancy = factor(seedlingSize_pred_noRE$dormancy,
                                         levels = c("Natural", "Anthropogenic"))


seedlingSize_plot_rain_density = ggplot(seedlingSize_pred_noRE, aes(prevwinterT_unscaled, seedlingSize,
                                                                                                                          col = abLarge_unscaled, group = abLarge_unscaled)) +
  facet_wrap(~ dormancy, scales = "free_x") +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = abLarge_unscaled), alpha = 0.1, col = NA) +
  scale_colour_manual(name = "Aboveground density\nof large individuals\n(in ind./m2)",
                      values = palette.colors(3, palette = "Dark 2")) +
  scale_fill_manual(name = "Aboveground density\nof large individuals\n(in ind./m2)",
                    values = palette.colors(3, palette = "Dark 2")) +
  xlab("Previous winter mean max. daily\ntemperature (Jan-Apr) (in ºC)") +
  ylab("Seedling size") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"))


## 6.2. Flowering as a function of rainfall per density (other covariates at average value) ----
# -----------------------------------------------------------------------------------------

flowering_pred_noRE = rbind(expand.grid(size = 0,
                                        prevfallR = 0,
                                        prevwinterR = (prevwinterR_values_anthropogenic - mean(prevwinterR_timeSeries_anthropogenic$prevwinterR_unscaled)) / (2 * sd(prevwinterR_timeSeries_anthropogenic$prevwinterR_unscaled)),
                                        prevwinterT = 0,
                                        abLarge = (c(2, 6, 10) - mean(yearly_density_per_square_anthropogenic$abLarge_unscaled)) / (2 * sd(yearly_density_per_square_anthropogenic$abLarge_unscaled)),
                                        TSFcont = 0,
                                        dormancy = "Anthropogenic"),
                            expand.grid(size = 0,
                                        prevfallR = (prevfallR_values_natural - mean(prevfallR_timeSeries_natural$prevfallR_unscaled)) / (2 * sd(prevfallR_timeSeries_natural$prevfallR_unscaled)),
                                        prevwinterR = 0,
                                        prevwinterT = 0,
                                        abLarge = (c(2, 6, 10) - mean(yearly_density_per_square_natural$abLarge_unscaled)) / (2 * sd(yearly_density_per_square_natural$abLarge_unscaled)),
                                        TSFcont = 0,
                                        dormancy = "Natural"))


flowering_pred_noRE$flowering = NA
flowering_pred_noRE$upr = NA
flowering_pred_noRE$lwr = NA


flowering_pred_noRE$flowering[which(flowering_pred_noRE$dormancy == "Natural")] = plogis(predict(flowering_natural, 
                                                                                                    newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                                    se.fit = T,
                                                                                                    exclude = c("s(time)", "s(site)"), 
                                                                                                    newdata.guaranteed = T)$fit)
flowering_pred_noRE$upr[which(flowering_pred_noRE$dormancy == "Natural")] = plogis(predict(flowering_natural, 
                                                                                              newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                              se.fit = T,
                                                                                              exclude = c("s(time)", "s(site)"), 
                                                                                              newdata.guaranteed = T)$fit + 
                                                                                        1.96 * predict(flowering_natural, 
                                                                                                       newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                                       se.fit = T,
                                                                                                       exclude = c("s(time)", "s(site)"), 
                                                                                                       newdata.guaranteed = T)$se.fit)
flowering_pred_noRE$lwr[which(flowering_pred_noRE$dormancy == "Natural")] = plogis(predict(flowering_natural, 
                                                                                              newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                              se.fit = T,
                                                                                              exclude = c("s(time)", "s(site)"), 
                                                                                              newdata.guaranteed = T)$fit - 
                                                                                        1.96 * predict(flowering_natural, 
                                                                                                       newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                                       se.fit = T,
                                                                                                       exclude = c("s(time)", "s(site)"), 
                                                                                                       newdata.guaranteed = T)$se.fit)

flowering_pred_noRE$flowering[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(flowering_anthropogenic, 
                                                                                                      newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                      se.fit = T,
                                                                                                      exclude = c("s(time)", "s(site)", "s(size,site)", "s(prevwinterR,site)", "s(abLarge,site)"), 
                                                                                                      newdata.guaranteed = T)$fit)
flowering_pred_noRE$upr[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(flowering_anthropogenic, 
                                                                                                newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                se.fit = T,
                                                                                                exclude = c("s(time)", "s(site)", "s(size,site)", "s(prevwinterR,site)", "s(abLarge,site)"), 
                                                                                                newdata.guaranteed = T)$fit + 
                                                                                          1.96 * predict(flowering_anthropogenic, 
                                                                                                         newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                         se.fit = T,
                                                                                                         exclude = c("s(time)", "s(site)", "s(size,site)", "s(prevwinterR,site)", "s(abLarge,site)"), 
                                                                                                         newdata.guaranteed = T)$se.fit)
flowering_pred_noRE$lwr[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(flowering_anthropogenic, 
                                                                                                newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                se.fit = T,
                                                                                                exclude = c("s(time)", "s(site)", "s(size,site)", "s(prevwinterR,site)", "s(abLarge,site)"), 
                                                                                                newdata.guaranteed = T)$fit - 
                                                                                          1.96 * predict(flowering_anthropogenic, 
                                                                                                         newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                         se.fit = T,
                                                                                                         exclude = c("s(time)", "s(site)", "s(size,site)", "s(prevwinterR,site)", "s(abLarge,site)"), 
                                                                                                         newdata.guaranteed = T)$se.fit)


flowering_pred_noRE$prevfallR_unscaled[which(flowering_pred_noRE$dormancy == "Natural")] = flowering_pred_noRE$prevfallR[which(flowering_pred_noRE$dormancy == "Natural")] * (2 * sd(prevfallR_timeSeries_natural$prevfallR_unscaled)) + mean(prevfallR_timeSeries_natural$prevfallR_unscaled)
flowering_pred_noRE$prevwinterR_unscaled[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = flowering_pred_noRE$prevwinterR[which(flowering_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(prevwinterR_timeSeries_anthropogenic$prevwinterR_unscaled)) + mean(prevwinterR_timeSeries_anthropogenic$prevwinterR_unscaled)

flowering_pred_noRE$abLarge_unscaled[which(flowering_pred_noRE$dormancy == "Natural")] = flowering_pred_noRE$abLarge[which(flowering_pred_noRE$dormancy == "Natural")] * (2 * sd(yearly_density_per_square_natural$abLarge_unscaled)) + mean(yearly_density_per_square_natural$abLarge_unscaled)
flowering_pred_noRE$abLarge_unscaled[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = flowering_pred_noRE$abLarge[which(flowering_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(yearly_density_per_square_anthropogenic$abLarge_unscaled)) + mean(yearly_density_per_square_anthropogenic$abLarge_unscaled)

flowering_pred_noRE$rainfall[which(flowering_pred_noRE$dormancy == "Natural")] = flowering_pred_noRE$prevfallR[which(flowering_pred_noRE$dormancy == "Natural")]
flowering_pred_noRE$rainfall_unscaled[which(flowering_pred_noRE$dormancy == "Natural")] = flowering_pred_noRE$prevfallR_unscaled[which(flowering_pred_noRE$dormancy == "Natural")]

flowering_pred_noRE$rainfall[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = flowering_pred_noRE$prevwinterR[which(flowering_pred_noRE$dormancy == "Anthropogenic")]
flowering_pred_noRE$rainfall_unscaled[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = flowering_pred_noRE$prevwinterR_unscaled[which(flowering_pred_noRE$dormancy == "Anthropogenic")]


flowering_pred_noRE$abLarge_unscaled = as.factor(flowering_pred_noRE$abLarge_unscaled)


flowering_plot_rain_density_natural = ggplot(flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], aes(rainfall_unscaled, flowering,
                                                                                                                          col = abLarge_unscaled, group = abLarge_unscaled)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = abLarge_unscaled), alpha = 0.1, col = NA) +
  scale_colour_manual(name = "Aboveground density\nof large individuals\n(in ind./m2)",
                      values = palette.colors(3, palette = "Dark 2")) +
  scale_fill_manual(name = "Aboveground density\nof large individuals\n(in ind./m2)",
                    values = palette.colors(3, palette = "Dark 2")) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  xlab("Previous fall cumulative\nrainfall (Sep-Nov) (in mm)") +
  ylab("Flowering\nprobability") +
  ggtitle("Natural") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(colour = colPlot, size = 8, family = font,
                                  margin = margin(t = 0, r = 0, b = 5, l = 0)))


flowering_plot_rain_density_anthropogenic = ggplot(flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], aes(rainfall_unscaled, flowering,
                                                                                                                              col = abLarge_unscaled, group = abLarge_unscaled)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = abLarge_unscaled), alpha = 0.1, col = NA) +
  scale_colour_manual(name = "Aboveground density\nof large individuals\n(in ind./m2)",
                      values = palette.colors(3, palette = "Dark 2")) +
  scale_fill_manual(name = "Aboveground density\nof large individuals\n(in ind./m2)",
                    values = palette.colors(3, palette = "Dark 2")) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  xlab("Previous winter cumulative\nrainfall (Jan-Apr) (in mm)") +
  ylab("Flowering probability") +
  ggtitle("Anthropogenic") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(colour = colPlot, size = 8, family = font,
                                  margin = margin(t = 0, r = 0, b = 5, l = 0))) +
  plot_layout(tag_level = 'new')


## 6.3. Growth as a function of rainfall per density (other covariates at average value) ----
# --------------------------------------------------------------------------------------

growth_pred_noRE = rbind(expand.grid(size = 0,
                                     fallR = 0,
                                     summerT = (summerT_values_anthropogenic - mean(summerT_timeSeries_anthropogenic$summerT_unscaled)) / (2 * sd(summerT_timeSeries_anthropogenic$summerT_unscaled)),
                                     abLarge = (c(2, 6, 10) - mean(yearly_density_per_square_anthropogenic$abLarge_unscaled)) / (2 * sd(yearly_density_per_square_anthropogenic$abLarge_unscaled)),
                                     TSFcont = 0,
                                     dormancy = "Anthropogenic"),
                         expand.grid(size = 0,
                                     fallR = (fallR_values_natural - mean(fallR_timeSeries_natural$fallR_unscaled)) / (2 * sd(fallR_timeSeries_natural$fallR_unscaled)),
                                     summerT = 0,
                                     abLarge = (c(2, 6, 10) - mean(yearly_density_per_square_natural$abLarge_unscaled)) / (2 * sd(yearly_density_per_square_natural$abLarge_unscaled)),
                                     TSFcont = 0,
                                     dormancy = "Natural"))


growth_pred_noRE$sizeNext = NA
growth_pred_noRE$upr = NA
growth_pred_noRE$lwr = NA


growth_pred_noRE$sizeNext[which(growth_pred_noRE$dormancy == "Natural")] = (predict(growth_natural, 
                                                                                                    newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Natural"), ], 
                                                                                                    se.fit = T,
                                                                                                    exclude = c("s(time)", "s(site)", "s(size,site)"), 
                                                                                                    newdata.guaranteed = T)$fit)
growth_pred_noRE$upr[which(growth_pred_noRE$dormancy == "Natural")] = (predict(growth_natural, 
                                                                                           newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Natural"), ], 
                                                                                           se.fit = T,
                                                                                           exclude = c("s(time)", "s(site)", "s(size,site)"), 
                                                                                           newdata.guaranteed = T)$fit + 
                                                                                     1.96 * predict(growth_natural, 
                                                                                                    newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Natural"), ], 
                                                                                                    se.fit = T,
                                                                                                    exclude = c("s(time)", "s(site)", "s(size,site)"), 
                                                                                                    newdata.guaranteed = T)$se.fit)
growth_pred_noRE$lwr[which(growth_pred_noRE$dormancy == "Natural")] = (predict(growth_natural, 
                                                                                           newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Natural"), ], 
                                                                                           se.fit = T,
                                                                                           exclude = c("s(time)", "s(site)", "s(size,site)"), 
                                                                                           newdata.guaranteed = T)$fit - 
                                                                                     1.96 * predict(growth_natural, 
                                                                                                    newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Natural"), ], 
                                                                                                    se.fit = T,
                                                                                                    exclude = c("s(time)", "s(site)", "s(size,site)"), 
                                                                                                    newdata.guaranteed = T)$se.fit)

growth_pred_noRE$sizeNext[which(growth_pred_noRE$dormancy == "Anthropogenic")] = (predict(growth_anthropogenic, 
                                                                                                      newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                      se.fit = T,
                                                                                                      exclude = c("s(time)", "s(site)", "s(size,site)", "s(summerT,site)"), 
                                                                                                      newdata.guaranteed = T)$fit)
growth_pred_noRE$upr[which(growth_pred_noRE$dormancy == "Anthropogenic")] = (predict(growth_anthropogenic, 
                                                                                             newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)", "s(site)", "s(size,site)", "s(summerT,site)"), 
                                                                                             newdata.guaranteed = T)$fit + 
                                                                                       1.96 * predict(growth_anthropogenic, 
                                                                                                      newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                      se.fit = T,
                                                                                                      exclude = c("s(time)", "s(site)", "s(size,site)", "s(summerT,site)"), 
                                                                                                      newdata.guaranteed = T)$se.fit)
growth_pred_noRE$lwr[which(growth_pred_noRE$dormancy == "Anthropogenic")] = (predict(growth_anthropogenic, 
                                                                                             newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)", "s(site)", "s(size,site)", "s(summerT,site)"), 
                                                                                             newdata.guaranteed = T)$fit - 
                                                                                       1.96 * predict(growth_anthropogenic, 
                                                                                                      newdata = growth_pred_noRE[which(growth_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                      se.fit = T,
                                                                                                      exclude = c("s(time)", "s(site)", "s(size,site)", "s(summerT,site)"), 
                                                                                                      newdata.guaranteed = T)$se.fit)


growth_pred_noRE$fallR_unscaled[which(growth_pred_noRE$dormancy == "Natural")] = growth_pred_noRE$fallR[which(growth_pred_noRE$dormancy == "Natural")] * (2 * sd(fallR_timeSeries_natural$fallR_unscaled)) + mean(fallR_timeSeries_natural$fallR_unscaled)
growth_pred_noRE$summerT_unscaled[which(growth_pred_noRE$dormancy == "Anthropogenic")] = growth_pred_noRE$summerT[which(growth_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(summerT_timeSeries_anthropogenic$summerT_unscaled)) + mean(summerT_timeSeries_anthropogenic$summerT_unscaled)

growth_pred_noRE$abLarge_unscaled[which(growth_pred_noRE$dormancy == "Natural")] = growth_pred_noRE$abLarge[which(growth_pred_noRE$dormancy == "Natural")] * (2 * sd(yearly_density_per_square_natural$abLarge_unscaled)) + mean(yearly_density_per_square_natural$abLarge_unscaled)
growth_pred_noRE$abLarge_unscaled[which(growth_pred_noRE$dormancy == "Anthropogenic")] = growth_pred_noRE$abLarge[which(growth_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(yearly_density_per_square_anthropogenic$abLarge_unscaled)) + mean(yearly_density_per_square_anthropogenic$abLarge_unscaled)

growth_pred_noRE$abLarge_unscaled = as.factor(growth_pred_noRE$abLarge_unscaled)
growth_pred_noRE$dormancy = factor(growth_pred_noRE$dormancy,
                                   levels = c("Natural", "Anthropogenic"))


growth_plot_rain_density_natural = ggplot(growth_pred_noRE[which(growth_pred_noRE$dormancy == "Natural"), ], aes(fallR_unscaled, sizeNext,
                                                                    col = abLarge_unscaled, group = abLarge_unscaled)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = abLarge_unscaled), alpha = 0.1, col = NA) +
  scale_colour_manual(name = "Aboveground density\nof large individuals\n(in ind./m2)",
                      values = palette.colors(3, palette = "Dark 2")) +
  scale_fill_manual(name = "Aboveground density\nof large individuals\n(in ind./m2)",
                    values = palette.colors(3, palette = "Dark 2")) +
  xlab("Next fall cumulative\nrainfall (Jan-Apr) (in ºC)") +
  ylab("Size next year") +
  ggtitle("Natural") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(colour = colPlot, size = 8, family = font,
                                  margin = margin(t = 0, r = 0, b = 5, l = 0)))

growth_plot_temp_density_anthropogenic = ggplot(growth_pred_noRE[which(growth_pred_noRE$dormancy == "Anthropogenic"), ], aes(summerT_unscaled, sizeNext,
                                                                                                                 col = abLarge_unscaled, group = abLarge_unscaled)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = abLarge_unscaled), alpha = 0.1, col = NA) +
  scale_colour_manual(name = "Aboveground density\nof large individuals\n(in ind./m2)",
                      values = palette.colors(3, palette = "Dark 2")) +
  scale_fill_manual(name = "Aboveground density\nof large individuals\n(in ind./m2)",
                    values = palette.colors(3, palette = "Dark 2")) +
  xlab("Next summer mean max. daily\ntemperature (Jan-Apr) (in ºC)") +
  ylab("Size next year") +
  ggtitle("Anthropogenic") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(colour = colPlot, size = 8, family = font,
                                  margin = margin(t = 0, r = 0, b = 5, l = 0))) +
  plot_layout(tag_level = 'new')


## 6.4. Final plot ----
# ----------------

flowering_plot_rain_density = flowering_plot_rain_density_natural + flowering_plot_rain_density_anthropogenic 

growth_plot_climate_density_natural = growth_plot_rain_density_natural + growth_plot_temp_density_anthropogenic 


png(filename = "Output/Plots/AppendixS1_Figure3.png", 
    width = 13, 
    height = 15, 
    units = "cm", 
    bg = "white", 
    res = 600, 
    type = "cairo")

seedlingSize_plot_rain_density + flowering_plot_rain_density + 
  growth_plot_climate_density_natural +
  plot_layout(guides = "collect", 
              ncol = 1) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") &
  theme(plot.tag = element_text(size = 10))
  

dev.off() 





###########################################################################
#
# 7. Plotting predictions - Appendix - TSF and size effects ----
#
###########################################################################

## 7.1. Survival as a function of TSF (other covariates at average value) ----
# -----------------------------------------------------------------------

surv_pred_noRE = expand.grid(size = 0,
                             fallR = 0,
                             summerT = 0,
                             abLarge = 0,
                             TSFcont = (TSFcont_values_natural - mean(droso_natural$TSFcont_unscaled, na.rm = T)) / (2 * sd(droso_natural$TSFcont_unscaled, na.rm = T)),
                             dormancy = "Natural")

surv_pred_noRE$surv = NA
surv_pred_noRE$upr = NA
surv_pred_noRE$lwr = NA


surv_pred_noRE$surv[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                     newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                     se.fit = T,
                                                                                     exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                     newdata.guaranteed = T)$fit)
surv_pred_noRE$upr[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                    newdata.guaranteed = T)$fit + 
                                                                              1.96 * predict(surv_natural, 
                                                                                             newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                             newdata.guaranteed = T)$se.fit)
surv_pred_noRE$lwr[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                    newdata.guaranteed = T)$fit - 
                                                                              1.96 * predict(surv_natural, 
                                                                                             newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                             newdata.guaranteed = T)$se.fit)

surv_pred_noRE$TSFcont_unscaled[which(surv_pred_noRE$dormancy == "Natural")] = surv_pred_noRE$TSFcont[which(surv_pred_noRE$dormancy == "Natural")] * (2 * sd(droso_natural$TSFcont_unscaled)) + mean(droso_natural$TSFcont_unscaled)

surv_plot_TSF = ggplot(surv_pred_noRE, aes(TSFcont_unscaled, surv)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, col = NA) +
  xlab("Time since fire (in years)") +
  ylab("Survival probability") +
  ggtitle("Natural") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(colour = colPlot, size = 8, family = font,
                                  margin = margin(t = 0, r = 0, b = 5, l = 0)))


## 7.2. Survival as a function of size (other covariates at average value) ----
# ------------------------------------------------------------------------

surv_pred_noRE = rbind(expand.grid(size = (seq(min(droso_anthropogenic$size_unscaled[which(!is.na(droso_anthropogenic$surv))]), max(droso_anthropogenic$size_unscaled[which(!is.na(droso_anthropogenic$surv))]), length.out = 50) - mean(droso_anthropogenic$size_unscaled, na.rm = T)) / (2 * sd(droso_anthropogenic$size_unscaled, na.rm = T)),
                                   fallR = 0,
                                   summerT = 0,
                                   abLarge = 0,
                                   TSFcont = 0,
                                   dormancy = "Anthropogenic"),
                            expand.grid(size = (seq(min(droso_natural$size_unscaled[which(!is.na(droso_natural$surv))]), max(droso_natural$size_unscaled[which(!is.na(droso_natural$surv))]), length.out = 50) - mean(droso_natural$size_unscaled, na.rm = T)) / (2 * sd(droso_natural$size_unscaled, na.rm = T)),
                                        fallR = 0,
                                        summerT = 0,
                                        abLarge = 0,
                                        TSFcont = 0,
                                        dormancy = "Natural"))

surv_pred_noRE$surv = NA
surv_pred_noRE$upr = NA
surv_pred_noRE$lwr = NA


surv_pred_noRE$surv[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                                    se.fit = T,
                                                                                                    exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                                    newdata.guaranteed = T)$fit)
surv_pred_noRE$upr[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                              newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                              se.fit = T,
                                                                                              exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                              newdata.guaranteed = T)$fit + 
                                                                                        1.96 * predict(surv_natural, 
                                                                                                       newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                                       se.fit = T,
                                                                                                       exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                                       newdata.guaranteed = T)$se.fit)
surv_pred_noRE$lwr[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                              newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                              se.fit = T,
                                                                                              exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                              newdata.guaranteed = T)$fit - 
                                                                                        1.96 * predict(surv_natural, 
                                                                                                       newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                                       se.fit = T,
                                                                                                       exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                                       newdata.guaranteed = T)$se.fit)


surv_pred_noRE$surv[which(surv_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                                          newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                          se.fit = T,
                                                                                                          exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                                          newdata.guaranteed = T)$fit)
surv_pred_noRE$upr[which(surv_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                    se.fit = T,
                                                                                                    exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                                    newdata.guaranteed = T)$fit + 
                                                                                              1.96 * predict(surv_anthropogenic, 
                                                                                                             newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                             se.fit = T,
                                                                                                             exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                                             newdata.guaranteed = T)$se.fit)
surv_pred_noRE$lwr[which(surv_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                    se.fit = T,
                                                                                                    exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                                    newdata.guaranteed = T)$fit - 
                                                                                              1.96 * predict(surv_anthropogenic, 
                                                                                                             newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                             se.fit = T,
                                                                                                             exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                                             newdata.guaranteed = T)$se.fit)

surv_pred_noRE$size_unscaled[which(surv_pred_noRE$dormancy == "Natural")] = surv_pred_noRE$size[which(flowering_pred_noRE$dormancy == "Natural")] * (2 * sd(droso_natural$size_unscaled)) + mean(droso_natural$size_unscaled)

surv_pred_noRE$size_unscaled[which(surv_pred_noRE$dormancy == "Anthropogenic")] = surv_pred_noRE$size[which(surv_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(droso_anthropogenic$size_unscaled)) + mean(droso_anthropogenic$size_unscaled)

surv_pred_noRE$dormancy = factor(surv_pred_noRE$dormancy, 
                                      levels = c("Natural", "Anthropogenic"))

surv_plot_size = ggplot(surv_pred_noRE, aes(size_unscaled, surv)) +
  facet_wrap(~dormancy, scales = "free_x") +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, col = NA) +
  xlab("Size") +
  ylab("Survival probability") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(colour = colPlot, size = 8, family = font,
                                  margin = margin(t = 0, r = 0, b = 5, l = 0)))



## 7.3. Flowering as a function of TSF (other covariates at average value) ----
# ------------------------------------------------------------------------

flowering_pred_noRE = expand.grid(size = 0,
                                  prevfallR = 0,
                                  prevwinterT = 0,
                                  abLarge = 0,
                                  TSFcont = (TSFcont_values_natural - mean(droso_natural$TSFcont_unscaled, na.rm = T)) / (2 * sd(droso_natural$TSFcont_unscaled, na.rm = T)),
                                  dormancy = "Natural")

flowering_pred_noRE$flowering = NA
flowering_pred_noRE$upr = NA
flowering_pred_noRE$lwr = NA


flowering_pred_noRE$flowering[which(flowering_pred_noRE$dormancy == "Natural")] = plogis(predict(flowering_natural, 
                                                                                                    newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                                    se.fit = T,
                                                                                                    exclude = c("s(time)", "s(site)"), 
                                                                                                    newdata.guaranteed = T)$fit)
flowering_pred_noRE$upr[which(flowering_pred_noRE$dormancy == "Natural")] = plogis(predict(flowering_natural, 
                                                                                              newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                              se.fit = T,
                                                                                              exclude = c("s(time)", "s(site)"), 
                                                                                              newdata.guaranteed = T)$fit + 
                                                                                        1.96 * predict(flowering_natural, 
                                                                                                       newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                                       se.fit = T,
                                                                                                       exclude = c("s(time)", "s(site)"), 
                                                                                                       newdata.guaranteed = T)$se.fit)
flowering_pred_noRE$lwr[which(flowering_pred_noRE$dormancy == "Natural")] = plogis(predict(flowering_natural, 
                                                                                              newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                              se.fit = T,
                                                                                              exclude = c("s(time)", "s(site)"), 
                                                                                              newdata.guaranteed = T)$fit - 
                                                                                        1.96 * predict(flowering_natural, 
                                                                                                       newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                                       se.fit = T,
                                                                                                       exclude = c("s(time)", "s(site)"), 
                                                                                                       newdata.guaranteed = T)$se.fit)

flowering_pred_noRE$TSFcont_unscaled[which(flowering_pred_noRE$dormancy == "Natural")] = flowering_pred_noRE$TSFcont[which(flowering_pred_noRE$dormancy == "Natural")] * (2 * sd(droso_natural$TSFcont_unscaled)) + mean(droso_natural$TSFcont_unscaled)

flowering_plot_TSF = ggplot(flowering_pred_noRE, aes(TSFcont_unscaled, flowering)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, col = NA) +
  xlab("Time since fire (in years)") +
  ylab("Flowering probability") +
  ggtitle("Natural") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(colour = colPlot, size = 8, family = font,
                                  margin = margin(t = 0, r = 0, b = 5, l = 0)))




## 7.4. Number of flowers as a function of TSF (other covariates at average value) ----
# --------------------------------------------------------------------------------

nbFlowers_pred_noRE = expand.grid(size = (mean_size_repro_natural - mean(droso_natural$size_unscaled, na.rm = T)) / (2 * sd(droso_natural$size_unscaled, na.rm = T)),
                                        prevfallR = 0,
                                        prevwinterR = 0,
                                        prevwinterT = 0,
                                        abLarge = 0,
                                        TSFcont = (TSFcont_values_natural - mean(droso_natural$TSFcont_unscaled, na.rm = T)) / (2 * sd(droso_natural$TSFcont_unscaled, na.rm = T)),
                                        dormancy = "Natural")

nbFlowers_pred_noRE$nbFlowers = NA
nbFlowers_pred_noRE$upr = NA
nbFlowers_pred_noRE$lwr = NA


nbFlowers_pred_noRE$nbFlowers[which(nbFlowers_pred_noRE$dormancy == "Natural")] = exp(predict(nbFlow_natural, 
                                                                                              newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Natural"), ], 
                                                                                              se.fit = T,
                                                                                              exclude = c("s(time)", "s(site)"), 
                                                                                              newdata.guaranteed = T)$fit)
nbFlowers_pred_noRE$upr[which(nbFlowers_pred_noRE$dormancy == "Natural")] = exp(predict(nbFlow_natural, 
                                                                                        newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Natural"), ], 
                                                                                        se.fit = T,
                                                                                        exclude = c("s(time)", "s(site)"), 
                                                                                        newdata.guaranteed = T)$fit + 
                                                                                  1.96 * predict(nbFlow_natural, 
                                                                                                 newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Natural"), ], 
                                                                                                 se.fit = T,
                                                                                                 exclude = c("s(time)", "s(site)"), 
                                                                                                 newdata.guaranteed = T)$se.fit)
nbFlowers_pred_noRE$lwr[which(nbFlowers_pred_noRE$dormancy == "Natural")] = exp(predict(nbFlow_natural, 
                                                                                        newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Natural"), ], 
                                                                                        se.fit = T,
                                                                                        exclude = c("s(time)", "s(site)"), 
                                                                                        newdata.guaranteed = T)$fit - 
                                                                                  1.96 * predict(nbFlow_natural, 
                                                                                                 newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Natural"), ], 
                                                                                                 se.fit = T,
                                                                                                 exclude = c("s(time)", "s(site)"), 
                                                                                                 newdata.guaranteed = T)$se.fit)

nbFlowers_pred_noRE$TSFcont_unscaled[which(nbFlowers_pred_noRE$dormancy == "Natural")] = nbFlowers_pred_noRE$TSFcont[which(nbFlowers_pred_noRE$dormancy == "Natural")] * (2 * sd(droso_natural$TSFcont_unscaled)) + mean(droso_natural$TSFcont_unscaled)

nbFlowers_plot_TSF = ggplot(nbFlowers_pred_noRE, aes(TSFcont_unscaled, nbFlowers)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, col = NA) +
  xlab("Time since fire (in years)") +
  ylab("Number of flowers") +
  ggtitle("Natural") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(colour = colPlot, size = 8, family = font,
                                  margin = margin(t = 0, r = 0, b = 5, l = 0)))




## 7.5. Flowering as a function of size (other covariates at average value) ----
# -------------------------------------------------------------------------

flowering_pred_noRE = rbind(expand.grid(size = (seq(min(droso_anthropogenic$size_unscaled[which(!is.na(droso_anthropogenic$fl))]), max(droso_anthropogenic$size_unscaled[which(!is.na(droso_anthropogenic$fl))]), length.out = 50) - mean(droso_anthropogenic$size_unscaled, na.rm = T)) / (2 * sd(droso_anthropogenic$size_unscaled, na.rm = T)),
                                        prevfallR = 0,
                                        prevwinterR = 0,
                                        prevwinterT = 0,
                                        abLarge = 0,
                                        TSFcont = 0,
                                        dormancy = "Anthropogenic"),
                            expand.grid(size = (seq(min(droso_natural$size_unscaled[which(!is.na(droso_natural$fl))]), max(droso_natural$size_unscaled[which(!is.na(droso_natural$fl))]), length.out = 50) - mean(droso_natural$size_unscaled, na.rm = T)) / (2 * sd(droso_natural$size_unscaled, na.rm = T)),
                                        prevfallR = 0,
                                        prevwinterR = 0,
                                        prevwinterT = 0,
                                        abLarge = 0,
                                        TSFcont = 0,
                                        dormancy = "Natural"))

flowering_pred_noRE$flowering = NA
flowering_pred_noRE$upr = NA
flowering_pred_noRE$lwr = NA


flowering_pred_noRE$flowering[which(flowering_pred_noRE$dormancy == "Natural")] = plogis(predict(flowering_natural, 
                                                                                                    newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                                    se.fit = T,
                                                                                                    exclude = c("s(time)", "s(site)"), 
                                                                                                    newdata.guaranteed = T)$fit)
flowering_pred_noRE$upr[which(flowering_pred_noRE$dormancy == "Natural")] = plogis(predict(flowering_natural, 
                                                                                              newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                              se.fit = T,
                                                                                              exclude = c("s(time)", "s(site)"), 
                                                                                              newdata.guaranteed = T)$fit + 
                                                                                        1.96 * predict(flowering_natural, 
                                                                                                       newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                                       se.fit = T,
                                                                                                       exclude = c("s(time)", "s(site)"), 
                                                                                                       newdata.guaranteed = T)$se.fit)
flowering_pred_noRE$lwr[which(flowering_pred_noRE$dormancy == "Natural")] = plogis(predict(flowering_natural, 
                                                                                              newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                              se.fit = T,
                                                                                              exclude = c("s(time)", "s(site)"), 
                                                                                              newdata.guaranteed = T)$fit - 
                                                                                        1.96 * predict(flowering_natural, 
                                                                                                       newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Natural"), ], 
                                                                                                       se.fit = T,
                                                                                                       exclude = c("s(time)", "s(site)"), 
                                                                                                       newdata.guaranteed = T)$se.fit)


flowering_pred_noRE$flowering[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(flowering_anthropogenic, 
                                                                                                      newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                      se.fit = T,
                                                                                                      exclude = c("s(time)", "s(site)", "s(size,site)", "s(prevwinterR,site)", "s(abLarge,site)"), 
                                                                                                      newdata.guaranteed = T)$fit)
flowering_pred_noRE$upr[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(flowering_anthropogenic, 
                                                                                                newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                se.fit = T,
                                                                                                exclude = c("s(time)", "s(site)", "s(size,site)", "s(prevwinterR,site)", "s(abLarge,site)"), 
                                                                                                newdata.guaranteed = T)$fit + 
                                                                                          1.96 * predict(flowering_anthropogenic, 
                                                                                                         newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                         se.fit = T,
                                                                                                         exclude = c("s(time)", "s(site)", "s(size,site)", "s(prevwinterR,site)", "s(abLarge,site)"), 
                                                                                                         newdata.guaranteed = T)$se.fit)
flowering_pred_noRE$lwr[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(flowering_anthropogenic, 
                                                                                                newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                se.fit = T,
                                                                                                exclude = c("s(time)", "s(site)", "s(size,site)", "s(prevwinterR,site)", "s(abLarge,site)"), 
                                                                                                newdata.guaranteed = T)$fit - 
                                                                                          1.96 * predict(flowering_anthropogenic, 
                                                                                                         newdata = flowering_pred_noRE[which(flowering_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                         se.fit = T,
                                                                                                         exclude = c("s(time)", "s(site)", "s(size,site)", "s(prevwinterR,site)", "s(abLarge,site)"), 
                                                                                                         newdata.guaranteed = T)$se.fit)

flowering_pred_noRE$size_unscaled[which(flowering_pred_noRE$dormancy == "Natural")] = flowering_pred_noRE$size[which(flowering_pred_noRE$dormancy == "Natural")] * (2 * sd(droso_natural$size_unscaled)) + mean(droso_natural$size_unscaled)

flowering_pred_noRE$size_unscaled[which(flowering_pred_noRE$dormancy == "Anthropogenic")] = flowering_pred_noRE$size[which(flowering_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(droso_anthropogenic$size_unscaled)) + mean(droso_anthropogenic$size_unscaled)

flowering_pred_noRE$dormancy = factor(flowering_pred_noRE$dormancy, 
                                      levels = c("Natural", "Anthropogenic"))

flowering_plot_size = ggplot(flowering_pred_noRE, aes(size_unscaled, flowering)) +
  facet_wrap(~dormancy, scales = "free_x") +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, col = NA) +
  xlab("Size") +
  ylab("Flowering probability") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(colour = colPlot, size = 8, family = font,
                                  margin = margin(t = 0, r = 0, b = 5, l = 0)))


## 7.6. Number of flowers as a function of size (other covariates at average value) ----
# ---------------------------------------------------------------------------------

nbFlowers_pred_noRE = rbind(expand.grid(size = (seq(min(droso_anthropogenic$size_unscaled[which(!is.na(droso_anthropogenic$nbFlow))]), max(droso_anthropogenic$size_unscaled[which(!is.na(droso_anthropogenic$nbFlow))]), length.out = 50) - mean(droso_anthropogenic$size_unscaled, na.rm = T)) / (2 * sd(droso_anthropogenic$size_unscaled, na.rm = T)),
                                        prevfallR = 0,
                                        prevwinterR = 0,
                                        prevwinterT = 0,
                                        abLarge = 0,
                                        TSFcont = 0,
                                        dormancy = "Anthropogenic"),
                            expand.grid(size = (seq(min(droso_natural$size_unscaled[which(!is.na(droso_natural$nbFlow))]), max(droso_natural$size_unscaled[which(!is.na(droso_natural$nbFlow))]), length.out = 50) - mean(droso_natural$size_unscaled, na.rm = T)) / (2 * sd(droso_natural$size_unscaled, na.rm = T)),
                                        prevfallR = 0,
                                        prevwinterR = 0,
                                        prevwinterT = 0,
                                        abLarge = 0,
                                        TSFcont = 0,
                                        dormancy = "Natural"))

nbFlowers_pred_noRE$nbFlowers = NA
nbFlowers_pred_noRE$upr = NA
nbFlowers_pred_noRE$lwr = NA


nbFlowers_pred_noRE$nbFlowers[which(nbFlowers_pred_noRE$dormancy == "Natural")] = exp(predict(nbFlow_natural, 
                                                                                                    newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Natural"), ], 
                                                                                                    se.fit = T,
                                                                                                    exclude = c("s(time)", "s(site)"), 
                                                                                                    newdata.guaranteed = T)$fit)
nbFlowers_pred_noRE$upr[which(nbFlowers_pred_noRE$dormancy == "Natural")] = exp(predict(nbFlow_natural, 
                                                                                              newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Natural"), ], 
                                                                                              se.fit = T,
                                                                                              exclude = c("s(time)", "s(site)"), 
                                                                                              newdata.guaranteed = T)$fit + 
                                                                                        1.96 * predict(nbFlow_natural, 
                                                                                                       newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Natural"), ], 
                                                                                                       se.fit = T,
                                                                                                       exclude = c("s(time)", "s(site)"), 
                                                                                                       newdata.guaranteed = T)$se.fit)
nbFlowers_pred_noRE$lwr[which(nbFlowers_pred_noRE$dormancy == "Natural")] = exp(predict(nbFlow_natural, 
                                                                                              newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Natural"), ], 
                                                                                              se.fit = T,
                                                                                              exclude = c("s(time)", "s(site)"), 
                                                                                              newdata.guaranteed = T)$fit - 
                                                                                        1.96 * predict(nbFlow_natural, 
                                                                                                       newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Natural"), ], 
                                                                                                       se.fit = T,
                                                                                                       exclude = c("s(time)", "s(site)"), 
                                                                                                       newdata.guaranteed = T)$se.fit)


nbFlowers_pred_noRE$nbFlowers[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic")] = exp(predict(nbFlow_anthropogenic, 
                                                                                                      newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                      se.fit = T,
                                                                                                      exclude = c("s(time)", "s(site)", "s(prevfallR,site)", "s(size,site)"), 
                                                                                                      newdata.guaranteed = T)$fit)
nbFlowers_pred_noRE$upr[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic")] = exp(predict(nbFlow_anthropogenic, 
                                                                                                newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                se.fit = T,
                                                                                                exclude = c("s(time)", "s(site)", "s(prevfallR,site)", "s(size,site)"), 
                                                                                                newdata.guaranteed = T)$fit + 
                                                                                          1.96 * predict(nbFlow_anthropogenic, 
                                                                                                         newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                         se.fit = T,
                                                                                                         exclude = c("s(time)", "s(site)", "s(prevfallR,site)", "s(size,site)"), 
                                                                                                         newdata.guaranteed = T)$se.fit)
nbFlowers_pred_noRE$lwr[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic")] = exp(predict(nbFlow_anthropogenic, 
                                                                                                newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                se.fit = T,
                                                                                                exclude = c("s(time)", "s(site)", "s(prevfallR,site)", "s(size,site)"), 
                                                                                                newdata.guaranteed = T)$fit - 
                                                                                          1.96 * predict(nbFlow_anthropogenic, 
                                                                                                         newdata = nbFlowers_pred_noRE[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                                         se.fit = T,
                                                                                                         exclude = c("s(time)", "s(site)", "s(prevfallR,site)", "s(size,site)"), 
                                                                                                         newdata.guaranteed = T)$se.fit)

nbFlowers_pred_noRE$size_unscaled[which(nbFlowers_pred_noRE$dormancy == "Natural")] = nbFlowers_pred_noRE$size[which(nbFlowers_pred_noRE$dormancy == "Natural")] * (2 * sd(droso_natural$size_unscaled)) + mean(droso_natural$size_unscaled)

nbFlowers_pred_noRE$size_unscaled[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic")] = nbFlowers_pred_noRE$size[which(nbFlowers_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(droso_anthropogenic$size_unscaled)) + mean(droso_anthropogenic$size_unscaled)

nbFlowers_pred_noRE$dormancy = factor(nbFlowers_pred_noRE$dormancy, 
                                      levels = c("Natural", "Anthropogenic"))

nbFlowers_plot_size = ggplot(nbFlowers_pred_noRE, aes(size_unscaled, nbFlowers)) +
  facet_wrap(~dormancy, scales = "free_x") +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, col = NA) +
  xlab("Size") +
  ylab("Number of flowers") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(colour = colPlot, size = 8, family = font,
                                  margin = margin(t = 0, r = 0, b = 5, l = 0)))



## 7.7. Survival as a function of rainfall per size (other covariates at average value) ----
# --------------------------------------------------------------------------------------------

surv_pred_noRE = rbind(expand.grid(size = (size_values_anthropogenic[seq(1, length(size_values_anthropogenic), length.out = 5)] - mean(droso_anthropogenic$size_unscaled)) / (2 * sd(droso_anthropogenic$size_unscaled)),
                                   fallR = (fallR_values_anthropogenic - mean(fallR_timeSeries_anthropogenic$fallR_unscaled)) / (2 * sd(fallR_timeSeries_anthropogenic$fallR_unscaled)),
                                   summerT = 0,
                                   abLarge = 0,
                                   TSFcont = 0,
                                   dormancy = "Anthropogenic"),
                       expand.grid(size = (size_values_natural[seq(1, length(size_values_natural), length.out = 5)] - mean(droso_natural$size_unscaled)) / (2 * sd(droso_natural$size_unscaled)),
                                   fallR = (fallR_values_natural - mean(fallR_timeSeries_natural$fallR_unscaled)) / (2 * sd(fallR_timeSeries_natural$fallR_unscaled)),
                                   summerT = 0,
                                   abLarge = 0,
                                   TSFcont = 0,
                                   dormancy = "Natural"))


surv_pred_noRE$surv = NA
surv_pred_noRE$upr = NA
surv_pred_noRE$lwr = NA


surv_pred_noRE$surv[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                     newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                     se.fit = T,
                                                                                     exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                     newdata.guaranteed = T)$fit)
surv_pred_noRE$upr[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                    newdata.guaranteed = T)$fit + 
                                                                              1.96 * predict(surv_natural, 
                                                                                             newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                             newdata.guaranteed = T)$se.fit)
surv_pred_noRE$lwr[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                    newdata.guaranteed = T)$fit - 
                                                                              1.96 * predict(surv_natural, 
                                                                                             newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                             newdata.guaranteed = T)$se.fit)


surv_pred_noRE$surv[which(surv_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                       newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                       se.fit = T,
                                                                                       exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                       newdata.guaranteed = T)$fit)
surv_pred_noRE$upr[which(surv_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                      newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                      se.fit = T,
                                                                                      exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                      newdata.guaranteed = T)$fit + 
                                                                                1.96 * predict(surv_anthropogenic, 
                                                                                               newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                               se.fit = T,
                                                                                               exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                               newdata.guaranteed = T)$se.fit)
surv_pred_noRE$lwr[which(surv_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                      newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                      se.fit = T,
                                                                                      exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                      newdata.guaranteed = T)$fit - 
                                                                                1.96 * predict(surv_anthropogenic, 
                                                                                               newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                               se.fit = T,
                                                                                               exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                               newdata.guaranteed = T)$se.fit)

surv_pred_noRE$size_unscaled[which(surv_pred_noRE$dormancy == "Natural")] = surv_pred_noRE$size[which(surv_pred_noRE$dormancy == "Natural")] * (2 * sd(droso_natural$size_unscaled)) + mean(droso_natural$size_unscaled)
surv_pred_noRE$size_unscaled[which(surv_pred_noRE$dormancy == "Anthropogenic")] = surv_pred_noRE$size[which(surv_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(droso_anthropogenic$size_unscaled)) + mean(droso_anthropogenic$size_unscaled)

surv_pred_noRE$fallR_unscaled[which(surv_pred_noRE$dormancy == "Natural")] = surv_pred_noRE$fallR[which(surv_pred_noRE$dormancy == "Natural")] * (2 * sd(fallR_timeSeries_natural$fallR_unscaled)) + mean(fallR_timeSeries_natural$fallR_unscaled)
surv_pred_noRE$fallR_unscaled[which(surv_pred_noRE$dormancy == "Anthropogenic")] = surv_pred_noRE$fallR[which(surv_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(fallR_timeSeries_anthropogenic$fallR_unscaled)) + mean(fallR_timeSeries_anthropogenic$fallR_unscaled)

surv_pred_noRE$dormancy = factor(surv_pred_noRE$dormancy,
                                 levels = c("Natural", "Anthropogenic"))

surv_plot_rain_size = ggplot(surv_pred_noRE, aes(fallR_unscaled, surv,
                                                 col = size_unscaled, group = size_unscaled)) +
  facet_wrap(~dormancy) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = size_unscaled), alpha = 0.1, col = NA) +
  scale_colour_viridis(name = "Individual size") +
  scale_fill_viridis(name = "Individual size") +
  xlab("Next fall cumulative\nrainfall (Sep-Nov) (in mm)") +
  ylab("Survival probability") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"))




## 7.8. Survival as a function of temperature per size (other covariates at average value) ----
# --------------------------------------------------------------------------------------------

surv_pred_noRE = rbind(expand.grid(size = (size_values_anthropogenic[seq(1, length(size_values_anthropogenic), length.out = 5)] - mean(droso_anthropogenic$size_unscaled)) / (2 * sd(droso_anthropogenic$size_unscaled)),
                                   fallR = 0,
                                   summerT = (summerT_values_anthropogenic - mean(summerT_timeSeries_anthropogenic$summerT_unscaled)) / (2 * sd(summerT_timeSeries_anthropogenic$summerT_unscaled)),
                                   abLarge = 0,
                                   TSFcont = 0,
                                   dormancy = "Anthropogenic"),
                       expand.grid(size = (size_values_natural[seq(1, length(size_values_natural), length.out = 5)] - mean(droso_natural$size_unscaled)) / (2 * sd(droso_natural$size_unscaled)),
                                   fallR = 0,
                                   summerT = (summerT_values_natural - mean(summerT_timeSeries_natural$summerT_unscaled)) / (2 * sd(summerT_timeSeries_natural$summerT_unscaled)),
                                   abLarge = 0,
                                   TSFcont = 0,
                                   dormancy = "Natural"))


surv_pred_noRE$surv = NA
surv_pred_noRE$upr = NA
surv_pred_noRE$lwr = NA


surv_pred_noRE$surv[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                     newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                     se.fit = T,
                                                                                     exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                     newdata.guaranteed = T)$fit)
surv_pred_noRE$upr[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                    newdata.guaranteed = T)$fit + 
                                                                              1.96 * predict(surv_natural, 
                                                                                             newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                             newdata.guaranteed = T)$se.fit)
surv_pred_noRE$lwr[which(surv_pred_noRE$dormancy == "Natural")] = plogis(predict(surv_natural, 
                                                                                    newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                    se.fit = T,
                                                                                    exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                    newdata.guaranteed = T)$fit - 
                                                                              1.96 * predict(surv_natural, 
                                                                                             newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Natural"), ], 
                                                                                             se.fit = T,
                                                                                             exclude = c("s(time)", "s(site)", "s(summerT,site)"), 
                                                                                             newdata.guaranteed = T)$se.fit)


surv_pred_noRE$surv[which(surv_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                       newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                       se.fit = T,
                                                                                       exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                       newdata.guaranteed = T)$fit)
surv_pred_noRE$upr[which(surv_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                      newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                      se.fit = T,
                                                                                      exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                      newdata.guaranteed = T)$fit + 
                                                                                1.96 * predict(surv_anthropogenic, 
                                                                                               newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                               se.fit = T,
                                                                                               exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                               newdata.guaranteed = T)$se.fit)
surv_pred_noRE$lwr[which(surv_pred_noRE$dormancy == "Anthropogenic")] = plogis(predict(surv_anthropogenic, 
                                                                                      newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                      se.fit = T,
                                                                                      exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                      newdata.guaranteed = T)$fit - 
                                                                                1.96 * predict(surv_anthropogenic, 
                                                                                               newdata = surv_pred_noRE[which(surv_pred_noRE$dormancy == "Anthropogenic"), ], 
                                                                                               se.fit = T,
                                                                                               exclude = c("s(time)", "s(site)", "s(size,site)", "s(fallR,site)"), 
                                                                                               newdata.guaranteed = T)$se.fit)

surv_pred_noRE$size_unscaled[which(surv_pred_noRE$dormancy == "Natural")] = surv_pred_noRE$size[which(surv_pred_noRE$dormancy == "Natural")] * (2 * sd(droso_natural$size_unscaled)) + mean(droso_natural$size_unscaled)
surv_pred_noRE$size_unscaled[which(surv_pred_noRE$dormancy == "Anthropogenic")] = surv_pred_noRE$size[which(surv_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(droso_anthropogenic$size_unscaled)) + mean(droso_anthropogenic$size_unscaled)

surv_pred_noRE$summerT_unscaled[which(surv_pred_noRE$dormancy == "Natural")] = surv_pred_noRE$summerT[which(surv_pred_noRE$dormancy == "Natural")] * (2 * sd(summerT_timeSeries_natural$summerT_unscaled)) + mean(summerT_timeSeries_natural$summerT_unscaled)
surv_pred_noRE$summerT_unscaled[which(surv_pred_noRE$dormancy == "Anthropogenic")] = surv_pred_noRE$summerT[which(surv_pred_noRE$dormancy == "Anthropogenic")] * (2 * sd(summerT_timeSeries_anthropogenic$summerT_unscaled)) + mean(summerT_timeSeries_anthropogenic$summerT_unscaled)

surv_pred_noRE$dormancy = factor(surv_pred_noRE$dormancy,
                                 levels = c("Natural", "Anthropogenic"))

surv_plot_temp_size = ggplot(surv_pred_noRE, aes(summerT_unscaled, surv,
                                                 col = size_unscaled, group = size_unscaled)) +
  facet_wrap(~dormancy) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = size_unscaled), alpha = 0.1, col = NA) +
  scale_colour_viridis(name = "Individual size") +
  scale_fill_viridis(name = "Individual size") +
  xlab("Next summer max. daily mean\ntemperature (Sep-Nov) (in mm)") +
  ylab("Survival probability") +
  theme_general() +
  theme(legend.key.width = unit(0.5, "cm"))


## 7.9. Final plot ----
# ----------------

png(filename = "Output/Plots/AppendixS1_Figure4.png", 
    width = 16, 
    height = 20, 
    units = "cm", 
    bg = "white", 
    res = 600, 
    type = "cairo")

surv_plot_TSF + surv_plot_size + 
  flowering_plot_TSF + nbFlowers_plot_TSF +
  flowering_plot_size + nbFlowers_plot_size +
  surv_plot_rain_size + surv_plot_temp_size +
  plot_layout(guides = "collect",
              ncol = 2) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") &
  theme(plot.tag = element_text(size = 10))

dev.off()
