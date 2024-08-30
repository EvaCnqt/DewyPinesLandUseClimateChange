############################################################################
#
# This script uses the validation data to plot the results.
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


## 1.3. Loading data ----
# ------------------

ibm_results_anthropogenic = read.csv("Output/Projections/ValidationResults_Anthropogenic.csv")
ibm_results_natural = read.csv("Output/Projections/ValidationResults_Natural.csv")

ibm_results = rbind(ibm_results_anthropogenic, ibm_results_natural[, which(colnames(ibm_results_natural) %in% colnames(ibm_results_anthropogenic))])


size_dist_natural = read.csv("Output/Projections/SizeDistribution_Natural.csv")
size_dist_anthropogenic = read.csv("Output/Projections/SizeDistribution_Anthropogenic.csv")

size_dist_df = rbind(size_dist_natural[, which(colnames(size_dist_natural) %in% colnames(size_dist_anthropogenic))], size_dist_anthropogenic)




###########################################################################
#
# 2. Plot results ----
#
###########################################################################

cbbPalette <- c("#111111", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## 2.1. Summarize metrics ----
# -----------------------

mean_log_meanChangeAboveground = aggregate(log_meanChangeAboveground ~ source + dormancy + population + year, 
                                           data = ibm_results[which(!is.infinite(ibm_results$log_meanChangeAboveground)), ], 
                                           FUN = function(x) c(mean(x, na.rm = T), 
                                                               quantile(x, 
                                                                        probs = c(0.025, 
                                                                                  0.975), 
                                                                        na.rm = T)))
mean_log_meanChangeAboveground = cbind(mean_log_meanChangeAboveground[, c("source", "dormancy", "population", "year")],
                                       mean_log_meanChangeAboveground$log_meanChangeAboveground)
colnames(mean_log_meanChangeAboveground) = c("source", "dormancy", "population", "year", 
                                             "mean", "lwr", "upr")


mean_log_lambda = aggregate(log_lambda ~ source + dormancy + population + year, 
                            data = ibm_results[which(!is.infinite(ibm_results$log_lambda)), ], 
                            FUN = function(x) c(mean(x, na.rm = T), 
                                              quantile(x, 
                                                       probs = c(0.025, 
                                                                 0.975), 
                                                       na.rm = T)))
mean_log_lambda = cbind(mean_log_lambda[, c("source", "dormancy", "population", "year")],
                        mean_log_lambda$log_lambda)
colnames(mean_log_lambda) = c("source", "dormancy", "population", "year", 
                            "mean", "lwr", "upr")


mean_pop_size = aggregate(pop_size ~ source + dormancy + population + year, 
                          data = ibm_results, 
                          FUN = function(x) c(mean(x, na.rm = T), 
                                              quantile(x, 
                                                       probs = c(0.025, 
                                                                 0.975), 
                                                       na.rm = T)))
mean_pop_size = cbind(mean_pop_size[, c("source", "dormancy", "population", "year")],
                      mean_pop_size$pop_size)
colnames(mean_pop_size) = c("source", "dormancy", "population", "year", 
                            "mean", "lwr", "upr")


mean_pop_size_repro = aggregate(pop_size_repro ~ source + dormancy + population + year, 
                                data = ibm_results, 
                                FUN = function(x) c(mean(x, na.rm = T),
                                                    quantile(x, 
                                                             probs = c(0.025, 
                                                                       0.975), 
                                                             na.rm = T)))
mean_pop_size_repro = cbind(mean_pop_size_repro[, c("source", "dormancy", "population", "year")],
                            mean_pop_size_repro$pop_size_repro)
colnames(mean_pop_size_repro) = c("source", "dormancy", "population", "year", 
                                  "mean", "lwr", "upr")


## 2.2. Plot results ----
# ------------------

# Log mean change in aboveground population abundance

mean_log_meanChangeAboveground$dormancy = factor(mean_log_meanChangeAboveground$dormancy, 
                                                 levels = c("Natural", "Anthropogenic"))

png(filename = "Output/Plots/Validation_LogMeanChangeAboveground.png", 
    width = 8,
    height = 13,
    units = "cm",
    bg = "white",
    res = 600)

ggplot(mean_log_meanChangeAboveground[which(mean_log_meanChangeAboveground$source == "Projected"), ], 
       aes(x = year, y = mean, colour = dormancy, fill = dormancy)) +
  facet_wrap(~ population, scales = "free_y",
             ncol = 2,
             labeller = labeller(population = c("Bujeo" = "Bujeo",
                                                "MonteraTorero" = "Montera del\nTorero",
                                                "Prisoneros" = "Prisoneros",
                                                "Retin" = "Sierra del\nRetín Disturbed",
                                                "SCarbDist" = "Sierra\nCarbonera\nDisturbed",
                                                "SierraCarboneraY5" = "Sierra\nCarbonera\nYoung",
                                                "SierraRetinY5" = "Sierra del\nRetín Young",
                                                "Vertedero" = "Vertedero"))) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, col = NA) +
  geom_line(linewidth = 0.5) +
  geom_point(data = mean_log_meanChangeAboveground[which(mean_log_meanChangeAboveground$source == "Observed"), ],
             aes(x = year, y = mean),
             size = 0.8) +
  scale_x_continuous(n.breaks = 4) + 
  scale_colour_manual(name = "Habitat type",
                      values = cbbPalette[1:2]) +
  scale_fill_manual(name = "Habitat type",
                    values = cbbPalette[1:2]) +
  labs(x = "Year", 
       y = "log(mean change in aboveground population size)") +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title.x = element_text(size = 8, colour = "black", 
                                    margin = margin(t = 0, r = 0, b = 10, l = 0)), 
        axis.title.y = element_text(size = 8, colour = "black", 
                                    margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 2, r = 0, b = 5, l = 0)), 
        axis.text.y = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        panel.grid = element_blank(),
        strip.background.x = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 8), 
        legend.position = "bottom", 
        legend.key.size = unit(1, "lines"),
        strip.background = element_blank())

dev.off()


# Population size

png(filename = "Output/Plots/Validation_PopulationSize.png", 
    width = 14,
    height = 10,
    units = "cm",
    bg = "white",
    res = 600)

ggplot(mean_pop_size[which(mean_pop_size$source == "Projected"), ], 
       aes(x = year, y = mean, colour = dormancy, fill = dormancy)) +
  facet_wrap(~ population, scales = "free_y",
             ncol = 3,
             labeller = labeller(population = c("Bujeo" = "Bujeo",
                                                "MonteraTorero" = "Montera del\nTorero",
                                                "Prisoneros" = "Prisoneros",
                                                "Retin" = "Sierra del\nRetín Disturbed",
                                                "SCarbDist" = "Sierra\nCarbonera\nDisturbed",
                                                "SierraCarboneraY5" = "Sierra\nCarbonera\nYoung",
                                                "SierraRetinY5" = "Sierra del\nRetín Young",
                                                "Vertedero" = "Vertedero"))) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, col = NA) +
  geom_line(linewidth = 0.5) +
  geom_point(data = mean_pop_size[which(mean_pop_size$source == "Observed"), ],
             aes(x = year, y = mean),
             size = 0.8) +
  scale_x_continuous(breaks = c(2012, 2016, 2020)) + 
  scale_colour_manual(name = "Habitat type",
                      values = cbbPalette[1:2]) +
  scale_fill_manual(name = "Habitat type",
                    values = cbbPalette[1:2]) +
  labs(x = "Year", 
       y = "Aboveground population size") +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title.x = element_text(size = 8, colour = "black", 
                                    margin = margin(t = 0, r = 0, b = 10, l = 0)), 
        axis.title.y = element_text(size = 8, colour = "black", 
                                    margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 2, r = 0, b = 5, l = 0)), 
        axis.text.y = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        panel.grid = element_blank(),
        strip.background.x = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 8), 
        legend.position = "bottom", 
        legend.key.size = unit(1, "lines"),
        strip.background = element_blank())

dev.off()


# Size distribution

png(filename = "Output/Plots/Validation_SizeDistribution.png", 
    width = 16,
    height = 20,
    units = "cm",
    bg = "white",
    res = 600)

ggplot(size_dist_df, 
       aes(x = size_unscaled, colour = srce, fill = srce)) +
  facet_grid(time ~ site, scales = "free",
             labeller = labeller(site = c("Bujeo" = "Bujeo",
                                                    "MonteraTorero" = "Montera del\nTorero",
                                                    "Prisoneros" = "Prisoneros",
                                                    "Retin" = "Sierra del\nRetín\nDisturbed",
                                                    "SCarbDist" = "Sierra\nCarbonera\nDisturbed",
                                                    "SierraCarboneraY5" = "Sierra\nCarbonera\nYoung",
                                                    "SierraRetinY5" = "Sierra del\nRetín Young",
                                                    "Vertedero" = "Vertedero"))) +
  geom_density(alpha = 0.2) +
  scale_colour_manual(name = "Distribution",
                      values = viridis(2, end = 1)) +
  scale_fill_manual(name = "Distribution",
                    values = viridis(2, end = 1)) +
  scale_x_continuous(n.breaks = 4) +
  labs(x = "Size", 
       y = "Density") +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title.x = element_text(size = 8, colour = "black", 
                                    margin = margin(t = , r = 0, b = 10, l = 0)), 
        axis.title.y = element_text(size = 8, colour = "black", 
                                    margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 2, r = 0, b = 5, l = 0)), 
        axis.text.y = element_text(size = 7, colour = "black", 
                                   margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        strip.background.x = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 8), 
        legend.position = "bottom", 
        legend.key.size = unit(1, "lines"),
        panel.grid = element_blank(),
        strip.background = element_blank())

dev.off()