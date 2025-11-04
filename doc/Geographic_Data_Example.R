## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  message   = FALSE,
  warning   = FALSE,
  fig.width = 10,
  fig_height = 10
)

## ----climate_data-est_SCB, eval = FALSE---------------------------------------
# library(patchwork)
# library(nlme)
# library(SCoRES)
# library(ggplot2)
# 
# # Load data
# data(climate_data)
# start_time <- Sys.time()
# 
# # construct confidence region for the increase of the mean temperature (June-August) in North America between the 20th and 21st centuries
# temp = SCB_gls_geospatial(sp_list = climate_data$Z,
#                        level = 2,
#                        data_fit = climate_data$X,
#                        w = c(1,0,0,0),
#                        correlation = climate_data$correlation,
#                        mask = climate_data$mask,
#                        alpha = 0.1)
# 
# end_time <- Sys.time()
# time_taken <- end_time - start_time
# print(time_taken)

## ----climate_data-plot_cs, eval = FALSE---------------------------------------
# par(mfrow = c(2, 3), mar = c(3, 3, 2, 1))
# p2 = plot_cs(list(scb_up = temp$scb_up, scb_low = temp$scb_low), levels = c(1.5, 2,2.5,3), x = temp$x, y = temp$y, mu_hat = temp$mu_hat, xlab = "Longitude", ylab = "Latitude", level_label = T, min.size = 40, palette = "Spectral", color_level_label = "black")
# p1 = plot_cs(list(scb_up = temp$scb_up, scb_low = temp$scb_low), levels = c(1.5,2,2.5), x = temp$x, y = temp$y, mu_hat = temp$mu_hat, xlab = "Longitude", ylab = "Latitude",together = F)
# p = p2/p1
# p

