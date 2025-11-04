## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  message   = FALSE,
  warning   = FALSE,
  fig.width = 6
)

## ----load pupil data, echo = TRUE---------------------------------------------
library(SCoRES)
library(mgcv)
data(pupil)

## ----fit gam-fpca model-1, echo = TRUE----------------------------------------
pupil_fpca <- SCoRES::prepare_pupil_fpca(pupil)
fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
            s(seconds, by = use, k=30, bs = "cr") +
            s(id, by = Phi1, bs="re") +
            s(id, by = Phi2, bs="re") +
            s(id, by = Phi3, bs="re") +
            s(id, by = Phi4, bs="re"),
            method = "fREML", data = pupil_fpca, discrete = TRUE)

## ----pupil_plot_cs_cma, echo = TRUE, eval = FALSE-----------------------------
# # CMA approach
# results_pupil_cma <- SCoRES::SCB_functional_outcome(
#                                           data_df = pupil,
#                                           object = fosr_mod,
#                                           method = "cma",
#                                           fitted = TRUE,
#                                           alpha = 0.05,
#                                           outcome = "percent_change",
#                                           domain = "seconds",
#                                           subset = c("use = 1"),
#                                           id = "id")
# 
# results_pupil_cma <- tibble::as_tibble(results_pupil_cma)
# plot_cs(results_pupil_cma,
#         levels = c(-18, -20, -22, -24),
#         x = results_pupil_cma$domain,
#         mu_hat = results_pupil_cma$mu_hat,
#         xlab = "Seconds",
#         ylab = "Percent_Outcome",
#         level_label = T,
#         min.size = 40,
#         palette = "Spectral",
#         color_level_label = "black")

## ----ccds_plot_cs_pupil_para, echo = TRUE, eval = FALSE-----------------------
# #CMA approach for parameter function
# results_pupil_cma_para <- SCoRES::SCB_functional_outcome(
#                                           data_df = pupil,
#                                           object = fosr_mod,
#                                           method = "cma",
#                                           fitted = FALSE,
#                                           alpha = 0.05,
#                                           outcome = "percent_change",
#                                           domain = "seconds",
#                                           subset = c("use = 1"),
#                                           id = "id")
# 
# results_pupil_cma_para <- tibble::as_tibble(results_pupil_cma_para)
# plot_cs(results_pupil_cma_para,
#         levels = c(4.5, 5, 5.5, 6),
#         x = results_pupil_cma_para$domain,
#         mu_hat = results_pupil_cma_para$mu_hat,
#         xlab = "Seconds",
#         ylab = "Percent_Outcome",
#         level_label = T,
#         min.size = 40,
#         palette = "Spectral",
#         color_level_label = "black")

## ----pupil_plot_cs_multiplier, echo = TRUE------------------------------------
# Multiplier-t Bootstrap
results_pupil_multiplier <- SCoRES::SCB_functional_outcome(
                                          data_df = pupil,
                                          object = fosr_mod, 
                                          method = "multiplier",
                                          fitted = TRUE, 
                                          alpha = 0.05, 
                                          outcome = "percent_change",
                                          domain = "seconds", 
                                          subset = c("use = 1"), 
                                          id = "id")

results_pupil_multiplier <- tibble::as_tibble(results_pupil_multiplier)
plot_cs(results_pupil_multiplier,
        levels = c(-18, -20, -22, -24), 
        x = results_pupil_multiplier$domain, 
        mu_hat = results_pupil_multiplier$mu_hat, 
        xlab = "Seconds", 
        ylab = "Percent_Outcome",
        level_label = T, 
        min.size = 40, 
        palette = "Spectral", 
        color_level_label = "black")

## ----pupil_plot_cs_multiplier_para, echo = TRUE-------------------------------
results_pupil_multiplier_para <- SCoRES::SCB_functional_outcome(
                                          data_df = pupil,
                                          object = fosr_mod, 
                                          method = "multiplier",
                                          fitted = FALSE, 
                                          alpha = 0.05, 
                                          outcome = "percent_change",
                                          domain = "seconds", 
                                          subset = c("use = 1"), 
                                          id = "id")

results_pupil_multiplier_para <- tibble::as_tibble(results_pupil_multiplier_para)
plot_cs(results_pupil_multiplier_para,
        levels = c(4.5, 5, 5.5, 6), 
        x = results_pupil_multiplier_para$domain, 
        mu_hat = results_pupil_multiplier_para$mu_hat, 
        xlab = "Seconds", 
        ylab = "Percent_Outcome",
        level_label = T, 
        min.size = 40, 
        palette = "Spectral", 
        color_level_label = "black")

## ----pupil_plot_cs_multiplier_overall, echo = TRUE, eval = FALSE--------------
# results_pupil_multiplier_overall <- SCoRES::SCB_functional_outcome(
#                                           data_df = pupil,
#                                           method = "multiplier",
#                                           fitted = TRUE,
#                                           alpha = 0.05,
#                                           outcome = "percent_change",
#                                           domain = "seconds",
#                                           subset = c("use = 1"),
#                                           id = "id")
# 
# results_pupil_multiplier_overall <- tibble::as_tibble(results_pupil_multiplier_overall)
# plot_cs(results_pupil_multiplier_overall,
#         levels = c(-18, -20, -22, -24),
#         x = results_pupil_multiplier_overall$domain,
#         mu_hat = results_pupil_multiplier_overall$mu_hat,
#         xlab = "Seconds",
#         ylab = "Percent_Outcome",
#         level_label = T,
#         min.size = 40,
#         palette = "Spectral",
#         color_level_label = "black")

## ----fit gam-fpca model-2, echo = TRUE----------------------------------------
pupil_fpca <- SCoRES::prepare_pupil_fpca(pupil, example = "extended")
fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
            s(seconds, by = use, k=30, bs = "cr") +
            s(seconds, by = age, k = 30, bs = "cr") +
            s(seconds, by = gender, k = 30, bs = "cr") +
            s(id, by = Phi1, bs="re") +
            s(id, by = Phi2, bs="re") +
            s(id, by = Phi3, bs="re") +
            s(id, by = Phi4, bs="re"),
            method = "fREML", data = pupil_fpca, discrete = TRUE)

## ----pupil_plot_cs_cma_2, echo = TRUE, eval = FALSE---------------------------
# # CMA approach
# results_pupil_cma <- SCoRES::SCB_functional_outcome(
#                                           data_df = pupil,
#                                           object = fosr_mod,
#                                           method = "cma",
#                                           fitted = TRUE,
#                                           alpha = 0.05,
#                                           outcome = "percent_change",
#                                           domain = "seconds",
#                                           subset = c("use = 1",
#                                                      "age = 40",
#                                                      "gender = 0"),
#                                           id = "id")
# 
# results_pupil_cma <- tibble::as_tibble(results_pupil_cma)
# plot_cs(results_pupil_cma,
#         levels = c(-18, -20, -22, -24),
#         x = results_pupil_cma$domain,
#         mu_hat = results_pupil_cma$mu_hat,
#         xlab = "Seconds",
#         ylab = "Percent_Outcome",
#         level_label = T,
#         min.size = 40,
#         palette = "Spectral",
#         color_level_label = "black")

## ----pupil_plot_cs_multiplier_2, echo = TRUE----------------------------------
# Multiplier-t Bootstrap
results_pupil_multiplier <- SCoRES::SCB_functional_outcome(
                                          data_df = pupil,
                                          object = fosr_mod,
                                          fitted = TRUE,
                                          method = "multiplier", 
                                          alpha = 0.05, 
                                          outcome = "percent_change", 
                                          domain = "seconds", 
                                          subset = c("use = 1", 
                                                     "age = 40", 
                                                     "gender = 0"), 
                                          id = "id")

results_pupil_multiplier <- tibble::as_tibble(results_pupil_multiplier)
plot_cs(results_pupil_multiplier,
        levels = c(-18, -20, -22, -24), 
        x = results_pupil_multiplier$domain, 
        mu_hat = results_pupil_multiplier$mu_hat, 
        xlab = "Seconds", 
        ylab = "Percent_Outcome",
        level_label = T, 
        min.size = 40, 
        palette = "Spectral", 
        color_level_label = "black")

