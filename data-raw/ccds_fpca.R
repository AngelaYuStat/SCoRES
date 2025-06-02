## code to prepare `ccds_fpca` dataset goes here

library(here)
library(dplyr)
library(mgcv)
library(tidyr)
library(tibble)
library(refund)

here::i_am(
  "data-raw/ccds_fpca.R"
)

load(here::here("data", "ccds.rda"))

# fit the initial mean model
mean_mod = mgcv::gam(percent_change ~ s(seconds, k=30, bs="cr") +
                       s(seconds, by=use, k=30, bs = "cr"),
                     data = ccds, method = "REML")

# obtain residuals from the mean model
resid_df = ccds %>%
  filter(!is.na(percent_change)) %>%
  dplyr::select(subject, seconds) %>%
  mutate(resid = mean_mod$residuals) %>%
  pivot_wider(names_from = seconds, values_from = resid, names_prefix = "resid.")

resid_mat = as.matrix(resid_df[,-1])
rownames(resid_mat) = resid_df$subject

# estimate the FPCA for GAMM FPCA model
fpca_results = fpca.face(resid_mat, argvals = unique(ccds$seconds), knots = 15)
eigenfunctions = as.data.frame(fpca_results$efunctions)
colnames(eigenfunctions) = paste0("Phi", seq(1, fpca_results$npc))
eigenfunctions$seconds = unique(ccds$seconds)
ccds_fpca =  ccds %>% left_join(., eigenfunctions, by = "seconds") %>%
  as_tibble() %>%
  arrange(subject, seconds) %>%
  mutate(subject = factor(subject))

usethis::use_data(ccds_fpca, overwrite = TRUE)
