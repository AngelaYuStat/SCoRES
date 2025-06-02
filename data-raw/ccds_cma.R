## code to prepare `ccds_cma` dataset goes here

library(here)
library(dplyr)
library(mgcv)
library(tidyr)
library(tibble)
library(refund)

here::i_am(
  "data-raw/ccds_cma.R"
)

load(here::here("data", "ccds_fpca.rda"))

# obtain the functional regression model object
fosr_mod <- mgcv::bam(percent_change ~ s(seconds, k=30, bs="cr") +
                        s(seconds, by = use, k=30, bs = "cr") +
                        s(subject, by = Phi1, bs="re") +
                        s(subject, by = Phi2, bs="re")+
                        s(subject, by = Phi3, bs="re") +
                        s(subject, by = Phi4, bs="re"),
                        method = "fREML", data = ccds_fpca, discrete = TRUE)

ccds_cma <- cma(fosr_mod, "gam", "seconds", groups = "use", subject = "subject")

usethis::use_data(ccds_cma, overwrite = TRUE)
