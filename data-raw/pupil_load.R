## code to load `ccds` dataset goes here

library(readr)
library(dplyr)

# read Rdata
dat_path <- system.file("extdata", "pupil.Rdata", package = "invSCI")
load(dat_path)
pupil_raw <- pupil

pupil <- pupil_raw %>%
  mutate(id = as.factor(id)) %>%
  mutate(gender = as.numeric(gender == "Female")) %>%
  mutate(smoker = ifelse(use_group %in% c("Daily - Concentrates",
                                          "Daily - Flower",
                                          "Occasional - Flower"), 1, 0)) %>%
  mutate(daily = ifelse(use_group %in% c("Daily - Concentrates",
                                         "Daily - Flower"), 1, 0)) %>%
  tibble::as_tibble()

# save as .rda
usethis::use_data(pupil, overwrite = TRUE)
