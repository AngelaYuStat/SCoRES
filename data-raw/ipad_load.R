## code to load `ipad` dataset goes here

library(readr)
library(dplyr)

# read Rdata
dat_path <- system.file("extdata", "earlytimepoint_recentUse_20240901.rds", package = "SCoRES")
ipad <- readRDS(dat_path)

ipad <- ipad %>%
  mutate(id = as.factor(id)) %>%
  filter(timept == 2) %>%
  tibble::as_tibble()

# save as .rda
usethis::use_data(ipad, overwrite = TRUE)
