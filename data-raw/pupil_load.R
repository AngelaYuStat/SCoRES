## code to load `ccds` dataset goes here

library(readr)
library(dplyr)

# read Rdata
dat_path <- system.file("extdata", "pupil.Rdata", package = "SCoRES")
load(dat_path)
pupil_raw <- pupil

pupil <- pupil_raw %>%
  mutate(id = as.factor(id)) %>%
  mutate(gender = as.numeric(gender == "Female")) %>%
  tibble::as_tibble()

pupil <- pupil %>%
  filter(!(id %in% c("003-1068", "003-112")))

# save as .rda
usethis::use_data(pupil, overwrite = TRUE)
