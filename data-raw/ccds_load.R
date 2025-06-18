## code to load `ccds` dataset goes here

library(readr)

# read csv
dat_path <- system.file("extdata", "ccds1_functional.csv", package = "invSCI")
ccds_raw <- read_csv(dat_path)

# select right eye and tp = post
ccds <- ccds_raw %>%
  filter(eye == "Right", tp == "post") %>%
  dplyr::select(subject_id, seconds, use, percent_change) %>%
  mutate(subject = as.factor(subject_id)) %>%
  mutate(use = as.numeric(use == "use")) %>%
  #drop_na(use, percent_change) %>%
  tibble::as_tibble()

# save as .rda
usethis::use_data(ccds, overwrite = TRUE)
