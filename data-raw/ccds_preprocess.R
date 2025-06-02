## code to process `ccds` dataset goes here

library(here)
library(dplyr)

here::i_am(
  "data-raw/ccds_preprocess.R"
)

load(here::here("data", "ccds_raw.rda"))

# select right eye and tp = post
ccds <- ccds_raw %>%
  filter(eye == "Right", tp == "post") %>%
  dplyr::select(subject_id, seconds, use, percent_change) %>%
  mutate(subject = as.factor(subject_id)) %>%
  mutate(use = as.numeric(use == "use")) %>%
  #drop_na(use, percent_change) %>%
  tibble()

# save as .rda
usethis::use_data(ccds, overwrite = TRUE)
