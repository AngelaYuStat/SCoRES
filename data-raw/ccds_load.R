## code to load `ccds` dataset goes here

library(readr)

# read csv
ccds_raw <- read_csv("D:/FDA/data/ccds1_functional.csv")

# save as .rda
usethis::use_data(ccds_raw, overwrite = TRUE)
