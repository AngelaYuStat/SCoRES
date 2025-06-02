## code to prepare `climate_data` dataset goes here

library(usethis)
library(readr)

# read .dat file
dat_path <- system.file("extdata", "climate_data.dat", package = "invSCI")
#load(dat_path)
#ls()

#climate_data <- list(epoch = epoch, temp = temp, orog = orog, mask = mask,
  #lat = lat, lon = lon, year = year, current = current, future = future)

#usethis::use_data(climate_data, overwrite = TRUE)
