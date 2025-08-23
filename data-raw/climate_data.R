## code to prepare `climate_data` dataset goes here

library(usethis)
library(readr)

# read .dat file
dat_path <- system.file("extdata", "climate_data.dat", package = "SCoRES")
attach(temp)
load(dat_path)
ls()

#climate_data <- list(epoch = epoch, temp = temp, orog = orog, mask = mask,
  #lat = lat, lon = lon, year = year, current = current, future = future)

mask = mask * ifelse(orog>0, 1, NA)
mask = mask[17:215,22:110]

#Define the t vectors and n.
ta = current
ta = ta - mean(ta)
tb = future
tb = tb - mean(tb)

na = length(ta)
nb = length(tb)
n = na + nb

#Define Design matrix X.
X1 = c(rep(0,na),rep(1,nb))
X2 = rep(1,n)
X3 = c(ta,rep(0,nb))
X4 = c(rep(0,na),tb)
X = cbind(X1,X2,X3,X4)

#Define the data Y.
lon <- lon[17:215]
lat <- lat[22:110]
Y = array(0,c(length(lon),length(lat),n))
for(j in 1:na) Y[,,j] = summer[17:215,22:110,j,1]
for(j in 1:nb) Y[,,na+j] = summer[17:215,22:110,j,2]
#for(j in 1:na) Y[,,j] = winter[17:215,22:110,j,1]
#for(j in 1:nb) Y[,,na+j] = winter[17:215,22:110,j,2]
Z = list(x = lon, y = lat, obs = Y)
correlation = "corAR1"

climate_data <- list(Z = Z, mask = mask, X = X, correlation = correlation)
usethis::use_data(climate_data, overwrite = TRUE)
