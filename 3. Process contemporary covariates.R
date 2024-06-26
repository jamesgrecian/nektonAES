#########################################
### Processing contemporary data sets ###
#########################################

# 2024-06-17

# Characterise contemporary environmental conditions
# These will be used to drive species distribution models
# Use climatological means over Austral summer: October - March

# load libraries
require(tidyverse)
require(sf)
sf::sf_use_s2(FALSE)
require(raster)


################
### WOA 2018 ###
################

# Temperature, salinity and mixed layer depth available from World Ocean Atlas 2018 database (Boyer et al. 2018).
# https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/
# https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DOC/woa18documentation.pdf
# t_an is objectively analyzed climatological mean
# level 1 is 0m depth
# download monthly files and then combine to create seasonal mean

# WOA temperature
fn <- list.files("~/Dropbox/ANTSIE climate data/WOA2018/temperature", full.names = T)
temp <- stack()
for (i in 1:length(fn)){
  foo <- brick(fn[i], varname = "t_an")
  foo <- subset(foo, 1)
  temp <- stack(temp, foo)
}
temp <- mean(temp)
plot(temp)
writeRaster(temp, "data/temperature_025_aus_mean_wgs84", overwrite = TRUE)

# WOA salinity
fn <- list.files("~/Dropbox/ANTSIE climate data/WOA2018/salinity", full.names = T)
sal <- stack()
for (i in 1:length(fn)){
  foo <- brick(fn[i], varname = "s_an")
  foo <- subset(foo, 1)
  sal <- stack(sal, foo)
}
sal <- mean(sal)
plot(sal)
writeRaster(sal, "data/salinity_025_aus_mean_wgs84", overwrite = TRUE)

# WOA mixed-layer depth
fn <- list.files("~/Dropbox/ANTSIE climate data/WOA2018/mixed layer depth", full.names = T)
mld <- stack()
for (i in 1:length(fn)){
  foo <- brick(fn[i], varname = "M_an")
  foo <- subset(foo, 1)
  mld <- stack(mld, foo)
}
mld <- mean(mld)
plot(mld)
writeRaster(mld, "data/mixedlayerdepth_025_aus_mean_wgs84", overwrite = TRUE)


#############
### AVISO ###
#############

# Sea surface height data were extracted from SSALTO/ DUACS and
# AVISO Mean sea level height above geoid
# https://www.aviso.altimetry.fr/index.php?id=1526

fn <- list.files("~/Dropbox/ANTSIE climate data/AVISO/", full.names = T)
ssh <- raster::stack(fn)
ssh <- rotate(ssh)
ssh <- mean(ssh, na.rm = T)
plot(ssh)
writeRaster(ssh, "data/surfaceheight_025_aus_mean_wgs84", overwrite = TRUE)


#############
### NSIDC ###
#############

# recalculate the sea ice climatology data
# use the 1981-2010 average to align with the WOA2018 data

#load libraries
require(tidyverse)
require(stringr)
require(raster)
source("R/monthly_ice_average.R")

# function downloads all the daily data for the years and month specified
# returns the average across those years for the focal month
nsidc_01 <- monthly_ice_average(year_min = 1981, year_max = 2010, m = 1)
nsidc_02 <- monthly_ice_average(year_min = 1981, year_max = 2010, m = 2)
nsidc_03 <- monthly_ice_average(year_min = 1981, year_max = 2010, m = 3)
nsidc_10 <- monthly_ice_average(year_min = 1981, year_max = 2010, m = 10)
nsidc_11 <- monthly_ice_average(year_min = 1981, year_max = 2010, m = 11)
nsidc_12 <- monthly_ice_average(year_min = 1981, year_max = 2010, m = 12)
# no data for December 1987 and January 1988 due to satellite problems

sic <- stack(nsidc_01, nsidc_02, nsidc_03, nsidc_10, nsidc_11, nsidc_12)
names(sic) <- c("Jan", "Feb", "Mar", "Oct", "Nov", "Dec")

# output netcdf file for Rahul with monthly averages
writeRaster(sic, "data/seaice_NSIDC_1981_2010_Oct_Mar.nc", format = "CDF")

# output R raster
writeRaster(sic, "data/seaice_025_aus_mean_wgs84_months", overwrite = TRUE)
writeRaster(mean(sic), "data/seaice_025_aus_mean_wgs84", overwrite = TRUE)


################################################
### Store contemporary covariates as a stack ###
################################################

# load libraries
require(tidyverse)
require(sf)
sf::sf_use_s2(FALSE)
require(raster)

# load contemporary environmental covariates
sst <- raster("data/temperature_025_aus_mean_wgs84")
sal <- raster("data/salinity_025_aus_mean_wgs84")
ssh <- raster("data/surfaceheight_025_aus_mean_wgs84")
mld <- raster("data/mixedlayerdepth_025_aus_mean_wgs84")
sic <- raster("data/seaice_025_aus_mean_wgs84") # sea ice is already polar projection
sst_grad <- terrain(sst, opt = "slope", unit = "radians")
ssh_grad <- terrain(ssh, opt = "slope", unit = "radians")

# mixed layer depth raster only extends to 78.5 S so extent to match
mld <- extend(mld, sst, value = NA)

# stack the covariates into single object - not sic as it's already projected
covs <- stack(sst, sst_grad, sal, ssh, ssh_grad, mld)
names(covs) <- c("sst", "sst_grad", "sal", "ssh", "ssh_grad", "mld")

# project to southern hemisphere polar projection
# need to clip to southern hemisphere first
ext <- extent(-180, 180, -90, -40)
covs <- crop(covs, ext)
template <- projectRaster(from = covs, to = sic, alignOnly = T) # make a template that aligns to the sic grid but wider for whole southern ocean
covs <- projectRaster(from = covs, to = template) 

sic <- extend(sic, covs, value = 0)
sic <- mask(sic, mask = subset(covs, 1), maskvalue = NA)
names(sic) <- "sic"
covs <- stack(covs, sic)

# download some bathymetry data as an example from marmap
# use 10 minute resolution as proxy of <0.25 degrees
# this can be downgraded to 25 km x 25 km
bathy <- marmap::getNOAA.bathy(lon1 = -180,
                               lon2 = 180,
                               lat1 = -90,
                               lat2 = -40,
                               resolution = 10)
bat <- marmap::as.raster(bathy)
bat <- projectRaster(from = bat, to = covs)
bat <- mask(bat, mask = subset(covs, 1), maskvalue = NA)
names(bat) <- "bat"
covs <- stack(covs, bat)
plot(covs)

# output the covariates for use in spatial predictions
saveRDS(covs, "data/covariate_stack.rds")

# ends