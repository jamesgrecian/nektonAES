#######################################################################
### Create sea ice climatology by averaging over years for each day ###
#######################################################################

# 2024-05-14

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



