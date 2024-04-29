#####################################################################################
### Append contemporary covariates to the combined species presence/ absence data ###
#####################################################################################

# 2024-04-25

# amended from ANTSIE script
# tried to update WOA data from 18 to 23 but mixed layer depth unavailable
# already dropped species by this stage
# should the sea ice data be clipped to different time period?

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

######################
### Append to data ###
######################

# load species presence data
dat <- readRDS("data/presence_absence_data_10k.rds")

# extract environmental covariates
pulled <- raster::extract(covs, as.matrix(dat[c("x", "y")])) %>% as_tibble()

# append to data frame
dat <- dat |> bind_cols(pulled)

# drop NAs
dat <- dat |> drop_na(sst, sst_grad, sal, ssh, ssh_grad, mld, sic, bat)

# output to file
saveRDS(dat, "data/presence_absence_data_10k_with_covariates_2024-04-25.rds")

# ends