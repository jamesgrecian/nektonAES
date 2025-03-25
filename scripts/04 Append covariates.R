#####################################################################################
### Append contemporary covariates to the combined species presence/ absence data ###
#####################################################################################

# 2024-06-17

# amended from ANTSIE script
# tried to update WOA data from 18 to 23 but mixed layer depth unavailable
# already dropped species by this stage
# should the sea ice data be clipped to different time period?

# load libraries
require(tidyverse)
require(sf)
sf::sf_use_s2(FALSE)
require(raster)

# source the covariates for use in spatial predictions
covs <- readRDS("data/covariate_stack.rds")

######################
### Append to data ###
######################

# load species presence data
dat <- readRDS("data/presence_absence_data_10k.rds")

# extract environmental covariates
pulled <- raster::extract(covs, as.matrix(dat[c("x", "y")])) |> as_tibble()

# append to data frame
dat <- dat |> bind_cols(pulled)

# drop NAs
dat <- dat |> drop_na(sst, sst_grad, sal, ssh, ssh_grad, mld, sic, bat)

# output to file
saveRDS(dat, "data/presence_absence_data_10k_with_covariates_2024-06-17.rds")

# ends