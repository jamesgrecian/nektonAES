##############################
### Change Factor Protocol ###
##############################

# 2023-08-03 first draft
# 2024-07-03 update for CMIP6 SSP data

# Function to generate the covariates for future climate projections
# Input is:
# 1. a covariate stack of observations - the one used to fit the distribution models
# 2. Name of a model - this will be pattern matched so partial should work
# Output will be a covariate stack that can be passed to model prediction functions

# Idea based on change factor protocol:
# Tabor & Williams (2010) Globally downscaled climate projections for assessing the conservation impacts of climate change. Ecol Appl 20, 554–565
# Lorenz et al. (2016) Downscaled and debiased climate simulations for North America from 21,000 years ago to 2100AD. Sci. Data 3, 160048

# Function will calculate the climate model anomaly
# the difference between scenario of interest and the Historical period:
# anomaly = ssp - historical

# This anomaly is then added to the observed environmental conditions:
# adjusted = observations + anomaly

change_factor <- function(covariate_stack, model, scenario = "ssp245"){

  # source GCM model data
  files <- list.files("data", pattern = model, full.names = T)
  historical <- readRDS(files[1])
  ssp245 <- readRDS(files[2])
  ssp585 <- readRDS(files[3])
  
  # some GCM sea ice concentration data has NA strip
  historical$siconc <- focal(historical$siconc, w = matrix(1, nrow = 3, ncol = 3), fun = mean, NAonly = TRUE, na.rm = TRUE) # replace NAs with focal neighborhood mean
  historical$siconc <- mask(historical$siconc, historical$mlotst) # interpolation leaks onto land so mask by a GCM without the NA issue
  ssp245$siconc <- focal(ssp245$siconc, w = matrix(1, nrow = 3, ncol = 3), fun = mean, NAonly = TRUE, na.rm = TRUE) # replace NAs with focal neighborhood mean
  ssp245$siconc <- mask(ssp245$siconc, ssp245$mlotst) # interpolation leaks onto land so mask by a GCM without the NA issue
  ssp585$siconc <- focal(ssp585$siconc, w = matrix(1, nrow = 3, ncol = 3), fun = mean, NAonly = TRUE, na.rm = TRUE) # replace NAs with focal neighborhood mean
  ssp585$siconc <- mask(ssp585$siconc, historical$mlotst) # interpolation leaks onto land so mask by a GCM without the NA issue
  
  # anomaly = ssp - historical
  if(scenario == "ssp245"){anomaly <- ssp245 - historical}
  if(scenario == "ssp585"){anomaly <- ssp585 - historical} 
  
  names(anomaly) <- c("mld", "sal", "sst", "ssh", "sic") # rename to match cov stack names
  
  # store bathymetry raster to reappend to stack after change factor protocol
  bat <- covariate_stack$bat
  
  # ensure stacks are in correct order before adding anomaly and observations together
  covariate_stack <- subset(covariate_stack, names(anomaly)) # also drop gradients as we need to recalculate these
  
  # align the grids
  anomaly <- projectRaster(anomaly, covariate_stack) # align the grids
  anomaly <- mask(anomaly, mask = covariate_stack$sic, maskvalue = NA)
  
  # add anomaly to the observations for change factor outcome
  covs <- covariate_stack + anomaly
  
  # bound the covariates
  covs$mld[covs$mld < 0] <- 0 # remove mld values above water
  covs$sic[covs$sic < 0] <- 0 # remove negative sea ice concentrations
  covs$sal[covs$sal < 30] <- NA # remove salinities below 30 ppt
  
  # create slope covariate
  sst_grad <- terrain(covs$sst, opt = "slope", unit = "radians")
  ssh_grad <- terrain(covs$ssh, opt = "slope", unit = "radians")
  
  # create output stack
  output_stack <- stack(covs$sst, sst_grad, covs$sal, covs$ssh, ssh_grad, covs$mld, covs$sic, bat)
  names(output_stack) <- c("sst", "sst_grad", "sal", "ssh", "ssh_grad", "mld", "sic", "bat")
  
  #plot(output_stack)
  
  return(output_stack)
}


