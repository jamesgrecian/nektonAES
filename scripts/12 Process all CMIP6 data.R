##################################################################
### Process CMIP6 data for historical runs and SSP projections ###
##################################################################

# 2024-07-01

# Inputs here will be the netcdf files processed using the cdo bash script
# Load each file and take a seasonal average
# Output as a raster stack of covariates for each CMIP6 GCM?

# load libraries
require(tidyverse)
require(raster)
require(ncdf4)

# helper function to calculate mean
source("R/process_hist.R")

make_gcm_stack <- function(files){
  mlotst <- raster::stack(files[1])
  mlotst <- process_hist(mlotst, project = T)
  sos <- raster::stack(files[2])
  sos <- process_hist(sos, project = T)
  tos <- raster::stack(files[3])
  tos <- process_hist(tos, project = T)
  zos <- raster::stack(files[4])
  zos <- process_hist(zos, project = T)
  siconc <- raster::stack(files[5])
  siconc <- process_hist(siconc, project = T)
  
  gcm_stack <- stack(mlotst, sos, tos, zos, siconc)
  names(gcm_stack) <- c("mlotst", "sos", "tos", "zos", "siconc")
  return(gcm_stack)
}

# create processed stack for each
files <- list.files("/Volumes/Extreme SSD/CMIP6/CMIP/BCC/BCC-CSM2-MR", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/BCC-CSM2-MR_historical_covariate_stack.rds")

files <- list.files("/Volumes/Extreme SSD/CMIP6/CMIP/CAS/FGOALS-g3", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/FGOALS-g3_historical_covariate_stack.rds")

files <- list.files("/Volumes/Extreme SSD/CMIP6/CMIP/CMCC/CMCC-CM2-SR5", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/CMCC-CM2-SR5_historical_covariate_stack.rds")

files <- list.files("/Volumes/Extreme SSD/CMIP6/CMIP/CSIRO-ARCCSS/ACCESS-CM2", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/ACCESS-CM2_historical_covariate_stack.rds")

files <- list.files("/Volumes/Extreme SSD/CMIP6/CMIP/IPSL/IPSL-CM6A-LR", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/IPSL-CM6A-LR_historical_covariate_stack.rds")

files <- list.files("/Volumes/Extreme SSD/CMIP6/CMIP/MRI/MRI-ESM2-0", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/MRI-ESM2-0_historical_covariate_stack.rds")

files <- list.files("/Volumes/Extreme SSD/CMIP6/CMIP/NCAR/CESM2-WACCM", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/CESM2-WACCM_historical_covariate_stack.rds")

files <- list.files("/Volumes/Extreme SSD/CMIP6/CMIP/NUIST/NESM3", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/NESM3_historical_covariate_stack.rds")


### SSP245 ###
source("R/process_ssp.r")

make_gcm_stack <- function(files){
  mlotst <- raster::stack(files[1])
  mlotst <- process_ssp(mlotst, project = T)
  sos <- raster::stack(files[2])
  sos <- process_ssp(sos, project = T)
  tos <- raster::stack(files[3])
  tos <- process_ssp(tos, project = T)
  zos <- raster::stack(files[4])
  zos <- process_ssp(zos, project = T)
  siconc <- raster::stack(files[5])
  siconc <- process_ssp(siconc, project = T)
  
  gcm_stack <- stack(mlotst, sos, tos, zos, siconc)
  names(gcm_stack) <- c("mlotst", "sos", "tos", "zos", "siconc")
  return(gcm_stack)
}

# create processed stack for each GCM and SSP
### BCC-CSM2-MR ###
files <- list.files("/Volumes/Extreme SSD/CMIP6/ScenarioMIP/BCC/BCC-CSM2-MR/ssp245", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
plot(gcm_stack)
saveRDS(gcm_stack, "data/BCC-CSM2-MR_ssp245_covariate_stack.rds")

files <- list.files("/Volumes/Extreme SSD/CMIP6/ScenarioMIP/BCC/BCC-CSM2-MR/ssp585", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
plot(gcm_stack)
saveRDS(gcm_stack, "data/BCC-CSM2-MR_ssp585_covariate_stack.rds")

### FGOALS-g3 ###
# more than one file per covariate
files <- list.files("/Volumes/Extreme SSD/CMIP6/ScenarioMIP/CAS/FGOALS-g3/ssp245", recursive = T, pattern = "bil_1x1.nc", full.names = T)
mlotst <- raster::stack(files[1:2])
mlotst <- process_ssp(mlotst, project = T)
sos <- raster::stack(files[3:4])
sos <- process_ssp(sos, project = T)
tos <- raster::stack(files[5:6])
tos <- process_ssp(tos, project = T)
zos <- raster::stack(files[7:8])
zos <- process_ssp(zos, project = T)
siconc <- raster::stack(files[9])
siconc <- process_ssp(siconc, project = T)
gcm_stack <- stack(mlotst, sos, tos, zos, siconc)
names(gcm_stack) <- c("mlotst", "sos", "tos", "zos", "siconc")
plot(gcm_stack)
saveRDS(gcm_stack, "data/FGOALS-g3_ssp245_covariate_stack.rds")

files <- list.files("/Volumes/Extreme SSD/CMIP6/ScenarioMIP/CAS/FGOALS-g3/ssp585", recursive = T, pattern = "bil_1x1.nc", full.names = T)
mlotst <- raster::stack(files[1:2])
mlotst <- process_ssp(mlotst, project = T)
sos <- raster::stack(files[3:4])
sos <- process_ssp(sos, project = T)
tos <- raster::stack(files[5:6])
tos <- process_ssp(tos, project = T)
zos <- raster::stack(files[7:8])
zos <- process_ssp(zos, project = T)
siconc <- raster::stack(files[9])
siconc <- process_ssp(siconc, project = T)
gcm_stack <- stack(mlotst, sos, tos, zos, siconc)
names(gcm_stack) <- c("mlotst", "sos", "tos", "zos", "siconc")
plot(gcm_stack)
saveRDS(gcm_stack, "data/FGOALS-g3_ssp585_covariate_stack.rds")

### CMCC-CM2-SR5 ###
files <- list.files("/Volumes/Extreme SSD/CMIP6/ScenarioMIP/CMCC/CMCC-CM2-SR5/ssp245", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/CMCC-CM2-SR5_ssp245_covariate_stack.rds")

files <- list.files("/Volumes/Extreme SSD/CMIP6/ScenarioMIP/CMCC/CMCC-CM2-SR5/ssp585", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/CMCC-CM2-SR5_ssp585_covariate_stack.rds")

### ACCESS-CM2 ###
files <- list.files("/Volumes/Extreme SSD/CMIP6/ScenarioMIP/CSIRO-ARCCSS/ACCESS-CM2/ssp245", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/ACCESS-CM2_ssp245_covariate_stack.rds")

files <- list.files("/Volumes/Extreme SSD/CMIP6/ScenarioMIP/CSIRO-ARCCSS/ACCESS-CM2/ssp585", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/ACCESS-CM2_ssp585_covariate_stack.rds")

### IPSL-CM6A-LR ###
files <- list.files("/Volumes/Extreme SSD/CMIP6/ScenarioMIP/IPSL/IPSL-CM6A-LR/ssp245", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/IPSL-CM6A-LR_ssp245_covariate_stack.rds")

files <- list.files("/Volumes/Extreme SSD/CMIP6/ScenarioMIP/IPSL/IPSL-CM6A-LR/ssp585", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/IPSL-CM6A-LR_ssp585_covariate_stack.rds")

### MRI-ESM2-0 ###
files <- list.files("/Volumes/Extreme SSD/CMIP6/ScenarioMIP/MRI/MRI-ESM2-0/ssp245", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/MRI-ESM2-0_ssp245_covariate_stack.rds")

files <- list.files("/Volumes/Extreme SSD/CMIP6/ScenarioMIP/MRI/MRI-ESM2-0/ssp585", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/MRI-ESM2-0_ssp585_covariate_stack.rds")

### CESM2-WACCM ###
# only 2065 to 2100 so use amended function
files <- list.files("/Volumes/Extreme SSD/CMIP6/ScenarioMIP/NCAR/CESM2-WACCM/ssp245", recursive = T, pattern = "bil_1x1.nc", full.names = T)
mlotst <- raster::stack(files[1])
mlotst <- process_ssp_short(mlotst, project = T)
sos <- raster::stack(files[2])
sos <- process_ssp_short(sos, project = T)
tos <- raster::stack(files[3])
tos <- process_ssp_short(tos, project = T)
zos <- raster::stack(files[4])
zos <- process_ssp_short(zos, project = T)
siconc <- raster::stack(files[5])
siconc <- process_ssp_short(siconc, project = T)
gcm_stack <- stack(mlotst, sos, tos, zos, siconc)
names(gcm_stack) <- c("mlotst", "sos", "tos", "zos", "siconc")
plot(gcm_stack)
saveRDS(gcm_stack, "data/CESM2-WACCM_ssp245_covariate_stack.rds")

# weirdly SSP585 has data for 2015-2100
files <- list.files("/Volumes/Extreme SSD/CMIP6/ScenarioMIP/NCAR/CESM2-WACCM/ssp585", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
plot(gcm_stack)
saveRDS(gcm_stack, "data/CESM2-WACCM_ssp585_covariate_stack.rds")

### NESM ###
files <- list.files("/Volumes/Extreme SSD/CMIP6/ScenarioMIP/NUIST/NESM3/ssp245", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/NESM3_ssp245_covariate_stack.rds")

files <- list.files("/Volumes/Extreme SSD/CMIP6/ScenarioMIP/NUIST/NESM3/ssp585", recursive = T, pattern = "bil_1x1.nc", full.names = T)
gcm_stack <- make_gcm_stack(files)
saveRDS(gcm_stack, "data/NESM3_ssp585_covariate_stack.rds")

# ends