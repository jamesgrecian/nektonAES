####################################################
### Process CMIP6 data for comparison with NSIDC ###
####################################################

# 2024-06-26

# Previously downloaded historical sea ice data from 15 GCMs
# These have been transformed to 1x1 degree grid using cdo bash script
# Extract data for 1980 to 2010 and take Austral summer mean

# load libraries
require(tidyverse)
require(raster)
require(ndcf4)

# helper function to calculate mean
source("R/process_hist.R")

### Process data from each GCM ###

# ACCESS-CM2
access_cm2 <- raster::stack("/Users/home/ANTSIE data/CMIP6/CMIP/CSIRO-ARCCSS/ACCESS-CM2/historical/r1i1p1f1/SImon/siconc/gn/20200817/siconc_SImon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc")
access_cm2 <- process_hist(access_cm2, project = T)

# ACCESS-ESM1-5
access_esm15 <- raster::stack("/Users/home/ANTSIE data/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/SImon/siconc/gn/20200817/siconc_SImon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc")
access_esm15 <- process_hist(access_esm15, project = T)

# BCC-CSM2-MR
bcc_csm2_mr <- raster::stack("/Users/home/ANTSIE data/CMIP6/CMIP/BCC/BCC-CSM2-MR/historical/r1i1p1f1/SImon/siconc/gn/20200218/siconc_SImon_BCC-CSM2-MR_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc")
bcc_csm2_mr <- process_hist(bcc_csm2_mr, project = T)

# CanESM5
canesm5 <- raster::stack("/Users/home/ANTSIE data/CMIP6/CMIP/CCCma/CanESM5/historical/r1i1p1f1/SImon/siconc/gn/20190429/siconc_SImon_CanESM5_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc")
canesm5 <- process_hist(canesm5, project = T)

# CESM2-WACCM
cesm2_waccm <- raster::stack("/Users/home/ANTSIE data/CMIP6/CMIP/NCAR/CESM2-WACCM/historical/r1i1p1f1/SImon/siconc/gn/20190227/siconc_SImon_CESM2-WACCM_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc")
cesm2_waccm <- process_hist(cesm2_waccm, project = T) # leaves a line of NAs?

# CMCC-CM2-SR5
cmcc_cm2_sr5 <- raster::stack("/Users/home/ANTSIE data/CMIP6/CMIP/CMCC/CMCC-CM2-SR5/historical/r1i1p1f1/SImon/siconc/gn/20200616/siconc_SImon_CMCC-CM2-SR5_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc")
cmcc_cm2_sr5 <- process_hist(cmcc_cm2_sr5, project = T)

#EC-Earth3
ec_earth3 <- list.files("/Users/home/ANTSIE data/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3/historical/r1i1p1f1/SImon/siconc/gn/20200918", pattern = "bil_1x1", full.names = T)
ec_earth3 <- stack(ec_earth3)
ec_earth3 <- process_hist(ec_earth3, project = T)

#EC-Earth3-Veg
ec_earth3_veg <- list.files("/Users/home/ANTSIE data/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3-Veg/historical/r1i1p1f1/SImon/siconc/gn/20211207", pattern = "bil_1x1", full.names = T)
ec_earth3_veg <- stack(ec_earth3_veg)
ec_earth3_veg <- process_hist(ec_earth3_veg, project = T)

# FGOALS-g3
fgoals_g3 <- raster::stack("/Users/home/ANTSIE data/CMIP6/CMIP/CAS/FGOALS-g3/historical/r1i1p1f1/SImon/siconc/gn/20210108/siconc_SImon_FGOALS-g3_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc")
fgoals_g3 <- process_hist(fgoals_g3, project = T)

# IPSL-CM6A-LR
ipsl_cm6a_lr <- raster::stack("/Users/home/ANTSIE data/CMIP6/CMIP/IPSL/IPSL-CM6A-LR/historical/r1i1p1f1/SImon/siconc/gn/20180803/siconc_SImon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc")
ipsl_cm6a_lr <- process_hist(ipsl_cm6a_lr, project = T)

# MIROC6
miroc6 <- list.files("/Users/home/ANTSIE data/CMIP6/CMIP/MIROC/MIROC6/historical/r1i1p1f1/SImon/siconc/gn/20181212", pattern = "bil_1x1", full.names = T)
miroc6 <- stack(miroc6)
miroc6 <- process_hist(miroc6, project = T)

# MPI-ESM1-2-LR
mpi_esm1_2_lr <- list.files("/Users/home/ANTSIE data/CMIP6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/SImon/siconc/gn/20190710", pattern = "bil_1x1", full.names = T)
mpi_esm1_2_lr <- raster::stack(mpi_esm1_2_lr)
mpi_esm1_2_lr <- process_hist(mpi_esm1_2_lr, project = T)

# MRI-ESM2-0
mri_esm2_0 <- raster::stack("/Users/home/ANTSIE data/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/SImon/siconc/gn/20190904/siconc_SImon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc")
mri_esm2_0 <- process_hist(mri_esm2_0, project = T)

# NESM3
nesm3 <- raster::stack("/Users/home/ANTSIE data/CMIP6/CMIP/NUIST/NESM3/historical/r1i1p1f1/SImon/siconc/gn/20190704/siconc_SImon_NESM3_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc")
nesm3 <- process_hist(nesm3, project = T)

# NorESM2-LM
noresm2_lm <- list.files("/Users/home/ANTSIE data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/SImon/siconc/gn/20190815", pattern = "bil_1x1", full.names = T)
noresm2_lm <- raster::stack(noresm2_lm)
noresm2_lm <- process_hist(noresm2_lm, project = T)

### stack rasters ###
historical_siconc <- raster::stack(access_cm2, access_esm15, bcc_csm2_mr, canesm5,
                                   cesm2_waccm, cmcc_cm2_sr5, ec_earth3, ec_earth3_veg,
                                   fgoals_g3, ipsl_cm6a_lr, miroc6, mpi_esm1_2_lr,
                                   mri_esm2_0, nesm3, noresm2_lm)
names(historical_siconc) <- c("ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5",
                              "CESM2-WACCM", "CMCC-CM2-SR5", "EC-Earth3", "EC-Earth3-Veg",
                              "FGOALS-G3", "IPSL-CM6A-LR", "MIROC6", "MPI_ESM1-2-LR",
                              "MRI_ESM2-0", "NESM3", "NorESM2-LM")

# load in NSIDC observations for comparison
covs <- readRDS("data/covariate_stack.rds")
sic <- subset(covs, "sic")
sic <- projectRaster(sic, historical_siconc) # project sea ice to match historical run

# mask historical run to same extent as sic data
historical_siconc <- mask(historical_siconc, mask = sic, maskvalue = NA)

######################
### Calculate RMSE ###
######################

rasterRMSE <- function(observed, predicted){
  rmse <- observed - predicted
  rmse <- rmse^2
  rmse <- cellStats(rmse, mean)
  rmse <- sqrt(rmse)
  return(rmse)
}

# calculate the RMSE for each layer compared to the observed sea ice data
rasterRMSE(sic, historical_siconc) 

# as a nice tibble
rmse <- tibble(model = names(historical_siconc),
               rmse = round(rasterRMSE(sic, historical_siconc), 2))
rmse |> arrange(rmse)

# create object for plotting
rmse <- tibble(model = names(historical_siconc),
               rmse = paste("RMSE =", round(rasterRMSE(sic, historical_siconc), 2)))

rmse <- rmse |> 
  mutate(model = recode_factor(model,
                               ACCESS.CM2 = "ACCESS-CM2",
                               ACCESS.ESM1.5 = "ACCESS-ESM1-5",
                               BCC.CSM2.MR = "BCC-CSM2-MR",
                               CanESM5 = "CanESM5",
                               CESM2.WACCM = "CESM2-WACCM",
                               CMCC.CM2.SR5 = "CMCC-CM2-SR5",
                               EC.Earth3 = "EC-Earth3",
                               EC.Earth3.Veg = "EC-Earth3-Veg",
                               FGOALS.G3 = "FGOALS-G3",
                               IPSL.CM6A.LR = "IPSL-CM6A-LR",
                               MIROC6 = "MIROC6",
                               MPI_ESM1.2.LR = "MPI-ESM1-2-LR",
                               MRI_ESM2.0 = "MRI-ESM2-0",
                               NESM3 = "NESM3",
                               NorESM2.LM = "NorESM2-LM"))

############################
### Generate pretty plot ###
############################

# load libraries
require(sf)
sf::sf_use_s2(FALSE)

# pretty gradient and polar mask functions
source("R/discrete_gradient.R")
source("R/polar_mask.R")

# format data for plotting
ice_df <- stack(sic, historical_siconc)
names(ice_df)[1] <- "NSIDC"

ice_df <- ice_df |> 
  rasterToPoints() |> 
  as_tibble() |>
  pivot_longer(3:18,
               names_to = "model",
               values_to = "concentration") |>
  mutate(model = factor(model, 
                        levels = c("NSIDC", "ACCESS.CM2", "ACCESS.ESM1.5", "BCC.CSM2.MR",
                                   "CanESM5", "CESM2.WACCM", "CMCC.CM2.SR5", "EC.Earth3",
                                   "EC.Earth3.Veg", "FGOALS.G3", "IPSL.CM6A.LR",  "MIROC6",
                                   "MPI_ESM1.2.LR", "MRI_ESM2.0", "NESM3", "NorESM2.LM")))

# relabel to aid interpretation in plot
ice_df <- ice_df |> 
  mutate(model = recode_factor(model,
                               NSIDC = "NSIDC",
                               ACCESS.CM2 = "ACCESS-CM2",
                               ACCESS.ESM1.5 = "ACCESS-ESM1-5",
                               BCC.CSM2.MR = "BCC-CSM2-MR",
                               CanESM5 = "CanESM5",
                               CESM2.WACCM = "CESM2-WACCM",
                               CMCC.CM2.SR5 = "CMCC-CM2-SR5",
                               EC.Earth3 = "EC-Earth3",
                               EC.Earth3.Veg = "EC-Earth3-Veg",
                               FGOALS.G3 = "FGOALS-G3",
                               IPSL.CM6A.LR = "IPSL-CM6A-LR",
                               MIROC6 = "MIROC6",
                               MPI_ESM1.2.LR = "MPI-ESM1-2-LR",
                               MRI_ESM2.0 = "MRI-ESM2-0",
                               NESM3 = "NESM3",
                               NorESM2.LM = "NorESM2-LM"))

# buffer for plot
polar_buffer <- polar_mask(radius_size = 5750000)

# define projection
crs_polar <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# load in shapefile for background mask, clip and project
world_shp <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
CP <- sf::st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 0), crs = 4326) |> sf::st_as_sfc()
world_shp <- world_shp |> sf::st_crop(CP)
world_shp <- world_shp |> sf::st_transform(crs_polar)
world_shp <- world_shp |> st_union()

p1 <- ggplot() + 
  theme_void(base_size = 10,
             base_family = "Helvetica Neue") +
  geom_raster(aes(x = x, y = y, fill = concentration), data = ice_df) +
  geom_sf(aes(), data = world_shp, colour = "grey60", fill = "grey60") +
  geom_sf(aes(), data = polar_buffer$mask, fill = "white", color = NA) +
  geom_sf(aes(), data = polar_buffer$buffer, fill = NA, color = "grey40", size = 0.5 / .pt) +
  coord_sf(xlim = c(-6400000, 6400000),
           ylim = c(-6400000, 6400000),
           expand = FALSE,
           crs = crs_polar,
           ndiscr = 1000) +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~model) +
  scale_fill_discrete_gradient("Sea Ice Concentration (%)",
                               colours = rev(viridis::viridis(10)),
                               bins = 10,
                               limits = c(0.1, 100),
                               breaks = c(0.1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
                               labels = seq(0, 100, 10),
                               na.value = "white",
                               guide = guide_colourbar(nbin = 500,
                                                       #raster = T,
                                                       frame.colour = "grey40",
                                                       ticks.colour = "grey40",
                                                       frame.linewidth = .1,
                                                       barwidth = 20,
                                                       barheight = .5,
                                                       direction = "horizontal",
                                                       title.position = "top", #or "right"
                                                       title.theme = element_text(
                                                         hjust = 0.5,
                                                         size = 10,
                                                         colour = "black"))
  ) +
  theme(plot.title = element_text(
    hjust = 0.5,
    size = 10,
    family = "Helvetica Neue"),
    legend.text = element_text(size = 10, family = "Helvetica Neue"),
    legend.title = element_text(size = 10, family = "Helvetica Neue"),
    strip.text.x = element_text(size = 10, family = "Helvetica Neue", margin = margin(b = 3)), # pad space below
    panel.background = element_rect(fill = "black"),
    legend.position = "bottom") +
  geom_text(aes(x = Inf, y = Inf, label = rmse),
            size    = 8 / .pt,
            data    = rmse,
            hjust   = 1.05,
            vjust   = 1.5
  )


quartz(width = 8, height = 9)
print(p1)
quartz.save(file = "figures/CMIP6 sea ice comparison v2.jpeg",
            type = "jpeg",
            dev = dev.cur(),
            dpi = 500)
dev.off()












# ends

#cesm2_waccm <- focal(cesm2_waccm, w = matrix(1, nrow = 3, ncol = 3), fun = mean, NAonly = TRUE, na.rm = TRUE) # replace NAs with focal neighborhood mean
#cesm2_waccm <- mask(cesm2_waccm, canesm5) # interpolation leaks onto land so mask by a GCM without the NA issue


