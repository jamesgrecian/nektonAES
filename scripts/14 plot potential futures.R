##############################
### Plot potential futures ###
##############################

### think about how to visualise this in supplementary figures?
### how about a mult-panel showing individual model differences in covariates
### then only the 8 GCM average SSP245 and SSP585 for each species?

# Maps need to look pretty
# could try this https://timogrossenbacher.ch/2016/12/beautiful-thematic-maps-with-ggplot2-only/

# libraries
require(raster)
require(tidyverse)
require(sf)
sf::sf_use_s2(FALSE)
require(patchwork)

# pretty gradient and polar mask functions
source("R/discrete_gradient.R")
source("R/polar_mask.R")
polar_buffer <- polar_mask(radius_size = 5750000)

# source function to format data frame
source("R/futures_df.R")
source("R/plot_futures.R")

# define projection
crs_polar <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# load in shapefile for background mask, clip and project
world_shp <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
CP <- sf::st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 0), crs = 4326) |> sf::st_as_sfc()

world_shp <- world_shp |> sf::st_crop(CP)
world_shp <- world_shp |> sf::st_transform(crs_polar)
world_shp <- world_shp |> st_union()

# load front shapefile
fronts <- orsifronts::parkfronts |> 
  st_as_sf() |> 
  filter(front %in% c("SAF", "PF", "SACCF")) 

# format data
future_tibble <- futures_df(species = 1)  
plot_futures(future_tibble, output_name = "species 01 potential futures.jpg")

future_tibble <- futures_df(species = 2)  
plot_futures(future_tibble, output_name = "species 02 potential futures.jpg")

future_tibble <- futures_df(species = 3)  
plot_futures(future_tibble, output_name = "species 03 potential futures.jpg")

future_tibble <- futures_df(species = 4)  
plot_futures(future_tibble, output_name = "species 04 potential futures.jpg")

future_tibble <- futures_df(species = 5)  
plot_futures(future_tibble, output_name = "species 05 potential futures.jpg")

future_tibble <- futures_df(species = 6)  
plot_futures(future_tibble, output_name = "species 06 potential futures.jpg")

future_tibble <- futures_df(species = 7)  
plot_futures(future_tibble, output_name = "species 07 potential futures.jpg")

future_tibble <- futures_df(species = 8)  
plot_futures(future_tibble, output_name = "species 08 potential futures.jpg")

future_tibble <- futures_df(species = 9)  
plot_futures(future_tibble, output_name = "species 09 potential futures.jpg")

future_tibble <- futures_df(species = 10)  
plot_futures(future_tibble, output_name = "species 10 potential futures.jpg")

future_tibble <- futures_df(species = 11)  
plot_futures(future_tibble, output_name = "species 11 potential futures.jpg")

future_tibble <- futures_df(species = 12)  
plot_futures(future_tibble, output_name = "species 12 potential futures.jpg")

future_tibble <- futures_df(species = 13)  
plot_futures(future_tibble, output_name = "species 13 potential futures.jpg")

future_tibble <- futures_df(species = 14)  
plot_futures(future_tibble, output_name = "species 14 potential futures.jpg")

future_tibble <- futures_df(species = 15)  
plot_futures(future_tibble, output_name = "species 15 potential futures.jpg")

future_tibble <- futures_df(species = 16)  
plot_futures(future_tibble, output_name = "species 16 potential futures.jpg")

future_tibble <- futures_df(species = 17)  
plot_futures(future_tibble, output_name = "species 17 potential futures.jpg")

future_tibble <- futures_df(species = 18)  
plot_futures(future_tibble, output_name = "species 18 potential futures.jpg")

future_tibble <- futures_df(species = 19)  
plot_futures(future_tibble, output_name = "species 19 potential futures.jpg")

future_tibble <- futures_df(species = 20)  
plot_futures(future_tibble, output_name = "species 20 potential futures.jpg")

future_tibble <- futures_df(species = 21)  
plot_futures(future_tibble, output_name = "species 21 potential futures.jpg")

future_tibble <- futures_df(species = 22)  
plot_futures(future_tibble, output_name = "species 22 potential futures.jpg")

future_tibble <- futures_df(species = 23)  
plot_futures(future_tibble, output_name = "species 23 potential futures.jpg")

future_tibble <- futures_df(species = 24)  
plot_futures(future_tibble, output_name = "species 24 potential futures.jpg")

future_tibble <- futures_df(species = 25)  
plot_futures(future_tibble, output_name = "species 25 potential futures.jpg")

future_tibble <- futures_df(species = 26)  
plot_futures(future_tibble, output_name = "species 26 potential futures.jpg")

future_tibble <- futures_df(species = 27)  
plot_futures(future_tibble, output_name = "species 27 potential futures.jpg")

future_tibble <- futures_df(species = 28)  
plot_futures(future_tibble, output_name = "species 28 potential futures.jpg")

# ends
