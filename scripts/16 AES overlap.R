##############################################
### AES overlap with MPAs and predator AES ###
##############################################

# 2024-05-23
# 2025-02-04 - change to max value for all species

# load libraries
require(tidyverse)
require(raster)
require(sf)
sf::sf_use_s2(FALSE)
require(orsifronts)
require(patchwork)
require(scales)
require(wesanderson)

# define projection
prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" 

############################
### Generate AES polygon ###
############################

# load scripts
source("R/aes_poly.R")

# load rasters
fn <- list.files("data/ensemble model outputs", pattern = "ensemble", full.names = T)
fn <- fn[grep("v2", fn)]
x <- lapply(fn, readRDS)
dat <- stack(x)
names(dat) <- c("species01", "species02", "species03", "species04", "species05",
                "species06", "species07", "species08", "species09", "species10",
                "species11", "species12", "species13", "species14", "species15",
                "species16", "species17", "species18", "species19", "species20",
                "species21", "species22", "species23", "species24", "species25",
                "species26", "species27", "species28")

# should we average across all species?
# or should we break into groups?
fish <- subset(dat, c(1:13))
krill <- subset(dat, 14)
squid <- subset(dat, c(15:28))

# average across groups for equal weighting
all_sp <- stack(mean(fish), krill, mean(squid))
all_sp <- max(all_sp)
nekton_raster <- all_sp
nekton_aes <- aes_poly(all_sp, prob = .9)


################################################
### Compare with Hindell et al. predator AES ###
################################################

predator_aes <- readRDS("data/polyW90.rds")
predator_aes <- predator_aes |> st_transform(prj)

# what proportion of nekton aes overlaps with predator aes..?
total_area <- nekton_aes |> st_area() |> sum() |> units::set_units("km^2")
overlap_area <- nekton_aes |> st_intersection(predator_aes) |> st_area() |> sum() |> units::set_units("km^2")
overlap_area/total_area

# how to compare whole are rather than simply the aes?
# spatial comparison of the cell values?
predator_raster <- raster("data/meanImportanceWeighted")
predator_raster <- predator_raster |> projectRaster(all_sp)

test <- stack(nekton_raster, predator_raster)
names(test) <- c("nekton", "predators")
test <- test |> rasterToPoints() |> as_tibble()

# with nclass NULL takes 45 min
mod_test_s <- SpatialPack::modified.ttest(test$nekton,
                                        test$predators,
                                        coords = test[c("x", "y")],
                                        nclass = NULL)
#R = 0.4739, F(1, 40.9734) = 11.8693, p = 0.0013 when using mean
#R = 0.4334, F(1, 57.4789) = 13.2928, p < 0.001 when using max

#########################
### compare with MPAs ###
#########################

# generate source files for comparison
# pretty gradient and polar mask functions
source("R/discrete_gradient.R")
source("R/polar_mask.R")
polar_buffer <- polar_mask(radius_size = 5750000)

# load in shapefile for background mask, clip and project
world_shp <-
  rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
CP <-
  sf::st_bbox(c(
    xmin = -180,
    xmax = 180,
    ymin = -90,
    ymax = 0
  ), crs = 4326) |> sf::st_as_sfc()
world_shp <- world_shp |> sf::st_crop(CP)
world_shp <- world_shp |> sf::st_transform(prj)
world_shp <- world_shp |> st_union()

# MPA files
# https://github.com/ryanreisinger/soPredatorRegions/blob/main/dat_out/mpa_designated.RDS
mpa_designated <- readRDS("data/mpa_designated.rds")
mpa_proposed <- readRDS("data/mpa_proposed.rds")
mpa_designated <- mpa_designated |> st_transform(prj)
mpa_proposed <- mpa_proposed |> st_transform(prj)

# trying to plot these mpa files will crash session
# simplify the vertices
mpa_designated <- mpa_designated |> st_simplify(dTolerance = 1000)
mpa_proposed <- mpa_proposed |> st_simplify(dTolerance = 1000)

# CCAMLR boundary
CCAMLR_bound <- SOmap::SOmap_data$CCAMLR_statistical_areas |> st_as_sf() |> st_union()
CCAMLR_bound <- CCAMLR_bound |> st_transform(prj)
CCAMLR_bound <- CCAMLR_bound |> st_simplify(dTolerance = 5000) # simplify the shape - this will get rid of thick lines

# EEZ boundaries
SOmap_EEZ <- SOmap::SOmap_data$EEZ |> st_as_sf()
SOmap_EEZ <- SOmap_EEZ |> st_transform(prj)
SOmap_EEZ <- SOmap_EEZ |> st_intersection(polar_buffer$buffer) # trim to 40S
SOmap_EEZ <- SOmap_EEZ |> st_union() # fix internal errors that were causing issues with cookie cutting
SOmap_EEZ <- SOmap_EEZ |> st_simplify(dTolerance = 5000) # simplify the shape - this will get rid of thick lines

# ice shelf
ice <- rnaturalearth::ne_download(scale = 50, "antarctic_ice_shelves_polys", category = "physical", returnclass = "sf")
ice <- ice |> st_transform(prj)
ice <- ice |> st_union() # fix internal errors that were causing issues with cookie cutting
  
# some sums

# total area below 40S
polar_buffer$buffer |> st_area() |> sum() |> units::set_units("km^2")

# total area below 40S minus land
polar_buffer$buffer |> st_difference(world_shp) |> st_area() |> sum() |> units::set_units("km^2")

# total area below 40S minus land and ice shelf
land_and_ice <- world_shp |> st_union(ice)
polar_buffer$buffer |> st_difference(land_and_ice) |> st_area() |> sum() |> units::set_units("km^2")

# what is the CCAMLR area excluding land and ice shelf?
CCAMLR_bound |> st_area() |> units::set_units("km^2")
CCAMLR_bound |> st_difference(land_and_ice) |> st_area() |> units::set_units("km^2")

# what is the EEZ area?
SOmap_EEZ |> st_difference(world_shp) |> st_area() |> sum() |> units::set_units("km^2") # EEZ area without land included...
SOmap_EEZ |> st_area() |> sum() |> units::set_units("km^2")

# total area below 40S minus land and ice shelf and EEZ and CCAMLR
# high seas
clip_ccamlr <- polar_mask(radius_size = 3000000)
high_seas <- polar_buffer$buffer |> 
  st_difference(land_and_ice) |> 
  st_difference(SOmap_EEZ) |> 
  st_difference(CCAMLR_bound) |>
  st_difference(clip_ccamlr$buffer)

##################################
### Calculate various overlaps ###
##################################

# high seas
hsa <- high_seas |> st_area() |> sum() |> units::set_units("km^2")
hsa_nekton <- high_seas |> st_intersection(nekton_aes) |> st_area() |> sum() |> units::set_units("km^2")
ggplot() + 
  geom_sf(aes(), data = high_seas) +
  geom_sf(aes(), data = nekton_aes)

# CCAMLR
cmr <- CCAMLR_bound |> st_area() |> sum() |> units::set_units("km^2")
cmr_nekton <- CCAMLR_bound |> st_intersection(nekton_aes) |> st_area() |> sum() |> units::set_units("km^2")
ggplot() + 
  geom_sf(aes(), data = CCAMLR_bound) +
  geom_sf(aes(), data = nekton_aes)

# EEZs
SOmap_EEZ_clip <- SOmap_EEZ |> st_difference(world_shp)
eez <- SOmap_EEZ_clip |> st_area() |> sum() |> units::set_units("km^2")
eez_nekton <- SOmap_EEZ_clip |> st_intersection(nekton_aes) |> st_area() |> sum() |> units::set_units("km^2")

# create dataframe for Juristiction plot
jur_df <- tibble(group = rep(c("inside", "outside"), 3),
                 value = c(eez_nekton, eez-eez_nekton, cmr_nekton, cmr-cmr_nekton, hsa_nekton, hsa-hsa_nekton),
                 area = rep(c("EEZ", "CCAMLR", "high seas"), each = 2)) 
jur_df <- jur_df |> 
  mutate(area = factor(area, levels = c("high seas", "CCAMLR", "EEZ")),
         group = factor(group, levels = c("outside", "inside")))

cmr_nekton/total_area # nekton area in CCAMLR /  total nekton area
eez_nekton/total_area # nekton area in eez / total nekton area
hsa_nekton/total_area # nekton area in high seas / total nekton area

# nekton AES overlap with MPA network
# total area of nekton AES
total_area <- nekton_aes |> st_area() |> sum() |> units::set_units("km^2")
# area of nekton AES that overlaps with designated MPAs
desig_area <- nekton_aes |> st_intersection(mpa_designated) |> st_area() |> sum() |> units::set_units("km^2")
# area of nekton AES that overlaps with proposed MPAs
propo_area <- nekton_aes |> st_intersection(mpa_proposed) |> st_area() |> sum() |> units::set_units("km^2")
# area of nekton AES that is unprotected
unprotected <- total_area - (desig_area + propo_area)

# contrast with proportion of ocean below 40S currently protected
# total area below 40S minus land and ice shelf
land_and_ice <- world_shp |> st_union(ice)
marine_40S <- polar_buffer$buffer |> st_difference(land_and_ice) |> st_area() |> sum() |> units::set_units("km^2")

# what is the size of the MPAs in the shapefile?
mpa_prop <- mpa_proposed |> st_area() |> sum() |> units::set_units("km^2")
mpa_desi <- mpa_designated |> st_area() |> sum() |> units::set_units("km^2")
marine_40S_unprotected <- marine_40S - (mpa_prop + mpa_desi)

# create dataframe for plot
mpa_df <- tibble(group = rep(c("unprotected", "designated", "proposed"), 2),
                 value = c(unprotected, desig_area, propo_area, marine_40S_unprotected, mpa_desi, mpa_prop),
                 area = rep(c("AES", "40 S"), each = 3)) 
mpa_df <- mpa_df |> mutate(group = factor(group, levels = c("unprotected", "proposed", "designated")))

desig_area/total_area
(propo_area + desig_area)/total_area
mpa_desi / marine_40S
(mpa_desi + mpa_prop) / marine_40S


marine_40S
########################
### Compare with CHI ###
########################

# load chi, crop and project
chi <- raster("data/cumulative_impact_2013.tif")
e <- as(extent(-18040095, 18040134, -9020047, 0), 'SpatialPolygons')
chi <- crop(chi, e)
chi <- projectRaster(chi, nekton_raster)
chi <- mask(chi, nekton_raster)

chi_overlap <- stack(nekton_raster, chi)
names(chi_overlap) <- c("nekton", "chi")
chi_overlap <- chi_overlap |> rasterToPoints() |> as_tibble()

 # with or without na's?
# nclass change?
# with nclass NULL takes 45 min
# mod_test_chi <- SpatialPack::modified.ttest(chi_overlap$nekton,
#                                            chi_overlap$chi,
#                                          coords = test[c("x", "y")],
#                                          nclass = NULL)
#R = -0.1779, F(1, 195.381) = 6.3884, p = 0.0123 

# some stats
chi_nekton <- raster::extract(chi, nekton_aes)
chi_nekton <- chi_nekton |> unlist() |> as_tibble()
chi_nekton <- chi_nekton |> mutate(area = "Inside AES")
chi_all <- chi |> rasterToPoints() |> as_tibble() |> pull(cumulative_impact_2013)
chi_all <- tibble(value = chi_all,
                  area = "Southern Ocean below 40S")
chi_ccamlr <-raster::extract(chi, CCAMLR_bound |> st_as_sf())
chi_ccamlr <- chi_ccamlr |> unlist() |> as_tibble()
chi_ccamlr <- chi_ccamlr |> mutate(area = "Inside CCAMLR boundary")
chi_df <- rbind(chi_all, chi_ccamlr, chi_nekton)

####################
### Pretty plots ###
####################

# load fronts
fronts <- parkfronts |> 
  st_as_sf() |> 
  filter(front %in% c("SAF", "PF", "SACCF")) |>
  st_transform(prj)

# nekton and predators #
p1 <- ggplot() + 
  theme_void(base_size = 7,
             base_family = "Helvetica") +
  geom_raster(aes(x, y, fill = MEAN), data = predator_raster |> rasterToPoints()) +
  geom_sf(aes(), colour = "grey60", fill = "grey60",
          data = world_shp) +
  geom_sf(aes(), colour = "black", fill = NA, linewidth = .2, data = predator_aes) +
  geom_sf(aes(), colour = "white", fill = NA, linewidth = .2, nekton_aes) +
  geom_sf(aes(), fill = "white", color = "grey40", linewidth = 0.5 / .pt,
          data = polar_buffer$mask) +
  coord_sf(xlim = c(-6400000, 6400000),
           ylim = c(-6400000, 6400000),
           expand = FALSE,
           crs = projection(krill),
           ndiscr = 1000) +
  xlab(NULL) + 
  ylab(NULL) +
  scale_fill_discrete_gradient("Predator Habitat Importance",
                               colours = viridis::viridis(10),
                               bins = 10,
                               limits = c(0, 100),
                               breaks = seq(0, 100, 20),
                               labels = seq(0, 100, 20),
                               guide = guide_colourbar(
                                 nbin = 500,
                                 raster = T,
                                 frame.colour = "grey40",
                                 ticks.colour = "grey40",
                                 frame.linewidth = .1,
                                 barwidth = 10,
                                 barheight = .5,
                                 direction = "horizontal",
                                 title.position = "top",
                                 #or "right"
                                 title.theme = element_text(
#                                   angle = 90,
                                   hjust = 0.5,
                                   size = 7)
                               )
  ) +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 7,
                                  family = "Helvetica"),
        legend.text = element_text(size = 7, family = "Helvetica"),
        legend.title = element_text(size = 7, family = "Helvetica"),
        panel.background = element_rect(fill = "black"),
        legend.position = "bottom") +
  geom_text(aes(x =  500000,  y = -4500000, label = "Campbell\nPlateau"),      size = 7/.pt, lineheight = .8, family = "Helvetica", angle = 0, colour = "ivory") +
  geom_text(aes(x =  -700000, y = -5000000, label = "Chatham\nRise"),          size = 7/.pt, lineheight = .8, family = "Helvetica", angle = 0, colour = "ivory") +
  geom_text(aes(x =  3500000, y =  1500000, label = "Kerguelen\nPlateau"),     size = 7/.pt, lineheight = .8, family = "Helvetica", angle = 0, colour = "ivory") +
  geom_text(aes(x =  3000000, y =  3500000, label = "Del Cano Rise"),          size = 7/.pt, lineheight = .8, family = "Helvetica", angle = 0, colour = "ivory") +
  geom_text(aes(x =  -700000, y = -1800000, label = "Ross\nSea"),              size = 7/.pt, lineheight = .8, family = "Helvetica", angle = 0, colour = "ivory") +
  geom_text(aes(x = -1000000, y = -3300000, label = "Pacific-Antartic Ridge"), size = 7/.pt, lineheight = .8, family = "Helvetica", angle = 0, colour = "ivory") +
  geom_text(aes(x =   300000, y =  1600000, label = "Dronning\nMaud Land"),    size = 7/.pt, lineheight = .8, family = "Helvetica", angle = 0, colour = "ivory") +
  geom_text(aes(x = -4200000, y =  2400000, label = "Patagonian\nShelf"),      size = 7/.pt, lineheight = .8, family = "Helvetica", angle = 0, colour = "ivory") +
  geom_text(aes(x =  1800000, y =   300000, label = "Prydz\nBay"),             size = 7/.pt, lineheight = .8, family = "Helvetica", angle = 0, colour = "ivory") +
  geom_text(aes(x = -3000000, y =  3500000, label = "South\nGeorgia"),         size = 7/.pt, lineheight = .8, family = "Helvetica", angle = 0, colour = "ivory") +
  geom_text(aes(x = -2700000, y =  2600000, label = "Scotia\nSea"),            size = 7/.pt, lineheight = .8, family = "Helvetica", angle = 0, colour = "ivory") +
  geom_text(aes(x = -2400000, y = -1000000, label = "Amundsen\nSea"),            size = 7/.pt, lineheight = .8, family = "Helvetica", angle = 0, colour = "ivory")


ggsave(filename = "Nekton and Predator AES v5.jpeg",
       plot = p1,
       path = "figures/",
       width = 88,
       height = 100,
       units = "mm",
       dpi = 500)





  
  #annotate("text", label = "INDIAN\nOCEAN", x = 37.00, y = -34.0, size = 4.0, angle = 0, colour = "ivory") +
  #annotate("text", label = "ATLANTIC\nOCEAN", x = 13.10, y = -34.0, size = 4.0, angle = 0, colour = "ivory") +

bathy <- marmap::getNOAA.bathy(lon1 = -180,
                               lon2 = 180,
                               lat1 = -90,
                               lat2 = -40,
                               resolution = 10)
bat <- marmap::as.raster(bathy)
bat <- projectRaster(from = bat, to = predator_raster)

ggplot() + 
  geom_raster(aes(x, y, fill = layer), data = bat |> rasterToPoints()) +
  geom_sf(aes(), colour = "black", fill = NA, linewidth = .2, data = predator_aes) +
  geom_sf(aes(), colour = "white", fill = NA, linewidth = .2, nekton_aes) +
  coord_sf(xlim = c(-6400000, 6400000),
           ylim = c(-6400000, 6400000),
           expand = FALSE,
           crs = projection(bat),
           ndiscr = 1000) +
  geom_text(aes(x = 500000, y = -4500000, label = "Campbell\nPlateau"), size = 2.5, family = "Helvetica Neue", angle = 0, colour = "ivory") +
  geom_text(aes(x = 3500000, y = 1500000, label = "Kerguelen\nPlateau"), size = 2.5, family = "Helvetica Neue", angle = 0, colour = "ivory") +
  geom_text(aes(x = 3700000, y = 3100000, label = "Del Cano Rise"), size = 2.5, family = "Helvetica Neue", angle = 0, colour = "ivory") +
  geom_text(aes(x = -700000, y = -1800000, label = "Ross\nSea"), size = 2.5, family = "Helvetica Neue", angle = 0, colour = "ivory") +
  geom_text(aes(x = -700000, y = -3300000, label = "Pacific-Antartic Ridge"), size = 2.5, family = "Helvetica Neue", angle = 0, colour = "ivory") +
  geom_text(aes(x = 200000, y = 1700000, label = "Dronning Maud\nLand"), size = 2.5, family = "Helvetica Neue", angle = 0, colour = "ivory") +
  geom_text(aes(x = 200000, y = 1700000, label = "South Georgia"), size = 2.5, family = "Helvetica Neue", angle = 0, colour = "ivory")

    

# add some names to the plot?

# eastern South Pacific?
# South Georgia and in the Scotia Sea




# map of nekton, MPAs and jurisdiction
p2 <- ggplot() +
  theme_void(base_size = 7,
             base_family = "Helvetica") +
  geom_raster(aes(x = x, y = y, fill = layer), data = all_sp |> rasterToPoints() |> as_tibble()) +
  geom_sf(aes(colour = status), fill = alpha("#CCBA72", .5), data = mpa_designated) +
  geom_sf(aes(colour = status), fill = alpha("#79402E", .5), data = mpa_proposed) +
  geom_sf(aes(), colour = "grey80", fill = NA, linewidth = .2, data = SOmap_EEZ) +
  geom_sf(aes(), colour = "grey80",fill = NA, linewidth = .2, data = CCAMLR_bound) +
  geom_sf(aes(), colour = "grey60", fill = "grey60", data = world_shp) +
  geom_sf(aes(), colour = "black", fill = "black", data = ice) +
  geom_sf(aes(), colour = "white", fill = NA, linewidth = .2, nekton_aes) +
  geom_sf(aes(), fill = "white", color = "grey40", linewidth = 0.5 / .pt,
          data = polar_buffer$mask) +
  coord_sf(xlim = c(-6400000, 6400000),
           ylim = c(-6400000, 6400000),
           expand = FALSE,
           crs = projection(krill),
           ndiscr = 1000) +
  xlab(NULL) + 
  ylab(NULL) +
  scale_fill_discrete_gradient("Habitat Importance",
                               colours = viridis::viridis(10),
                               bins = 10,
                               limits = c(0, 1),
                               breaks = seq(0, 1, 0.2),
                               labels = seq(0, 1, 0.2),
                               guide = guide_colourbar(
                                 nbin = 500,
                                 raster = T,
                                 frame.colour = "grey40",
                                 ticks.colour = "grey40",
                                 frame.linewidth = .1,
                                 barwidth = 10,
                                 barheight = .5,
                                 direction = "horizontal",
                                 title.position = "top",
                                 #or "right"
                                 title.theme = element_text(
#                                   angle = 90,
                                   hjust = 0.5,
                                   size = 7)
                               )) +
  theme(legend.position = "bottom",
        legend.box="vertical",
        plot.title = element_text(hjust = 0.5,
                                  size = 7,
                                  family = "Helvetica"),
        legend.text = element_text(size = 7, family = "Helvetica"),
        legend.title = element_text(size = 7, family = "Helvetica"),
        panel.background = element_rect(fill = "black"),
        legend.key.size = unit(.4, 'cm'),
        legend.box.spacing = unit(0, "pt")) +
  scale_colour_manual(values = c("#CCBA72", "#79402E"),
                      labels = c("Designated MPA", "Proposed MPA"),
                      guide = guide_legend(title = NULL,
                                           override.aes = list(alpha = 1)))

# plot for MPA overlap
p3 <- ggplot() + 
  theme_bw(base_size = 7,
           base_family = "Helvetica") +
  geom_col(aes(x = area, y = as.numeric(value/1000000), fill = group), position="stack", data = mpa_df) +
  coord_flip() +
  ylab(expression(Area~(million~km^2))) + xlab(NULL) +
  scale_fill_manual(name = NULL, values = wes_palette("IsleofDogs1")) +
  geom_text(aes(x = "AES",  y = 15, label = "Current MPAs: 35.3%"),            vjust = -1, hjust = 0, size = 6/.pt, family = "Helvetica") +
  geom_text(aes(x = "AES",  y = 15, label = "Current & Proposed MPAs: 50.7%"), vjust =  2, hjust = 0, size = 6/.pt, family = "Helvetica") +
  geom_text(aes(x = "40 S", y = 15, label = "Current MPAs: 11.6%"),            vjust = -1, hjust = 0, size = 6/.pt, family = "Helvetica") +
  geom_text(aes(x = "40 S", y = 15, label = "Current & Proposed MPAs: 15.9%"), vjust =  2, hjust = 0, size = 6/.pt, family = "Helvetica") +
  scale_x_discrete(labels = c("> 40°S", "AES")) +
  theme(legend.key.size = unit(.4, 'cm'),
        legend.justification = "left",
        legend.text = element_text(size = 7, family = "Helvetica"),
        legend.title = element_text(size = 7, family = "Helvetica"),
        axis.text = element_text(size = 7, family = "Helvetica"),
        axis.title = element_text(size = 7, family = "Helvetica"))


# plot for jurisdictions
p4 <- ggplot() + 
  theme_bw(base_size = 7,
           base_family = "Helvetica") +
  geom_col(aes(x = area, y = as.numeric(value/1000000), fill = group), position="stack", data = jur_df) +
  coord_flip() +
  ylab(expression(Area~(million~km^2))) + xlab(NULL) +
  scale_fill_manual(name = NULL, values = wes_palette("FrenchDispatch")) +
  geom_text(aes(x = "EEZ",       y = as.numeric(jur_df$value[1])/1000000, label = "15.8% of AES"), hjust = -.05, size = 6/.pt, family = "Helvetica") +
  geom_text(aes(x = "CCAMLR",    y = as.numeric(jur_df$value[3])/1000000, label = "96.2% of AES"), hjust = -.05, size = 6/.pt, family = "Helvetica") +
  geom_text(aes(x = "high seas", y = as.numeric(jur_df$value[5])/1000000, label = "3.8% of AES"),  hjust = -.05, size = 6/.pt, family = "Helvetica") +
  scale_x_discrete(labels = c("High Seas", "CCAMLR", "ANJ")) +
  theme(legend.key.size = unit(.4, 'cm'),
        legend.justification = "left",
        legend.text = element_text(size = 7, family = "Helvetica"),
        legend.title = element_text(size = 7, family = "Helvetica"),
        axis.text = element_text(size = 7, family = "Helvetica"),
        axis.title = element_text(size = 7, family = "Helvetica"))


#p_out <- p2 + (plot_spacer()/p3/p4/plot_spacer() + plot_layout(heights = c(1, 3, 3, 1))) + plot_layout(widths = c(2, 1))

p_out <- p2 + (plot_spacer()/p3/p4/plot_spacer() + plot_layout(heights = c(1, 4, 4, 1))) + 
  plot_layout(widths = c(2, 1)) +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 9,
                                                                    family = "Helvetica Neue",
                                                                    face = "bold"))

ggsave(filename = "Nekton MPA and Jurisdiction v4.jpeg",
       plot = p_out,
       path = "figures/",
       width = 180,
       height = 120,
       units = "mm",
       dpi = 500)

# nekton and chi #
p5 <- ggplot() + 
  theme_void(base_size = 8,
             base_family = "Helvetica Neue") +
  geom_raster(aes(x, y, fill = cumulative_impact_2013), data = chi |> rasterToPoints()) +
#  geom_sf(aes(), colour = "light grey", linewidth = .1,
#          data = fronts) + 
  geom_sf(aes(), fill = NA, linewidth = .2, data = CCAMLR_bound) +
#  geom_text(aes(x = 2700000, y = -2100000, label = "SACCF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
#  geom_text(aes(x = 3050000, y = -2550000, label = "PF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
#  geom_text(aes(x = 3450000, y = -2950000, label = "SAF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_sf(aes(), colour = "grey60", fill = "grey60",
          data = world_shp) +
  geom_sf(aes(), colour = "black", fill = NA, linewidth = .2, nekton_aes) +
  geom_sf(aes(), fill = "white", color = "grey40", linewidth = 0.5 / .pt,
          data = polar_buffer$mask) +
  coord_sf(xlim = c(-6400000, 6400000),
           ylim = c(-6400000, 6400000),
           expand = FALSE,
           crs = projection(krill),
           ndiscr = 1000) +
  xlab(NULL) + 
  ylab(NULL) +
  scale_fill_discrete_gradient("Global cumulative human impact",
                               colours = rev(pals::brewer.rdylbu(n = 18)),
                               bins = 18,
                               transform = "sqrt",
                               limits = c(0, 4),
                               breaks = c(0, 1, 2.08, 3.16, 4),
                               labels = seq(0, 4),
                               oob = scales::squish,
                               guide = guide_colourbar(
                                 nbin = 500,
                                 frame.colour = "grey40",
                                 ticks.colour = "grey40",
                                 frame.linewidth = .1,
                                 barwidth = 10,
                                 barheight = .5,
                                 direction = "horizontal",
                                 title.position = "top",
                                 #or "right"
                                 title.theme = element_text(
                                   #angle = 90,
                                   hjust = 0.5,
                                   size = 8)
                               )) +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 8,
                                  family = "Helvetica Neue"),
        legend.text = element_text(size = 8, family = "Helvetica Neue"),
        legend.title = element_text(size = 8, family = "Helvetica Neue"),
        panel.background = element_rect(fill = "black"),
        legend.position = "bottom")

p6 <- ggplot() +
  theme_bw(base_size = 8,
           base_family = "Helvetica Neue") +
  geom_density(aes(x = value, y = ..scaled.., fill = area), data = chi_df, outline.type = "full") +
  scale_fill_manual(name = NULL, values = wes_palette("Chevalier1")) +
  scale_x_continuous("Cumulative Human Impact", limits = c(0, 4)) +
  scale_y_continuous("Scaled density", limits = c(0, 1)) +
  facet_wrap(~area, nrow = 3) +
  theme(strip.background = element_rect(fill = NA, size = 0),
        strip.text = element_text(size = 8),
        legend.position = "none")

ggsave(filename = "Nekton and CHI v2.jpeg",
       plot = p5,
       path = "figures/",
       width = 90,
       height = 120,
       units = "mm",
       dpi = 500)





# ends










# other ideas...

# require(mgcv)
# m1 <- gam(predators ~ nekton + te(x,y), data = test)
# summary(m1)
# m2 <- glm(predators ~ 1, data = test)
# anova(m1, m2, test = "F")
# summary(m1)

# in the biological conservation paper
# bootstrap R 10,000 times by subsampling 10% data
# cor_test <- NA
# for (i in 1:10000){
#  samp <- test |> sample_n(nrow(test)*.1)
#  cor_test[i] <- cor(samp$nekton, samp$predators, method = "spearman")
# }
# summary(cor_test)

# Ewan suggested Bhattacharyya's affinity...?
# Bhattacharyya's affinity manually
# norm_nekton <- nekton_raster/cellStats(nekton_raster, sum)
# norm_predator <- predator_raster/cellStats(predator_raster, sum)
# sum(sqrt(prob1 * prob2))
# cellStats(sqrt(norm_nekton) * sqrt(norm_predator), sum)
# sum(((raster1^0.5)*(raster2^0.5)))
# cellStats((norm_nekton^.5) * (norm_predator^.5), stat = "sum", na.rm = T)
# sum(sqrt(ud_i) * sqrt(ud_j))

# cor.test(test$nekton, test$predators)

# correlation_func <- function(data, indices) {
#  # Select data based on bootstrapped indices
#  boot_data <- data[indices, ]
#  # Calculate correlation coefficient
#  cor(boot_data[, 1], boot_data[, 2])
# }
# require(boot)
# boot_results <- boot(test[,3:4], correlation_func, R = 1000, parallel = "multicore") 

# require(terra)
# r.cor <- rasterCorrelation(nekton_raster |> rast(),
#                           predator_raster |> rast(),
#                           s = 5,
#                           type = "pearson")
# plot(r.cor)
