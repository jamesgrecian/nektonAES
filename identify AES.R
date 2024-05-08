###################################################
### Identify nekton AES from predicted surfaces ###
###################################################

# 2024-05-01

# load individual species rasters and combine to estimate nekton AES

# load libraries
require(tidyverse)
require(raster)
require(sf)
sf::sf_use_s2(FALSE)
require(orsifronts)
require(patchwork)

source("R/aes_poly.R")

# load rasters
fn <- list.files("data/ensemble model outputs", pattern = "ensemble", full.names = T)
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
all_sp <- mean(all_sp)

fish_aes <- aes_poly(mean(fish), prob = .9)
krill_aes <- aes_poly(krill, prob = .9)
squid_aes <- aes_poly(mean(squid), prob = .9)
nekton_aes <- aes_poly(all_sp, prob = .9)

# define projection
prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" 

# load fronts
fronts <- parkfronts |> 
  st_as_sf() |> 
  filter(front %in% c("SAF", "PF", "SACCF")) |>
  st_transform(prj)

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

##################################
### Generate some pretty plots ###
##################################

# fish plot
p1 <- ggplot() +
  theme_void(base_size = 10,
             base_family = "Helvetica Neue") +
  ggtitle("Fish") +
  geom_raster(aes(x = x, y = y, fill = layer), data = mean(fish) |> rasterToPoints() |> as_tibble()) +
  geom_sf(aes(), colour = "light grey", linewidth = .1,
          data = fronts) + 
  geom_text(aes(x = 2700000, y = -2100000, label = "SACCF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_text(aes(x = 3050000, y = -2550000, label = "PF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_text(aes(x = 3450000, y = -2950000, label = "SAF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_sf(aes(), colour = "grey60", fill = "grey60",
          data = world_shp) +
  geom_sf(aes(), colour = "white", fill = NA, linewidth = .5, fish_aes) +
  geom_sf(aes(), fill = "white", color = "grey40", size = 0.5 / .pt,
          data = polar_buffer$mask) +
  coord_sf(xlim = c(-6400000, 6400000),
           ylim = c(-6400000, 6400000),
           expand = FALSE,
           crs = prj,
           ndiscr = 1000) +
  xlab(NULL) + 
  ylab(NULL) +
  scale_fill_discrete_gradient("Habitat Suitability",
                               colours = viridis::viridis(10),
                               bins = 10,
                               limits = c(0, 1),
                               breaks = seq(0, 1, 0.2),
                               labels = seq(0, 1, 0.2),
                               oob = squish,
                               guide = guide_colourbar(
                                 nbin = 500,
                                 frame.colour = "grey40",
                                 ticks.colour = "grey40",
                                 frame.linewidth = .1,
                                 barwidth = .5,
                                 barheight = 10,
                                 direction = "vertical",
                                 title.position = "right",
                                 #or "right"
                                 title.theme = element_text(
                                   angle = 90,
                                   hjust = 0.5,
                                   size = 10)
                               )
  ) +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10,
                                  family = "Helvetica Neue"),
        legend.text = element_text(size = 10, family = "Helvetica Neue"),
        legend.title = element_text(size = 10, family = "Helvetica Neue"),
        panel.background = element_rect(fill = "black"))

# krill plot
p2 <- ggplot() +
  theme_void(base_size = 10,
             base_family = "Helvetica Neue") +
  ggtitle("Krill") +
  geom_raster(aes(x = x, y = y, fill = species14), data = krill |> rasterToPoints() |> as_tibble()) +
  geom_sf(aes(), colour = "light grey", linewidth = .1,
          data = fronts) + 
  geom_text(aes(x = 2700000, y = -2100000, label = "SACCF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_text(aes(x = 3050000, y = -2550000, label = "PF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_text(aes(x = 3450000, y = -2950000, label = "SAF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_sf(aes(), colour = "grey60", fill = "grey60",
          data = world_shp) +
  geom_sf(aes(), colour = "white", fill = NA, linewidth = .5, krill_aes) +
  geom_sf(aes(), fill = "white", color = "grey40", size = 0.5 / .pt,
          data = polar_buffer$mask) +
  coord_sf(xlim = c(-6400000, 6400000),
           ylim = c(-6400000, 6400000),
           expand = FALSE,
           crs = projection(krill),
           ndiscr = 1000) +
  xlab(NULL) + 
  ylab(NULL) +
  scale_fill_discrete_gradient("Habitat Suitability",
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
                                 barwidth = .5,
                                 barheight = 10,
                                 direction = "vertical",
                                 title.position = "right",
                                 #or "right"
                                 title.theme = element_text(
                                   angle = 90,
                                   hjust = 0.5,
                                   size = 10)
                               )
  ) +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10,
                                  family = "Helvetica Neue"),
        legend.text = element_text(size = 10, family = "Helvetica Neue"),
        legend.title = element_text(size = 10, family = "Helvetica Neue"),
        panel.background = element_rect(fill = "black"))

# squid plot
p3 <- ggplot() +
  theme_void(base_size = 10,
             base_family = "Helvetica Neue") +
  ggtitle("Squid") +
  geom_raster(aes(x = x, y = y, fill = layer), data = mean(squid) |> rasterToPoints() |> as_tibble()) +
  geom_sf(aes(), colour = "light grey", linewidth = .1,
          data = fronts) + 
  geom_text(aes(x = 2700000, y = -2100000, label = "SACCF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_text(aes(x = 3050000, y = -2550000, label = "PF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_text(aes(x = 3450000, y = -2950000, label = "SAF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_sf(aes(), colour = "grey60", fill = "grey60",
          data = world_shp) +
  geom_sf(aes(), colour = "white", fill = NA, linewidth = .5, squid_aes) +
  geom_sf(aes(), fill = "white", color = "grey40", size = 0.5 / .pt,
          data = polar_buffer$mask) +
  coord_sf(xlim = c(-6400000, 6400000),
           ylim = c(-6400000, 6400000),
           expand = FALSE,
           crs = projection(krill),
           ndiscr = 1000) +
  xlab(NULL) + 
  ylab(NULL) +
  scale_fill_discrete_gradient("Habitat Suitability",
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
                                 barwidth = .5,
                                 barheight = 10,
                                 direction = "vertical",
                                 title.position = "right",
                                 #or "right"
                                 title.theme = element_text(
                                   angle = 90,
                                   hjust = 0.5,
                                   size = 10)
                               )
  ) +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10,
                                  family = "Helvetica Neue"),
        legend.text = element_text(size = 10, family = "Helvetica Neue"),
        legend.title = element_text(size = 10, family = "Helvetica Neue"),
        panel.background = element_rect(fill = "black"))

# all nekton species
p4 <- ggplot() +
  theme_void(base_size = 10,
             base_family = "Helvetica Neue") +
  ggtitle("All Species") +
  geom_raster(aes(x = x, y = y, fill = layer), data = all_sp |> rasterToPoints() |> as_tibble()) +
  geom_sf(aes(), colour = "light grey", linewidth = .1,
          data = fronts) + 
  geom_text(aes(x = 2700000, y = -2100000, label = "SACCF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_text(aes(x = 3050000, y = -2550000, label = "PF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_text(aes(x = 3450000, y = -2950000, label = "SAF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_sf(aes(), colour = "grey60", fill = "grey60",
          data = world_shp) +
  geom_sf(aes(), colour = "white", fill = NA, linewidth = .5, nekton_aes) +
  geom_sf(aes(), fill = "white", color = "grey40", size = 0.5 / .pt,
          data = polar_buffer$mask) +
  coord_sf(xlim = c(-6400000, 6400000),
           ylim = c(-6400000, 6400000),
           expand = FALSE,
           crs = projection(krill),
           ndiscr = 1000) +
  xlab(NULL) + 
  ylab(NULL) +
  scale_fill_discrete_gradient("Habitat Suitability",
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
                                 barwidth = .5,
                                 barheight = 10,
                                 direction = "vertical",
                                 title.position = "right",
                                 #or "right"
                                 title.theme = element_text(
                                   angle = 90,
                                   hjust = 0.5,
                                   size = 10)
                               )
  ) +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10,
                                  family = "Helvetica Neue"),
        legend.text = element_text(size = 10, family = "Helvetica Neue"),
        legend.title = element_text(size = 10, family = "Helvetica Neue"),
        panel.background = element_rect(fill = "black"))

# stick plots together with patchwork
# save to output path
quartz(width = 9, height = 8)
p4 + (p1 / p2 / p3) + plot_layout(guides = "collect",
                                  widths = c(3, 1))
quartz.save(
  file = "figures/Nekton AES.jpeg",
  type = "jpeg",
  dev = dev.cur(),
  dpi = 500
)
dev.off()

#########################
### compare with MPAs ###
#########################

# https://github.com/ryanreisinger/soPredatorRegions/blob/main/dat_out/mpa_designated.RDS
mpa_designated <- readRDS("~/Downloads/mpa_designated.rds")
mpa_proposed <- readRDS("~/Downloads/mpa_proposed.rds")
mpa_designated <- mpa_designated |> st_transform(prj)
mpa_proposed <- mpa_proposed |> st_transform(prj)

# all nekton species
p5 <- ggplot() +
  theme_void(base_size = 10,
             base_family = "Helvetica Neue") +
  ggtitle("All Species") +
  geom_raster(aes(x = x, y = y, fill = layer), data = all_sp |> rasterToPoints() |> as_tibble()) +
  geom_sf(aes(), colour = "light grey", linewidth = .1,
          data = fronts) + 
  geom_text(aes(x = 2700000, y = -2100000, label = "SACCF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_text(aes(x = 3050000, y = -2550000, label = "PF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_text(aes(x = 3450000, y = -2950000, label = "SAF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_sf(aes(colour = status), fill = alpha("#EDCE76", .5), data = mpa_designated) +
  geom_sf(aes(colour = status), fill = alpha("#CA4F4E", .5), data = mpa_proposed) +
  geom_sf(aes(), colour = "grey60", fill = "grey60",
          data = world_shp) +
  geom_sf(aes(), colour = "white", fill = NA, linewidth = .5, nekton_aes) +
  geom_sf(aes(), fill = "white", color = "grey40", size = 0.5 / .pt,
          data = polar_buffer$mask) +
  coord_sf(xlim = c(-6400000, 6400000),
           ylim = c(-6400000, 6400000),
           expand = FALSE,
           crs = projection(krill),
           ndiscr = 1000) +
  xlab(NULL) + 
  ylab(NULL) +
  scale_fill_discrete_gradient("Habitat Suitability",
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
                                 barwidth = .5,
                                 barheight = 10,
                                 direction = "vertical",
                                 title.position = "right",
                                 #or "right"
                                 title.theme = element_text(
                                   angle = 90,
                                   hjust = 0.5,
                                   size = 10)
                               )
  ) +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10,
                                  family = "Helvetica Neue"),
        legend.text = element_text(size = 10, family = "Helvetica Neue"),
        legend.title = element_text(size = 10, family = "Helvetica Neue"),
        panel.background = element_rect(fill = "black")) +
  scale_colour_manual(values = c("#EDCE76", "#CA4F4E"),
                      labels = c("Designated MPA", "Proposed MPA"),
                      guide = guide_legend(title = NULL,
                                           override.aes = list(alpha = 1)))



quartz(width = 8, height = 8)
p5
quartz.save(
  file = "figures/Nekton MPA.jpeg",
  type = "jpeg",
  dev = dev.cur(),
  dpi = 500
)
dev.off()

# ends



# overlap?
total_area <- nekton_aes |> st_area() |> sum() |> units::set_units("km^2")
desig_area <- nekton_aes |> st_intersection(mpa_designated) |> st_area() |> sum() |> units::set_units("km^2")
propo_area <- nekton_aes |> st_intersection(mpa_proposed) |> st_area() |> sum() |> units::set_units("km^2")
unprotected <- total_area - (desig_area + propo_area)
foo <- tibble(group = c("unprotected", "designated", "proposed"),
              value = c(unprotected, desig_area, propo_area),
              area = rep("AES", 3))


CCAMLR_bound <- SOmap::SOmap_data$CCAMLR_statistical_areas |> st_as_sf() |> st_union()
CCAMLR_bound <- CCAMLR_bound |> st_transform(prj)
CCAMLR_bound |> st_area() |> units::set_units("km^2")
nekton_aes |> st_intersection(CCAMLR_bound) |> st_area() |> units::set_units("km^2")


ggplot() + 
  theme_bw() +
  geom_bar(aes(x = area, y = as.numeric(value/1000000), fill = group), position="stack", stat = "identity", data = foo) +
  coord_flip() +
  ylab("Area (million km2)") + xlab(NULL)



ggplot() + 
  geom_bar(aes(fill=condition, y=value, x=specie), position="stack", stat="identity", foo)


specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)
