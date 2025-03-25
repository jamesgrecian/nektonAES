###################################################
### Identify nekton AES from predicted surfaces ###
###################################################

# 2024-05-01
# 2025-01-31 average across species within the guild, but then take max value across guilds

# load individual species rasters and combine to estimate nekton AES

# load libraries
require(tidyverse)
require(raster)
require(sf)
sf::sf_use_s2(FALSE)
require(orsifronts)
require(patchwork)
require(scales)

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
# but take maximum value from each guild to retain high importance areas
all_sp <- max(all_sp)

# output guild aes for figure
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

# CCAMLR boundary
CCAMLR_bound <- SOmap::SOmap_data$CCAMLR_statistical_areas |> st_as_sf() |> st_union()
CCAMLR_bound <- CCAMLR_bound |> st_transform(prj)
CCAMLR_bound <- CCAMLR_bound |> st_simplify(dTolerance = 5000) # simplify the shape - this will get rid of thick lines

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
  theme_void(base_size = 7,
             base_family = "Helvetica") +
  ggtitle("Fish") +
  geom_raster(aes(x = x, y = y, fill = layer), data = mean(fish) |> rasterToPoints() |> as_tibble()) +
  geom_sf(aes(), colour = "light grey", linewidth = .1,
          data = fronts) + 
  geom_sf(aes(), colour = "grey60", fill = "grey60",
          data = world_shp) +
  geom_sf(aes(), colour = "white", fill = NA, linewidth = .2, nekton_aes) +
  geom_sf(aes(), fill = "white", color = "white",
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
                                   size = 7)
                               )
  ) +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 7,
                                  family = "Helvetica",
                                  colour = "black"),
        legend.text = element_text(size = 7, family = "Helvetica"),
        legend.title = element_text(size = 7, family = "Helvetica"),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.background = element_rect(fill = "black", colour = "black"),
        legend.position = "none")

# krill plot
p2 <- ggplot() +
  theme_void(base_size = 7,
             base_family = "Helvetica") +
  ggtitle("Krill") +
  geom_raster(aes(x = x, y = y, fill = species14), data = krill |> rasterToPoints() |> as_tibble()) +
  geom_sf(aes(), colour = "light grey", linewidth = .1,
          data = fronts) + 
  geom_sf(aes(), colour = "grey60", fill = "grey60",
          data = world_shp) +
  geom_sf(aes(), colour = "white", fill = NA, linewidth = .2, nekton_aes) +
  geom_sf(aes(), fill = "white", color = "white",
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
                                   size = 7)
                               )
  ) +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 7,
                                  family = "Helvetica",
                                  colour = "black"),
        legend.text = element_text(size = 7, family = "Helvetica"),
        legend.title = element_text(size = 7, family = "Helvetica"),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.background = element_rect(fill = "black", colour = "black"),
        legend.position = "none")

# squid plot
p3 <- ggplot() +
  theme_void(base_size = 7,
             base_family = "Helvetica") +
  ggtitle("Squid") +
  geom_raster(aes(x = x, y = y, fill = layer), data = mean(squid) |> rasterToPoints() |> as_tibble()) +
  geom_sf(aes(), colour = "light grey", linewidth = .1,
          data = fronts) + 
  geom_sf(aes(), colour = "grey60", fill = "grey60",
          data = world_shp) +
  geom_sf(aes(), colour = "white", fill = NA, linewidth = .2, nekton_aes) +
  geom_sf(aes(), fill = "white", color = "white",
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
                                   size = 7)
                               )
  ) +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 7,
                                  family = "Helvetica",
                                  colour = "black"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        legend.title = element_text(size = 7, family = "Helvetica"),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.background = element_rect(fill = "black", colour = "black"),
        legend.position = "none")

# all nekton species
p4 <- ggplot() +
  theme_void(base_size = 7,
             base_family = "Helvetica") +
  ggtitle("All Species") +
  geom_raster(aes(x = x, y = y, fill = layer), data = all_sp |> rasterToPoints() |> as_tibble()) +
  geom_sf(aes(), colour = "light grey", linewidth = .1,
          data = fronts) + 
  geom_text(aes(x = 2500000, y = 2650000, label = "SACCF"), size = 7/.pt , family = "Helvetica", angle = -45, colour = "grey") +
  geom_text(aes(x = 2950000, y = 3150000, label = "PF"), size = 7/.pt, family = "Helvetica", angle = -45, colour = "grey") +
  geom_text(aes(x = 3250000, y = 3500000, label = "SAF"), size = 7/.pt, family = "Helvetica", angle = -45, colour = "grey") +
#  geom_text(aes(x = 2700000, y = -2100000, label = "SACCF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
#  geom_text(aes(x = 3050000, y = -2550000, label = "PF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
#  geom_text(aes(x = 3450000, y = -2950000, label = "SAF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
  geom_sf(aes(), colour = "grey80",fill = NA, linewidth = .2, data = CCAMLR_bound) +
  geom_sf(aes(), colour = "grey60", fill = "grey60",
          data = world_shp) +
  geom_sf(aes(), colour = "white", fill = NA, linewidth = .3, nekton_aes) +
  geom_sf(aes(), fill = "white", color = "white",
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
                                   hjust = 0.5,
                                   size = 7)
                               )
  ) +
  theme(plot.title = element_text(hjust = 0.5,
                                  vjust = -15,
                                  size = 7,
                                  family = "Helvetica",
                                  colour = "black"),
        legend.text = element_text(size = 7, family = "Helvetica", colour = "black"),
        legend.title = element_text(size = 7, family = "Helvetica", colour = "black"),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.background = element_rect(fill = "black", colour = "black"),
        legend.box.spacing = unit(0, "pt"),
        legend.box.margin = margin(t = -40, r = 0, b= 0, l = 0, unit = "pt"),
        legend.position = "bottom")

# stick plots together with patchwork
# save to output path
p <- p4 + (p1 / p2 / p3) + 
  plot_layout(widths = c(2, 1)) +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 9,
                                                                    family = "Helvetica",
                                                                    face = "bold",
                                                                    colour = "black"),
                                            plot.background = element_rect(colour = NA, fill = NA),
                                            panel.background = element_rect(colour = NA, fill = "black"))

ggsave(filename = "Nekton AES v6 white.jpeg",
       plot = p,
       path = "figures/",
       device = "jpeg",
       width = 160,
       height = 140,
       units = "mm")


