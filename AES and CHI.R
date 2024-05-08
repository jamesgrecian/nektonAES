
require(raster)
require(tidyverse)
require(sf)
sf::sf_use_s2(FALSE)

source("~/nektonAES/R/polar_mask.R")
polar_buffer <- polar_mask(radius_size = 5750000)

source("R/discrete_gradient.R")


foo <- raster("data/cumulative_impact_2013.tif")

prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

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

e <- as(extent(-18040095, 18040134, -9020047, 0), 'SpatialPolygons') 
r <- crop(foo, e)
foo2 <- projectRaster(r, res = 25000, crs = prj)
plot(foo2)

foo3 <- foo2 |> mask(as_Spatial(polar_buffer$buffer))
plot(foo3)

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


CCAMLR_bound <- SOmap::SOmap_data$CCAMLR_statistical_areas |> st_as_sf() |> st_union()

p1 <- ggplot() + 
  theme_void(base_size = 10,
             base_family = "Helvetica Neue") +
  geom_raster(aes(x = x, y = y, fill = cumulative_impact_2013), data = foo3 |> rasterToPoints() |> as_tibble()) +
  geom_sf(aes(), colour = "black", fill = NA, linewidth = .5, nekton_aes) +
  geom_sf(aes(), fill = NA, colour = "grey60", linewidth = .5, data = CCAMLR_bound) +
  geom_sf(aes(), colour = "grey60", fill = "grey60",
          data = world_shp) +
  geom_sf(aes(), fill = "white", color = "grey40", size = 0.5 / .pt,
          data = polar_buffer$mask) +
  coord_sf(xlim = c(-6400000, 6400000),
           ylim = c(-6400000, 6400000),
           expand = FALSE,
           crs = prj,
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
                                 barwidth = 15,
                                 barheight = 1,
                                 direction = "horizontal",
                                 title.position = "top",
                                 #or "right"
                                 title.theme = element_text(
                                   #angle = 90,
                                   hjust = 0.5,
                                   size = 10)
                               )) +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10,
                                  family = "Helvetica Neue"),
        legend.text = element_text(size = 10, family = "Helvetica Neue"),
        legend.title = element_text(size = 10, family = "Helvetica Neue"),
        panel.background = element_rect(fill = "black"),
        legend.position = "bottom")


# save to output path
quartz(width = 9, height = 8)
p1
quartz.save(
  file = "Figures/Nekton AES and CHI v2.jpeg",
  type = "jpeg",
  dev = dev.cur(),
  dpi = 500
)
dev.off()


# some stats
chi_nekton <-raster::extract(foo3, nekton_aes)
chi_nekton <- chi_nekton |> unlist() |> as_tibble()
chi_nekton <- chi_nekton |> mutate(area = "Inside AES")
chi_all <- foo3 |> rasterToPoints() |> as_tibble() |> pull(cumulative_impact_2013)
chi_all <- tibble(value = chi_all,
                  area = "Southern Ocean below 30S")
chi_ccamlr <-raster::extract(foo3, SOmap::SOmap_data$CCAMLR_statistical_areas)
chi_ccamlr <- chi_ccamlr |> unlist() |> as_tibble()
chi_ccamlr <- chi_ccamlr |> mutate(area = "Inside CCAMLR boundary")
chi_df <- rbind(chi_all, chi_ccamlr, chi_nekton)

p2 <- ggplot() +
  theme_bw(base_size = 10,
           base_family = "Helvetica Neue") +
  geom_density(aes(x = value, y = ..scaled.., fill = area), data = chi_df, outline.type = "full") +
  scale_x_continuous("Cumulative Human Impact", limits = c(0, 4)) +
  scale_y_continuous("Scaled density", limits = c(0, 1)) +
  facet_wrap(~area, nrow = 3) +
  theme(strip.background = element_rect(fill = NA, size = 0),
        strip.text = element_text(size = 10),
        legend.position = "none")


require(patchwork)

quartz(width = 9, height = 7)
p1 + (plot_spacer() / p2 / plot_spacer() + plot_layout(heights = c(1, 2, 1))) + plot_layout(widths = c(3, 1))
quartz.save(
  file = "Figures/Nekton AES and CHI v3.jpeg",
  type = "jpeg",
  dev = dev.cur(),
  dpi = 500
)
dev.off()

