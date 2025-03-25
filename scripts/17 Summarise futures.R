###################################################################
### Calculate, summarise and plot patterns in potential futures ###
###################################################################

# 2024-09-23
# 2025-02-05 - change from mean across guilds to max across guilds

# load in all potential futures
# calculate AES for each GCM
# then explore shifts...

# output should be two stacks (or similar) with 8 layers each showing change in AES

# libraries
require(raster)
require(tidyverse)
require(sf)
sf::sf_use_s2(FALSE)
require(patchwork)
source("R/aes_poly.R")

# pretty gradient and polar mask functions
source("R/discrete_gradient.R")
source("R/polar_mask.R")
polar_buffer <- polar_mask(radius_size = 5750000)

# define projection
crs_polar <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# load in shapefile for background mask, clip and project
world_shp <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
CP <- sf::st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 0), crs = 4326) |> sf::st_as_sfc()

world_shp <- world_shp |> sf::st_crop(CP)
world_shp <- world_shp |> sf::st_transform(crs_polar)
world_shp <- world_shp |> st_union()

### Process contemporary data ###
files_now <- list.files("data/ensemble model outputs", pattern = "ensemble", full.names = T)
files_now <- files_now[grep("v2", files_now)]
stack_now <- lapply(files_now, readRDS) |> stack()
names(stack_now) <- c("species01", "species02", "species03", "species04", "species05",
                      "species06", "species07", "species08", "species09", "species10",
                      "species11", "species12", "species13", "species14", "species15",
                      "species16", "species17", "species18", "species19", "species20",
                      "species21", "species22", "species23", "species24", "species25",
                      "species26", "species27", "species28")
guilds <- c(rep(1, times = 13), 2, rep(3, times = 14))
stack_now_guilds <- stackApply(stack_now, guilds, mean)
#stack_now_mean <- mean(stack_now_guilds)
stack_now_max <- max(stack_now_guilds) # take max value across the guilds

### Process SSP2-4.5 data ###
files_245 <- list.files("data/ensemble model outputs", pattern = "ssp245", full.names = T)
stack_245 <- lapply(files_245, readRDS) 
stack_245_1 <- lapply(stack_245, subset, "ACCESS.CM2") |> stack() |> stackApply(guilds, mean)
stack_245_2 <- lapply(stack_245, subset, "BCC.CSM2.MR") |> stack() |> stackApply(guilds, mean)
stack_245_3 <- lapply(stack_245, subset, "CESM2.WACCM") |> stack() |> stackApply(guilds, mean)
stack_245_4 <- lapply(stack_245, subset, "CMCC.CM2.SR5") |> stack() |> stackApply(guilds, mean)
stack_245_5 <- lapply(stack_245, subset, "FGOALS.g3") |> stack() |> stackApply(guilds, mean)
stack_245_6 <- lapply(stack_245, subset, "IPSL.CM6A.LR") |> stack() |> stackApply(guilds, mean)
stack_245_7 <- lapply(stack_245, subset, "MRI.ESM2.0") |> stack() |> stackApply(guilds, mean)
stack_245_8 <- lapply(stack_245, subset, "NESM3") |> stack() |> stackApply(guilds, mean)

# this is now the max across the guilds for each ssp
aes_245 <- stack(max(stack_245_1), max(stack_245_2), max(stack_245_3), max(stack_245_4),
                 max(stack_245_5), max(stack_245_6), max(stack_245_7), max(stack_245_8))
names(aes_245) <- c("ACCESS.CM2", "BCC.CSM2.MR", "CESM2.WACCM", "CMCC.CM2.SR5",
                    "FGOALS.g3", "IPSL.CM6A.LR", "MRI.ESM2.0", "NESM3")

### Process SSP5-8.5 data ###
files_585 <- list.files("data/ensemble model outputs", pattern = "ssp585", full.names = T)
stack_585 <- lapply(files_585, readRDS) 
stack_585_1 <- lapply(stack_585, subset, "ACCESS.CM2") |> stack() |> stackApply(guilds, mean)
stack_585_2 <- lapply(stack_585, subset, "BCC.CSM2.MR") |> stack() |> stackApply(guilds, mean)
stack_585_3 <- lapply(stack_585, subset, "CESM2.WACCM") |> stack() |> stackApply(guilds, mean)
stack_585_4 <- lapply(stack_585, subset, "CMCC.CM2.SR5") |> stack() |> stackApply(guilds, mean)
stack_585_5 <- lapply(stack_585, subset, "FGOALS.g3") |> stack() |> stackApply(guilds, mean)
stack_585_6 <- lapply(stack_585, subset, "IPSL.CM6A.LR") |> stack() |> stackApply(guilds, mean)
stack_585_7 <- lapply(stack_585, subset, "MRI.ESM2.0") |> stack() |> stackApply(guilds, mean)
stack_585_8 <- lapply(stack_585, subset, "NESM3") |> stack() |> stackApply(guilds, mean)

# this is now the max across the guilds for each ssp
aes_585 <- stack(max(stack_585_1), max(stack_585_2), max(stack_585_3), max(stack_585_4),
                 max(stack_585_5), max(stack_585_6), max(stack_585_7), max(stack_585_8))
names(aes_585) <- c("ACCESS.CM2", "BCC.CSM2.MR", "CESM2.WACCM", "CMCC.CM2.SR5",
                    "FGOALS.g3", "IPSL.CM6A.LR", "MRI.ESM2.0", "NESM3")

### Does AES area change? ###
# don't use prob = .9 - calculate the value of the threshold for contemporary and then use that
# how does current threshold then become unsuitable?

contemporary_aes <- aes_poly(stack_now_max, prob = .9)

quant <- quantile(stack_now_max, p = .9) # now how to pass that to the GCM outputs

current_aes_size <- contemporary_aes |> st_area() |> sum() |> units::set_units("km^2")

future_aes_size <- tibble(GCM = names(aes_245),
                          SSP245 = c(aes_poly(subset(aes_245, 1), q = quant) |> st_area() |> sum() |> units::set_units("km^2"),
                                     aes_poly(subset(aes_245, 2), q = quant) |> st_area() |> sum() |> units::set_units("km^2"),
                                     aes_poly(subset(aes_245, 3), q = quant) |> st_area() |> sum() |> units::set_units("km^2"),
                                     aes_poly(subset(aes_245, 4), q = quant) |> st_area() |> sum() |> units::set_units("km^2"),
                                     aes_poly(subset(aes_245, 5), q = quant) |> st_area() |> sum() |> units::set_units("km^2"),
                                     aes_poly(subset(aes_245, 6), q = quant) |> st_area() |> sum() |> units::set_units("km^2"),
                                     aes_poly(subset(aes_245, 7), q = quant) |> st_area() |> sum() |> units::set_units("km^2"),
                                     aes_poly(subset(aes_245, 8), q = quant) |> st_area() |> sum() |> units::set_units("km^2")),
                          SSP585 = c(aes_poly(subset(aes_585, 1), q = quant) |> st_area() |> sum() |> units::set_units("km^2"),
                                     aes_poly(subset(aes_585, 2), q = quant) |> st_area() |> sum() |> units::set_units("km^2"),
                                     aes_poly(subset(aes_585, 3), q = quant) |> st_area() |> sum() |> units::set_units("km^2"),
                                     aes_poly(subset(aes_585, 4), q = quant) |> st_area() |> sum() |> units::set_units("km^2"),
                                     aes_poly(subset(aes_585, 5), q = quant) |> st_area() |> sum() |> units::set_units("km^2"),
                                     aes_poly(subset(aes_585, 6), q = quant) |> st_area() |> sum() |> units::set_units("km^2"),
                                     aes_poly(subset(aes_585, 7), q = quant) |> st_area() |> sum() |> units::set_units("km^2"),
                                     aes_poly(subset(aes_585, 8), q = quant) |> st_area() |> sum() |> units::set_units("km^2"))
                          )

future_aes_size <- future_aes_size |> mutate(change_245 = ((SSP245 - current_aes_size) / current_aes_size) * 100,
                                             change_585 = ((SSP585 - current_aes_size) / current_aes_size) * 100)
future_aes_size |> summarise(mean(change_245), sd(change_245),
                             mean(change_585), sd(change_585))

# calculate absolute change in aes habitat suitability
aes_245_change <- (aes_245 - stack_now_max)
aes_585_change <- (aes_585 - stack_now_max)
names(aes_245_change) <- names(aes_245)
names(aes_585_change) <- names(aes_585)

# calculate guild-level change....
aes_245_fish  <- stack(subset(stack_245_1, 1), subset(stack_245_2, 1), subset(stack_245_3, 1), subset(stack_245_4, 1),
                       subset(stack_245_5, 1), subset(stack_245_6, 1), subset(stack_245_7, 1), subset(stack_245_8, 1))
aes_245_krill <- stack(subset(stack_245_1, 2), subset(stack_245_2, 2), subset(stack_245_3, 2), subset(stack_245_4, 2),
                       subset(stack_245_5, 2), subset(stack_245_6, 2), subset(stack_245_7, 2), subset(stack_245_8, 2))
aes_245_squid <- stack(subset(stack_245_1, 3), subset(stack_245_2, 3), subset(stack_245_3, 3), subset(stack_245_4, 3),
                       subset(stack_245_5, 3), subset(stack_245_6, 3), subset(stack_245_7, 3), subset(stack_245_8, 3))
aes_585_fish  <- stack(subset(stack_585_1, 1), subset(stack_585_2, 1), subset(stack_585_3, 1), subset(stack_585_4, 1),
                       subset(stack_585_5, 1), subset(stack_585_6, 1), subset(stack_585_7, 1), subset(stack_585_8, 1))
aes_585_krill <- stack(subset(stack_585_1, 2), subset(stack_585_2, 2), subset(stack_585_3, 2), subset(stack_585_4, 2),
                       subset(stack_585_5, 2), subset(stack_585_6, 2), subset(stack_585_7, 2), subset(stack_585_8, 2))
aes_585_squid <- stack(subset(stack_585_1, 3), subset(stack_585_2, 3), subset(stack_585_3, 3), subset(stack_585_4, 3),
                       subset(stack_585_5, 3), subset(stack_585_6, 3), subset(stack_585_7, 3), subset(stack_585_8, 3))
names(aes_245_fish) <- names(aes_245)
names(aes_245_krill) <- names(aes_245)
names(aes_245_squid) <- names(aes_245)
names(aes_585_fish) <- names(aes_245)
names(aes_585_krill) <- names(aes_245)
names(aes_585_squid) <- names(aes_245)

# guild-level absolute change
aes_245_fish_change <- aes_245_fish - subset(stack_now_guilds, 1)
aes_585_fish_change <- aes_585_fish - subset(stack_now_guilds, 1)

aes_245_krill_change <- aes_245_krill - subset(stack_now_guilds, 2)
aes_585_krill_change <- aes_585_krill - subset(stack_now_guilds, 2)

aes_245_squid_change <- aes_245_squid - subset(stack_now_guilds, 3)
aes_585_squid_change <- aes_585_squid - subset(stack_now_guilds, 3)

future_aes_guild <- stack(mean(aes_245_fish_change),
                          mean(aes_245_krill_change),
                          mean(aes_245_squid_change),
                          mean(aes_585_fish_change),
                          mean(aes_585_krill_change),
                          mean(aes_585_squid_change))
names(future_aes_guild) <- c("Fish 245", "Krill 245", "Squid 245", "Fish 585", "Krill 585", "Squid 585")
future_aes_guild_df <- future_aes_guild |> rasterToPoints() |>  as_tibble() |> pivot_longer(3:8, names_to = "group", values_to = "preds")
future_aes_guild_df <- future_aes_guild_df |> separate_wider_delim(group, ".", names = c("guild", "scenario"))

future_aes_guild_df <- future_aes_guild_df |> 
  mutate(scenario = case_when(scenario == "245" ~ "SSP2-4.5",
                              scenario == "585" ~ "SSP5-8.5"))

########################
### Output the plots ###
########################

# some plots
p1 <- ggplot() + 
  ggtitle("SSP2-4.5") +
  theme_void(base_size = 7,
             base_family = "Helvetica Neue") +
  geom_raster(aes(x = x, y = y, fill = layer), data = mean(aes_245_change) |> rasterToPoints()) +
  geom_sf(aes(), data = world_shp, colour = "grey60", fill = "grey60") +
  geom_sf(aes(), data = polar_buffer$mask, fill = "white", color = "grey40", size = 0.5 / .pt) +
  geom_sf(aes(), data = contemporary_aes, fill = NA, colour = "black") +
  coord_sf(xlim = c(-6400000, 6400000), ylim = c(-6400000, 6400000), expand = FALSE, crs = crs_polar, ndiscr = 1000) +
  scale_fill_discrete_gradient("Absolute Change in Habitat Importance",
                               colours =  RColorBrewer::brewer.pal(10, "RdYlBu") [c(3:10)],
                               bins = 10,
                               limits = c(-.4, .6),
                               breaks = seq(-.4, .6, by = .2),
                               labels = seq(-.4, .6, by = .2),
                               oob = scales::squish,
                               guide = guide_colourbar(
                                 nbin = 500,
                                 display = "raster",
                                 frame.colour = "grey40",
                                 ticks.colour = "grey40",
                                 frame.linewidth = .1,
                                 barwidth = .5,
                                 barheight = 10,
                                 direction = "vertical",
                                 title.position = "right",
                                 title.theme = element_text(hjust = 0.5, size = 7, angle = 90))) +
  xlab(NULL) + ylab(NULL) +
  theme(plot.title = element_text(hjust = 0.5, size = 7, family = "Helvetica"),
        legend.text = element_text(size = 7, family = "Helvetica"),
        legend.title = element_text(size = 7, family = "Helvetica"),
        panel.background = element_rect(fill = "black"))


p2 <- ggplot() + 
  ggtitle("SSP5-8.5") +
  theme_void(base_size = 7,
             base_family = "Helvetica") +
  geom_raster(aes(x = x, y = y, fill = layer), data = mean(aes_585_change) |> rasterToPoints()) +
  geom_sf(aes(), data = world_shp, colour = "grey60", fill = "grey60") +
  geom_sf(aes(), data = polar_buffer$mask, fill = "white", color = "grey40", size = 0.5 / .pt) +
  geom_sf(aes(), data = contemporary_aes, fill = NA, colour = "black") +
  coord_sf(xlim = c(-6400000, 6400000), ylim = c(-6400000, 6400000), expand = FALSE, crs = crs_polar, ndiscr = 1000) +
  scale_fill_discrete_gradient("Absolute Change in Habitat Importance",
                               colours =  RColorBrewer::brewer.pal(10, "RdYlBu") [c(3:10)],
                               bins = 10,
                               limits = c(-.4, .6),
                               breaks = seq(-.4, .6, by = .2),
                               labels = seq(-.4, .6, by = .2),
                               oob = scales::squish,
                               guide = guide_colourbar(
                                 nbin = 500,
                                 display = "raster",
                                 frame.colour = "grey40",
                                 ticks.colour = "grey40",
                                 frame.linewidth = .1,
                                 barwidth = .5,
                                 barheight = 10,
                                 direction = "vertical",
                                 title.position = "right",
                                 title.theme = element_text(hjust = 0.5, size = 7, angle = 90))) +
  xlab(NULL) + ylab(NULL) +
  theme(plot.title = element_text(hjust = 0.5, size = 7, family = "Helvetica"),
        legend.text = element_text(size = 7, family = "Helvetica"),
        legend.title = element_text(size = 7, family = "Helvetica"),
        panel.background = element_rect(fill = "black"))

p3 <- ggplot() + 
  theme_void(base_size = 7,
             base_family = "Helvetica") +
  geom_raster(aes(x = x, y = y, fill = preds), data = future_aes_guild_df) +
  facet_grid(guild ~ scenario) +
  geom_sf(aes(), data = world_shp, colour = "grey60", fill = "grey60") +
  geom_sf(aes(), data = polar_buffer$mask, fill = "white", color = "grey40", size = 0.5 / .pt) +
#  geom_sf(aes(), data = contemporary_aes, fill = NA, colour = "black") +
  coord_sf(xlim = c(-6400000, 6400000), ylim = c(-6400000, 6400000), expand = FALSE, crs = crs_polar, ndiscr = 1000) +
  scale_fill_discrete_gradient("Absolute Change in Habitat Suitability",
                               colours =  RColorBrewer::brewer.pal(10, "RdYlBu") [c(3:10)],
                               bins = 10,
                               limits = c(-.4, .6),
                               breaks = seq(-.4, .6, by = .2),
                               labels = seq(-.4, .6, by = .2),
                               oob = scales::squish,
                               guide = guide_colourbar(
                                 nbin = 500,
                                 display = "raster",
                                 frame.colour = "grey40",
                                 ticks.colour = "grey40",
                                 frame.linewidth = .1,
                                 barwidth = 10,
                                 barheight = .5,
                                 direction = "horizontal",
                                 title.position = "top",
                                 title.theme = element_text(hjust = 0.5, size = 7))) +
  xlab(NULL) + ylab(NULL) +
  theme(plot.title = element_text(hjust = 0.5, size = 7, family = "Helvetica"),
        legend.text = element_text(size = 7, family = "Helvetica"),
        legend.title = element_text(size = 7, family = "Helvetica"),
        panel.background = element_rect(fill = "black"),
        legend.position = "bottom")

# horizontal panel 3
p3 <- ggplot() + 
  theme_void(base_size = 7,
             base_family = "Helvetica") +
  geom_raster(aes(x = x, y = y, fill = preds), data = future_aes_guild_df) +
  facet_grid(scenario ~ guild, switch = "y") +
  geom_sf(aes(), data = world_shp, colour = "grey60", fill = "grey60") +
  geom_sf(aes(), data = polar_buffer$mask, fill = "white", color = "grey40", size = 0.5 / .pt) +
  coord_sf(xlim = c(-6400000, 6400000), ylim = c(-6400000, 6400000), expand = FALSE, crs = crs_polar, ndiscr = 1000) +
  scale_fill_discrete_gradient("Absolute Change in Habitat Suitability",
                               colours =  RColorBrewer::brewer.pal(10, "RdYlBu") [c(3:10)],
                               bins = 10,
                               limits = c(-.4, .6),
                               breaks = seq(-.4, .6, by = .2),
                               labels = seq(-.4, .6, by = .2),
                               oob = scales::squish,
                               guide = guide_colourbar(
                                 nbin = 500,
                                 display = "raster",
                                 frame.colour = "grey40",
                                 ticks.colour = "grey40",
                                 frame.linewidth = .1,
                                 barwidth = .5,
                                 barheight = 10,
                                 direction = "vertical",
                                 title.position = "right",
                                 title.theme = element_text(hjust = 0.5, size = 7, angle = 90)))  +
  xlab(NULL) + ylab(NULL) +
  theme(plot.title = element_text(hjust = 0.5, size = 7, family = "Helvetica"),
        legend.text = element_text(size = 7, family = "Helvetica"),
        legend.title = element_text(size = 7, family = "Helvetica"),
        strip.text.x = element_text(size = 7, family = "Helvetica"),
        strip.text.y.left = element_text(size = 7, family = "Helvetica", angle = 90),
        panel.background = element_rect(fill = "black"),
        legend.position = "right")

#####################################
### Range and suitability changes ###
#####################################

# function to calculate raster thresholds from ROC curves
raster_threshold <- function(input_raster_path, sp){
  
  # load raster
  input_raster <- readRDS(input_raster_path)
  
  # load original data
  dat <- readRDS("data/presence_absence_data_10k_with_covariates_2024-06-17.rds")
  dat <- dat |> filter(species == unique(species)[sp]) # focal species
  
  # extract probabilities from raster and append to data
  est_prob <- raster::extract(input_raster, cbind(dat$x, dat$y))
  
  roc_obj <- pROC::roc(dat$PresAbs, est_prob)   # calculate ROC curve
  roc_threshold <- pROC::coords(roc_obj, x = "best")$threshold  # pull out threshold based on Youden
  
  return(roc_threshold)
}

# apply function to each species
thresholds <- map2(files_now, 1:28, raster_threshold)
thresholds <- thresholds |> unlist()

# create binary raster stack based on species specific thresholding
x <- stack_now > thresholds

# calculate the distance to the south pole for each raster cell
distance_now <- x |> 
  rasterToPoints() |> 
  as_tibble() |> 
  pivot_longer(3:30, names_to = "species", values_to = "preds") |>
  filter(preds > 0) |>
  mutate(distance = sqrt(x^2 + y^2)) |>
  group_by(species) |>
  summarise(distance = mean(distance)/1000)

# replica of the futures stacks
stack_245_threshold <- stack_245
stack_585_threshold <- stack_585

# threshold each species and GCM combo
for(i in 1:length(thresholds)){
  stack_245_threshold[[i]] <- stack_245[[i]] > thresholds[i]
  stack_585_threshold[[i]] <- stack_585[[i]] > thresholds[i]
}

# convert stack of rasters and GCMS to a tibble
range_shift_245 <- lapply(stack_245_threshold, function(x){ x |> rasterToPoints() |> as_tibble() |> 
    pivot_longer(3:10, names_to = "model", values_to = "preds") |>
    filter(preds > 0) |>
    group_by(model)
}
)

range_shift_585 <- lapply(stack_585_threshold, function(x){ x |> rasterToPoints() |> as_tibble() |> 
    pivot_longer(3:10, names_to = "model", values_to = "preds") |>
    filter(preds > 0) |>
    group_by(model)
}
)

range_shift_245 <- range_shift_245 |> 
  tibble() |>
  mutate(species = names(x)) |>
  unnest(cols = c(range_shift_245))

range_shift_585 <- range_shift_585 |> 
  tibble() |>
  mutate(species = names(x)) |>
  unnest(cols = c(range_shift_585))

distance_now <- distance_now |> 
  mutate(guild = case_when(species %in% c("species01", "species02", "species03", "species04",
                                          "species05", "species06", "species07", "species08",
                                          "species09", "species10", "species11", "species12",
                                          "species13") ~ "Fish",
                           species == "species14" ~ "Krill",
                           species %in% c("species15", "species16", "species17", "species18",
                                          "species19", "species20", "species21", "species22",
                                          "species23", "species24", "species25", "species26", 
                                          "species27", "species28") ~ "Squid"))

range_shift_245 <- range_shift_245 |> 
  mutate(guild = case_when(species %in% c("species01", "species02", "species03", "species04",
                                          "species05", "species06", "species07", "species08",
                                          "species09", "species10", "species11", "species12",
                                          "species13") ~ "Fish",
                           species == "species14" ~ "Krill",
                           species %in% c("species15", "species16", "species17", "species18",
                                          "species19", "species20", "species21", "species22",
                                          "species23", "species24", "species25", "species26", 
                                          "species27", "species28") ~ "Squid"))

range_shift_585 <- range_shift_585 |> 
  mutate(guild = case_when(species %in% c("species01", "species02", "species03", "species04",
                                          "species05", "species06", "species07", "species08",
                                          "species09", "species10", "species11", "species12",
                                          "species13") ~ "Fish",
                           species == "species14" ~ "Krill",
                           species %in% c("species15", "species16", "species17", "species18",
                                          "species19", "species20", "species21", "species22",
                                          "species23", "species24", "species25", "species26", 
                                          "species27", "species28") ~ "Squid"))

species_names <- c("B. antarcticus", "N. coatsorum", "P. antarctica", "E. antarctica", "E. carlsbergi", "G. bolini", "G. braueri",
                   "G. fraseri", "G. nicholsi", "G. opisthopterus", "K. anderssoni", "P. bolini", "P. tenisoni", "E. superba",
                   "A. antarcticus", "B. abyssicola", "G. glacialis", "G. antarcticus", "H. atlantica", "H. eltaninae", "K. longimana",
                   "M. hyadesi", "M. hamiltoni", "M. ingens", "M. robsoni", "P. glacialis", "S. circumantarctica", "T. filippovae")     

distance_now <- distance_now |> mutate(species = factor(species))
range_shift_245 <- range_shift_245 |> mutate(species = factor(species))
range_shift_585 <- range_shift_585 |> mutate(species = factor(species))

levels(distance_now$species) <- species_names
levels(range_shift_245$species) <- species_names
levels(range_shift_585$species) <- species_names

range_shift_245 <- range_shift_245 |> mutate(scenario = "SSP2-4.5")
range_shift_585 <- range_shift_585 |> mutate(scenario = "SSP5-8.5")

range_shift_245 <- range_shift_245 |> 
  mutate(distance = sqrt(x^2 + y^2)) |>
  group_by(species, model) |>
  summarise(distance = mean(distance)/1000)

range_shift_585 <- range_shift_585 |> 
  mutate(distance = sqrt(x^2 + y^2)) |>
  group_by(species, model) |>
  summarise(distance = mean(distance)/1000)

# combine future estimates of mean distance to present day
# difference between two is equivalent to a range shift
range_shift_245 <- range_shift_245 |> left_join(distance_now, by = "species")
range_shift_585 <- range_shift_585 |> left_join(distance_now, by = "species")
range_shift_245 <- range_shift_245 |> mutate(shift = distance.x - distance.y)
range_shift_585 <- range_shift_585 |> mutate(shift = distance.x - distance.y)

range_shift_245 <- range_shift_245 |> mutate(scenario = "SSP2-4.5")
range_shift_585 <- range_shift_585 |> mutate(scenario = "SSP5-8.5")
range_shift <- rbind(range_shift_245, range_shift_585)

range_shift_2 <- range_shift |> group_by(guild, model, scenario) |> mutate(shift = mean(shift))
range_shift_2 <- range_shift_2 |> dplyr:::select(guild, model, scenario, shift) |> unique()
range_shift_2 <- range_shift_2 |> ungroup()
range_shift_2 <- range_shift_2 |> arrange(guild, model, scenario)

p4 <- ggplot() + 
  theme_bw(base_size = 7,
           base_family = "Helvetica") +
  geom_boxplot(aes(x = guild, y = shift, colour = scenario, fill = scenario), 
               width = .2, outliers = F, position = position_dodge(-.7),
               data = range_shift_2) +
  scale_fill_manual(NULL, labels = c("SSP2-4.5", "SSP5-8.5"), values = wesanderson::wes_palette("AsteroidCity1")[2:3]) +
  scale_colour_manual(NULL, labels = c("SSP2-4.5", "SSP5-8.5"), values = wesanderson::wes_palette("AsteroidCity1")[2:3]) +
  geom_point(aes(x = guild, y = shift, group = scenario), 
             position = position_dodge(-.7), size = .5,
             data = range_shift_2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = "Mean range shift (km)",
                     limits = c(-550, 25), expand = c(0, 0)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7, family = "Helvetica"),
        legend.key.size = unit(.4, 'cm'),
        axis.title = element_text(size = 7, family = "Helvetica"),
        axis.text = element_text(size = 7, family = "Helvetica"))


##################################
### change in suitable habitat ###
##################################

# count number of thresholded cells
habitat_245 <- lapply(stack_245_threshold, function(x){ x |> rasterToPoints() |> as_tibble() |> 
    pivot_longer(3:10, names_to = "model", values_to = "preds") |>
    filter(preds > 0) |>
    group_by(model) |>
    tally()
}
)

habitat_585 <- lapply(stack_585_threshold, function(x){ x |> rasterToPoints() |> as_tibble() |> 
    pivot_longer(3:10, names_to = "model", values_to = "preds") |>
    filter(preds > 0) |>
    group_by(model) |>
    tally()
}
)

habitat_now <- x |> 
  rasterToPoints() |> 
  as_tibble() |> 
  pivot_longer(3:30, names_to = "species", values_to = "preds") |>
  group_by(species) |>
  filter(preds > 0) |>
  tally()

habitat_245 <- habitat_245 |> 
  tibble() |>
  mutate(species = names(x),
         size = habitat_now$n) |>
  unnest(cols = c(habitat_245))

habitat_585 <- habitat_585 |> 
  tibble() |>
  mutate(species = names(x),
         size = habitat_now$n) |>
  unnest(cols = c(habitat_585))

habitat_245 <- habitat_245 |> 
  mutate(guild = case_when(species %in% c("species01", "species02", "species03", "species04",
                                          "species05", "species06", "species07", "species08",
                                          "species09", "species10", "species11", "species12",
                                          "species13") ~ "Fish",
                           species == "species14" ~ "Krill",
                           species %in% c("species15", "species16", "species17", "species18",
                                          "species19", "species20", "species21", "species22",
                                          "species23", "species24", "species25", "species26", 
                                          "species27", "species28") ~ "Squid"))

habitat_585 <- habitat_585 |> 
  mutate(guild = case_when(species %in% c("species01", "species02", "species03", "species04",
                                          "species05", "species06", "species07", "species08",
                                          "species09", "species10", "species11", "species12",
                                          "species13") ~ "Fish",
                           species == "species14" ~ "Krill",
                           species %in% c("species15", "species16", "species17", "species18",
                                          "species19", "species20", "species21", "species22",
                                          "species23", "species24", "species25", "species26", 
                                          "species27", "species28") ~ "Squid"))

species_names <- c("B. antarcticus", "N. coatsorum", "P. antarctica", "E. antarctica", "E. carlsbergi", "G. bolini", "G. braueri",
                   "G. fraseri", "G. nicholsi", "G. opisthopterus", "K. anderssoni", "P. bolini", "P. tenisoni", "E. superba",
                   "A. antarcticus", "B. abyssicola", "G. glacialis", "G. antarcticus", "H. atlantica", "H. eltaninae", "K. longimana",
                   "M. hyadesi", "M. hamiltoni", "M. ingens", "M. robsoni", "P. glacialis", "S. circumantarctica", "T. filippovae")     

habitat_245 <- habitat_245 |> mutate(species = factor(species))
habitat_585 <- habitat_585 |> mutate(species = factor(species))
levels(habitat_245$species) <- species_names
levels(habitat_585$species) <- species_names
habitat_245 <- habitat_245 |> mutate(scenario = "SSP2-4.5")
habitat_585 <- habitat_585 |> mutate(scenario = "SSP5-8.5")
future_habitat <- rbind(habitat_245, habitat_585)

plot_df <- future_habitat |> group_by(guild, model, scenario) |> summarise(change = mean((n-size)/size))

p5 <- ggplot() + 
  theme_bw(base_size = 7,
           base_family = "Helvetica") +
  geom_boxplot(aes(x = guild, y = change*100, colour = scenario, fill = scenario), 
               width = .2, outliers = F, position = position_dodge(-.7),
               data = plot_df) +
  scale_fill_manual(NULL, values = wesanderson::wes_palette("AsteroidCity1")[2:3]) +
  scale_colour_manual(NULL, values = wesanderson::wes_palette("AsteroidCity1")[2:3]) +
  geom_point(aes(x = guild, y = change*100, group = scenario), 
             position = position_dodge(-.7), size = .5,
             data = plot_df) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = "Change in Area of \nSuitable Habitat (%)",
                     limits = c(-50, 50), expand = c(0, 0)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7, family = "Helvetica"),
        legend.key.size = unit(.4, 'cm'),
        axis.title = element_text(size = 7, family = "Helvetica"),
        axis.text = element_text(size = 7, family = "Helvetica"))

p <- (p1 + p2 + plot_layout(guides = "collect")) / (p3 + (p4 / p5 + plot_layout(guides = "collect") & theme(legend.position = "bottom")) + plot_layout(widths = c(3, 2))) +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 9,
                                                                  family = "Helvetica Neue",
                                                                  face = "bold",
                                                                  colour = "black"))

ggsave(filename = "Figure 4 absolute change 2025-02-26.jpeg",
       plot = p,
       path = "figures/",
       width = 140,
       height = 180,
       units = "mm",
       dpi = 500)



foo <- (((p1 & theme(legend.position = "none")) / (p2 & theme(legend.position = "none"))) | (p3 / (p4 + p5 + plot_layout(guides = "collect")) + plot_layout(heights = c(2, 1)))) + 
  plot_layout(widths = c(2, 3)) +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 9,
                                                                    family = "Helvetica Neue",
                                                                    face = "bold",
                                                                    colour = "black"))
ggsave(filename = "Figure 4 landscape 2025-02-26.jpeg",
       plot = foo,
       path = "figures/",
       width = 180,
       height = 120,
       units = "mm",
       dpi = 500)


(p1 / p2 + plot_layout(guides = "collect")) + (p3 / (p4 + p5 + plot_layout(guides = "collect"))) & theme(legend.position = "bottom"))) 
                                               
                                               + plot_layout(widths = c(3, 2))) +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 9,
                                                                    family = "Helvetica Neue",
                                                                    face = "bold",
                                                                    colour = "black"))




#################################################################
### Output species-specific plots for supplementary materials ###
#################################################################

future_habitat |>
  mutate(change = ((n-size)/size)*100) |>
  group_by(species, scenario) |>
  summarise(median(change)) |>
  print(n = Inf)

p1 <- ggplot() + 
  ggtitle("Fish") +
  theme_bw(base_size = 8,
           base_family = "Helvetica Neue") +
  geom_boxplot(aes(x = ((n-size)/size)*100, y = species, colour = scenario, fill = scenario), 
               width = .5, outliers = F, position = position_dodge(-.7),
               data = future_habitat |> filter(guild == "Fish")) +
  scale_fill_manual(values = wesanderson::wes_palette("AsteroidCity1")[2:3]) +
  scale_colour_manual(values = wesanderson::wes_palette("AsteroidCity1")[2:3]) +
  geom_point(aes(x = ((n-size)/size)*100, y = species, group = scenario), 
             position = position_dodge(-.7), size = .5,
             data = future_habitat |> filter(guild == "Fish")) +
  scale_y_discrete(name = NULL, limits = rev) +
  scale_x_continuous(name = NULL,
                     limits = c(-100, 40), expand = c(0, 0)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"),
        plot.title = element_text(size = 8,
                                  family = "Helvetica Neue"))

p2 <- ggplot() + 
  ggtitle("Krill") +
  theme_bw(base_size = 8,
           base_family = "Helvetica Neue") +
  geom_boxplot(aes(x = ((n-size)/size)*100, y = species, colour = scenario, fill = scenario), 
               width = .5, outliers = F, position = position_dodge(-.7),
               data = future_habitat |> filter(guild == "Krill")) +
  scale_fill_manual(values = wesanderson::wes_palette("AsteroidCity1")[2:3]) +
  scale_colour_manual(values = wesanderson::wes_palette("AsteroidCity1")[2:3]) +
  geom_point(aes(x = ((n-size)/size)*100, y = species, group = scenario), 
             position = position_dodge(-.7), size = .5,
             data = future_habitat |> filter(guild == "Krill")) +
  scale_y_discrete(name = NULL, limits = rev) +
  scale_x_continuous(name = "Change in Area of Suitable Habitat (%)",
                     limits = c(-100, 40), expand = c(0, 0)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"),
        plot.title = element_text(size = 8,
                                  family = "Helvetica Neue"))

#my_order <- trial |> 
#  filter(guild == "squid") |> 
#  filter(!species %in% c("T. filippovae", "S. circumantarctica")) |> 
#  mutate(change = ((n-size)/size)*100) |> 
#  group_by(species) |> 
#  summarise(m_c = mean(change)) |> 
#  arrange(m_c) |> 
#  pull(species)

p3 <- ggplot() + 
  ggtitle("Squid") +
  theme_bw(base_size = 8,
           base_family = "Helvetica Neue") +
  geom_boxplot(aes(x = ((n-size)/size)*100, y = species, colour = scenario, fill = scenario), 
               width = .5, outliers = F, position = position_dodge(-.7),
               data = future_habitat |> filter(guild == "Squid") |> filter(!species %in% c("T. filippovae", "S. circumantarctica"))) +
  scale_fill_manual(values = wesanderson::wes_palette("AsteroidCity1")[2:3]) +
  scale_colour_manual(values = wesanderson::wes_palette("AsteroidCity1")[2:3]) +
  geom_point(aes(x = ((n-size)/size)*100, y = species, group = scenario), 
             position = position_dodge(-.7), size = .5,
             data = future_habitat |> filter(guild == "Squid") |> filter(!species %in% c("T. filippovae", "S. circumantarctica"))) +
  scale_y_discrete(name = NULL, limits = rev) +
  scale_x_continuous(name = NULL, limits = c(-60, 85), expand = c(0, 0)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"),
        plot.title = element_text(size = 8,
                                  family = "Helvetica Neue"))

p4 <- ggplot() + 
  ggtitle("Squid cont.") +
  theme_bw(base_size = 8,
           base_family = "Helvetica Neue") +
  geom_boxplot(aes(x = ((n-size)/size)*100, y = species, colour = scenario, fill = scenario), 
               width = .5, outliers = F, position = position_dodge(-.7),
               data = future_habitat |> filter(species %in% c("T. filippovae", "S. circumantarctica"))) +
  scale_fill_manual(values = wesanderson::wes_palette("AsteroidCity1")[2:3]) +
  scale_colour_manual(values = wesanderson::wes_palette("AsteroidCity1")[2:3]) +
  geom_point(aes(x = ((n-size)/size)*100, y = species, group = scenario), 
             position = position_dodge(-.7), size = .5,
             data = future_habitat |> filter(species %in% c("T. filippovae", "S. circumantarctica"))) +
  scale_y_discrete(name = NULL, limits = rev) +
  scale_x_continuous(name = "Change in Area of Suitable Habitat (%)",
                     limits = c(0, 360), expand = c(0, 0)) +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"),
        plot.title = element_text(size = 8,
                                  family = "Helvetica Neue"))

p5 <- p1 + p2 + plot_layout(heights = c(14, 1))
p6 <- p3 + p4 + plot_layout(heights = c(13, 2))

p <- p5 | p6

ggsave(filename = "all in boxplot.jpeg",
       plot = p,
       path = "figures/",
       width = 180,
       height = 140,
       units = "mm",
       dpi = 500)


# end






# calculate percentage change in aes habitat suitability
#aes_245_percent_change <- ((aes_245 - stack_now_max) / stack_now_max) * 100
#aes_585_percent_change <- ((aes_585 - stack_now_max) / stack_now_max) * 100
#names(aes_245_percent_change) <- names(aes_245)
#names(aes_585_percent_change) <- names(aes_585)
#plot(mean(aes_245_percent_change))
#plot(mean(aes_585_percent_change))
#aes_245_fish_percent_change <- ((aes_245_fish - subset(stack_now_guilds, 1)) / subset(stack_now_guilds, 1)) * 100
#aes_585_fish_percent_change <- ((aes_585_fish - subset(stack_now_guilds, 1)) / subset(stack_now_guilds, 1)) * 100
#aes_245_krill_percent_change <- ((aes_245_krill - subset(stack_now_guilds, 2)) / subset(stack_now_guilds, 2)) * 100
#aes_585_krill_percent_change <- ((aes_585_krill - subset(stack_now_guilds, 2)) / subset(stack_now_guilds, 2)) * 100
#aes_245_squid_percent_change <- ((aes_245_squid - subset(stack_now_guilds, 3)) / subset(stack_now_guilds, 3)) * 100
#aes_585_squid_percent_change <- ((aes_585_squid - subset(stack_now_guilds, 3)) / subset(stack_now_guilds, 3)) * 100
# krill percentage change is a bit wild...
#foo <- mean(aes_245_krill_percent_change)
#y <- focal(foo, w=matrix(1, 3, 3), mean, na.rm = T)
#plot(y)


