################################
### Identify climate refugia ###
################################

# 2025-03-02

# Are there regions that remain relatively buffered from the effects of climate change?
# These could act as refugia for micronekton

# Based on exisiting analysis, examine overlap between present day and future AESs
# If overlap is high, then AESs represent micronekton habitat that is stable in the future

# Alternative is to calculate "in situ" or "ex situ" climate refugia
# What regions that are currently good remain good?
# What regions that are currently good become unsuitable?
# What regions that are currently unsuitable become suitable?
# How much of this area is within micronekton AES?

# libraries
require(raster)
require(tidyverse)
require(sf)
sf::sf_use_s2(FALSE)
require(patchwork)
source("R/aes_poly.R")

# define projection
prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" 

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

###
### Approach 1 - Calculate overlap between present and future AES
###

# given the contemporary raster and stack of future GCM inputs
# output the overlap
future_overlap <- function(contemporary_raster, future_raster){
  
  # calculate contemporary AES
  contemporary_aes <- aes_poly(contemporary_raster, prob = .9)
  
  # extract top decile value used for contemporary AES
  # Use that for the future so if habitat gets worse AES will reduce
  quant <- quantile(contemporary_raster, p = .9)
  
  # guild specific AES in future based on decile value now
  future_aes <- aes_poly(future_raster, q = quant)
  
  # how much of the contemporary AES is within future AESs?
  overlap <- st_area(st_intersection(contemporary_aes, future_aes)) / st_area(contemporary_aes)
  
  return(overlap)
}

# create output table containing the estimated overlap between AES area at end of century and present day
aes_refugia <- tibble(scenario = c(rep("SSP245", 8), rep("SSP585", 8)),
                      model =  rep(c("ACCESS.CM2", "BCC.CSM2.MR", "CESM2.WACCM", "CMCC.CM2.SR5",
                                     "FGOALS.g3", "IPSL.CM6A.LR", "MRI.ESM2.0", "NESM3"), 2),
                      overlap = c(future_overlap(stack_now_max, subset(aes_245, 1)),
                                  future_overlap(stack_now_max, subset(aes_245, 2)),
                                  future_overlap(stack_now_max, subset(aes_245, 3)),
                                  future_overlap(stack_now_max, subset(aes_245, 4)),
                                  future_overlap(stack_now_max, subset(aes_245, 5)),
                                  future_overlap(stack_now_max, subset(aes_245, 6)),
                                  future_overlap(stack_now_max, subset(aes_245, 7)),
                                  future_overlap(stack_now_max, subset(aes_245, 8)),
                                  future_overlap(stack_now_max, subset(aes_585, 1)),
                                  future_overlap(stack_now_max, subset(aes_585, 2)),
                                  future_overlap(stack_now_max, subset(aes_585, 3)),
                                  future_overlap(stack_now_max, subset(aes_585, 4)),
                                  future_overlap(stack_now_max, subset(aes_585, 5)),
                                  future_overlap(stack_now_max, subset(aes_585, 6)),
                                  future_overlap(stack_now_max, subset(aes_585, 7)),
                                  future_overlap(stack_now_max, subset(aes_585, 8))))

aes_refugia |> group_by(scenario) |> summarise(mean(overlap), sd(overlap))

# overlap between future AES and contemporary AES is high

###
### Approach 2 - in situ climate refugia
###

# identify which cells estimated to be presence locations now
# are still presence locations at the end of the century

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
x <- stack(x)

# replica of the futures stacks
stack_245_threshold <- stack_245
stack_585_threshold <- stack_585

# threshold each species and GCM combo
for(i in 1:length(thresholds)){
  stack_245_threshold[[i]] <- stack_245[[i]] > thresholds[i]
  stack_585_threshold[[i]] <- stack_585[[i]] > thresholds[i]
}

# function
# identify cells that are above threshold now and at end of century
# tally agreement between GCMs
# return taster layer of 0s and 1s showing presence cells
refugia_tally <- function(dist_now, dist_future, agreement){
  dist_future <- sum(dist_future)
  dist_future <- dist_future >= agreement
  dist_future[dist_now == 0] <- 0
  return(dist_future)
}

# use help function and loop through all 28 species
# could probably vectorise this
# output is a stack representing model agreement for each species
refugia_245 <- stack()
for(i in 1:28){
  stack_loop <- refugia_tally(dist_now = subset(x, i), dist_future = stack_245_threshold[[i]], agreement = 8)
  refugia_245 <- stack(refugia_245, stack_loop)
}
names(refugia_245) <- names(x)

refugia_585 <- stack()
for(i in 1:28){
  stack_loop <- refugia_tally(dist_now = subset(x, i), dist_future = stack_585_threshold[[i]], agreement = 8)
  refugia_585 <- stack(refugia_585, stack_loop)
}
names(refugia_585) <- names(x)

# now have a raster layer containing cells that are good now and at the end of the century
# Are there any cells that contain all three guilds now and at the end of the century?

guilds <- c(rep(1, times = 13), 2, rep(3, times = 14))

foo <- stackApply(x, guilds, sum)
foo[foo > 0] <- 1
foo <- sum(foo)
foo[foo == 1] <- 0
foo[foo == 2] <- 0
foo[foo == 3] <- 1
plot(foo)

# what percentage of the AES is also currently candidate refugia?
# so what percentage of AES is covered by green?
foo_mask <- mask(foo, as_Spatial(contemporary_aes))
foo_mask_values <- getValues(foo_mask)
foo_mask_values <- foo_mask_values[!is.na(foo_mask_values)] # drop na
sum(foo_mask_values)/length(foo_mask_values) #74.2% of AES is refugia now

foo_values <- getValues(foo)
foo_values <- foo_values[!is.na(foo_values)]

sum(foo_mask_values)/sum(foo_values)

refugia_245_guilds_mask_values <- getValues(refugia_245_guilds_mask)
refugia_585_guilds_mask_values <- getValues(refugia_585_guilds_mask)
refugia_245_guilds_mask_values <- refugia_245_guilds_mask_values[!is.na(refugia_245_guilds_mask_values)] # drop na
refugia_585_guilds_mask_values <- refugia_585_guilds_mask_values[!is.na(refugia_585_guilds_mask_values)] # drop na
sum(refugia_245_guilds_mask_values)/length(refugia_245_guilds_mask_values) #74.2% of AES is refugia now
sum(refugia_585_guilds_mask_values)/length(refugia_585_guilds_mask_values) #74.2% of AES is refugia now

plot(refugia_585_guilds)
plot(refugia_245_guilds)

refugia_245_guilds <- stackApply(refugia_245, guilds, sum)
refugia_585_guilds <- stackApply(refugia_585, guilds, sum)

refugia_245_guilds[refugia_245_guilds > 0] <- 1
refugia_245_guilds <- sum(refugia_245_guilds)
refugia_245_guilds[refugia_245_guilds == 1] <- 0
refugia_245_guilds[refugia_245_guilds == 2] <- 0
refugia_245_guilds[refugia_245_guilds == 3] <- 1

refugia_585_guilds[refugia_585_guilds > 0] <- 1
refugia_585_guilds <- sum(refugia_585_guilds)
refugia_585_guilds[refugia_585_guilds == 1] <- 0
refugia_585_guilds[refugia_585_guilds == 2] <- 0
refugia_585_guilds[refugia_585_guilds == 3] <- 1

# could then check what percentage of this region is within contemporary AESs
refugia_245_guilds_mask <- mask(refugia_245_guilds, as_Spatial(contemporary_aes))
refugia_585_guilds_mask <- mask(refugia_585_guilds, as_Spatial(contemporary_aes))

sum(getValues(refugia_245_guilds_mask), na.rm = T) / sum(getValues(refugia_245_guilds), na.rm = T)
sum(getValues(refugia_585_guilds_mask), na.rm = T) / sum(getValues(refugia_585_guilds), na.rm = T)


p1 <- ggplot() + 
  theme_bw() +
  geom_raster(aes(x = x, y = y, fill = layer), data = refugia_245_guilds |> rasterToPoints() |> as_tibble()) +
  geom_sf(aes(), data = contemporary_aes, fill = NA)

p2 <- ggplot() + 
  theme_bw() +
  geom_raster(aes(x = x, y = y, fill = layer), data = refugia_585_guilds |> rasterToPoints() |> as_tibble()) +
  geom_sf(aes(), data = contemporary_aes, fill = NA)

p1 + p2


# ends
