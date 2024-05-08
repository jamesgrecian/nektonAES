#################################################
### Use spatialsample package to assign folds ###
#################################################

# 2023-04-10
# modified 2023-04-17 to make PresAbs a factor
# modified 2023-05-29 to include weights

# More info on package here: https://spatialsample.tidymodels.org/index.html
# spatialsample::spatial_block_cv function seems to work similarlily to block_CV
# however the function has fewer issues partitioning the data
# this includes being able to assign k = 10 folds rather than k = 6 folds
# literature suggests k = 10 folds is better...
# for now focus on 500 km block size - this could be varied later
# also question over whether block size should vary per species... adds complexity
# blockCV paper suggests matching block size to spatial autocorrelation...?

# spatal_block_CV will not partition data correctly on first go
# sometimes one block will contain no presences
# this will cause issues later down the line - can't get AUC if there are no presences to test against
# solution is to run while loop and only stop when each block has `some` presences in it...

# load libraries
require(tidyverse)
require(spatialsample)
require(sf)

# define projection
prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# load data
dat <- readRDS("data/presence_absence_data_10k_with_covariates_2024-04-25.rds")

# make PresAbs a factor
dat <- dat |> mutate(PresAbs = factor(PresAbs))

# create dummy weights column to pass to recipes later
dat <- dat |> mutate(cwts = hardhat::importance_weights(NA))

# make sf object
dat_sf <- dat |> st_as_sf(coords = c("x", "y"), crs = prj)

# initialise list to store folds
folds <- NULL 

# need to ensure presence and absence points are in each fold
# create the fold, then tally presence and absence points
# make sure it's summing to 20
# when PresAbs is a factor values are stored in check as 1 2 not 0 1
# -1 to revert to 0 1 so that summing to 20 still works - not summing to 60!

for(k in c(1:28)){ 
  print(k)
  check <- NULL
  while(sum(check-1) < 20){ # 10 folds and need presences in both analysis and assessment groups 10 * 2 = 20
    check <- NULL
    sub_df <- dat_sf %>% filter(species == unique(species)[k])
    block_folds <- spatial_block_cv(sub_df, v = 10, square = F, cellsize = 50000, method = "snake")
    for (j in 1:10){
      check <- c(check,
                 block_folds$splits[[j]] %>% analysis %>% pull(PresAbs) %>% unique,
                 block_folds$splits[[j]] %>% assessment %>% pull(PresAbs) %>% unique)
    }
  }
  folds[[k]] <- block_folds
}

# save folds object
saveRDS(folds, "data/folds_weights_50_snake.rds")

# ends