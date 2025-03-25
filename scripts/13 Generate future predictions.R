####################################################################
### Generate end of century spatial predictions for each species ###
####################################################################

# 2024-07-08

# There are 28 species
# There are 8 GCMs
# There are 2 future scenarios
# Each ensemble model is 4 models x 10 folds
# (4 x 10) x (2 x 8) = 640 surfaces per species...

# Calculate future for each fold and each model
# Don't present individual folds, average over all
# Also average over ESDM models
# Then each esdm is only one field
# Generate that field for each GCM and each scenario
# Save a raster stack with 2 scenarios x 8 GCMs

# load libraries
require(tidyverse)
require(raster)
require(tidymodels)
require(sf)
sf::sf_use_s2(FALSE)
require(spatialsample)
require(DALEXtra)
require(orsifronts)

source("R/predict_futures.R")

# Mikes hack function
quantile.hardhat_importance_weights <- \(x, ...) rep(NA, length(x))

# load folds
folds_all <- readRDS("data/folds_weights_v2.rds")

# use custom function to generate a prediction stack for each SSP
# containing estimated distributions for each of the 8 GCMs from ensemble
# don't store the individual folds or individual SDM models
predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_01_v2.rds"),
                folds = folds_all[[1]],
                file_path_1 = "data/ensemble model outputs/futures_raster_01_ssp245_v2.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_01_ssp585_v2.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_02_v2.rds"),
                folds = folds_all[[2]],
                file_path_1 = "data/ensemble model outputs/futures_raster_02_ssp245_v2.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_02_ssp585_v2.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_03_v2.rds"),
                folds = folds_all[[3]],
                file_path_1 = "data/ensemble model outputs/futures_raster_03_ssp245_v2.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_03_ssp585_v2.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_04_v2.rds"),
                folds = folds_all[[4]],
                file_path_1 = "data/ensemble model outputs/futures_raster_04_ssp245_v2.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_04_ssp585_v2.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_05_v2.rds"),
                folds = folds_all[[5]],
                file_path_1 = "data/ensemble model outputs/futures_raster_05_ssp245_v2.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_05_ssp585_v2.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_06_v2.rds"),
                folds = folds_all[[6]],
                file_path_1 = "data/ensemble model outputs/futures_raster_06_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_06_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_07_v2.rds"),
                folds = folds_all[[7]],
                file_path_1 = "data/ensemble model outputs/futures_raster_07_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_07_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_08_v2.rds"),
                folds = folds_all[[8]],
                file_path_1 = "data/ensemble model outputs/futures_raster_08_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_08_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_09_v2.rds"),
                folds = folds_all[[9]],
                file_path_1 = "data/ensemble model outputs/futures_raster_09_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_09_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_10_v2.rds"),
                folds = folds_all[[10]],
                file_path_1 = "data/ensemble model outputs/futures_raster_10_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_10_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_11_v2.rds"),
                folds = folds_all[[11]],
                file_path_1 = "data/ensemble model outputs/futures_raster_11_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_11_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_12_v2.rds"),
                folds = folds_all[[12]],
                file_path_1 = "data/ensemble model outputs/futures_raster_12_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_12_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_13_v2.rds"),
                folds = folds_all[[13]],
                file_path_1 = "data/ensemble model outputs/futures_raster_13_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_13_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_14_v2.rds"),
                folds = folds_all[[14]],
                file_path_1 = "data/ensemble model outputs/futures_raster_14_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_14_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_15_v2.rds"),
                folds = folds_all[[15]],
                file_path_1 = "data/ensemble model outputs/futures_raster_15_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_15_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_16_v2.rds"),
                folds = folds_all[[16]],
                file_path_1 = "data/ensemble model outputs/futures_raster_16_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_16_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_17_v2.rds"),
                folds = folds_all[[17]],
                file_path_1 = "data/ensemble model outputs/futures_raster_17_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_17_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_18_v2.rds"),
                folds = folds_all[[18]],
                file_path_1 = "data/ensemble model outputs/futures_raster_18_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_18_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_19_v2.rds"),
                folds = folds_all[[19]],
                file_path_1 = "data/ensemble model outputs/futures_raster_19_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_19_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_20_v2.rds"),
                folds = folds_all[[20]],
                file_path_1 = "data/ensemble model outputs/futures_raster_20_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_20_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_21_v2.rds"),
                folds = folds_all[[21]],
                file_path_1 = "data/ensemble model outputs/futures_raster_21_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_21_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_22_v2.rds"),
                folds = folds_all[[22]],
                file_path_1 = "data/ensemble model outputs/futures_raster_22_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_22_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_23_v2.rds"),
                folds = folds_all[[23]],
                file_path_1 = "data/ensemble model outputs/futures_raster_23_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_23_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_24_v2.rds"),
                folds = folds_all[[24]],
                file_path_1 = "data/ensemble model outputs/futures_raster_24_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_24_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_25_v2.rds"),
                folds = folds_all[[25]],
                file_path_1 = "data/ensemble model outputs/futures_raster_25_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_25_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_26_v2.rds"),
                folds = folds_all[[26]],
                file_path_1 = "data/ensemble model outputs/futures_raster_26_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_26_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_27_v2.rds"),
                folds = folds_all[[27]],
                file_path_1 = "data/ensemble model outputs/futures_raster_27_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_27_ssp585.rds")

predict_futures(model_df = readRDS("data/ensemble model outputs/model_df_28_v2.rds"),
                folds = folds_all[[28]],
                file_path_1 = "data/ensemble model outputs/futures_raster_28_ssp245.rds",
                file_path_2 = "data/ensemble model outputs/futures_raster_28_ssp585.rds")

# ends
