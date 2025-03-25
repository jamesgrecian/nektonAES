##################
### futures_df ###
##################

# Helper function to generate plotting dataframe for potential futures
# input is simply species number
# output is a tibble that can be passed to geom_raster

# 2024-07-15

futures_df <- function(species = 1){
  
  # list available files
  files_now <- list.files("data/ensemble model outputs", pattern = "ensemble", full.names = T)
  files_now <- files_now[grep("v2", files_now)]
  files_245 <- list.files("data/ensemble model outputs", pattern = "ssp245", full.names = T)
  files_585 <- list.files("data/ensemble model outputs", pattern = "ssp585", full.names = T)
  
  # load files for focal species
  dat_now <- readRDS(files_now[species])
  dat_245 <- readRDS(files_245[species])
  dat_585 <- readRDS(files_585[species])
  
  # calculate estimated change at end of century
  dat_245_change <- dat_245 - dat_now
  dat_585_change <- dat_585 - dat_now
  names(dat_245_change) <- names(dat_245)
  names(dat_585_change) <- names(dat_585)
  
  # ensemble mean change
  ensemble_mean_245 <- mean(dat_245_change)
  ensemble_mean_585 <- mean(dat_585_change)
  
  # transform to tibbles for use with geom_raster
  dat_245_change_df <- dat_245_change |>
    rasterToPoints() |>
    as_tibble() |>
    pivot_longer(3:10, names_to = "model", values_to = "preds") |>
    mutate(scenario = "ssp245") |>
    dplyr::select(x, y, scenario, model, preds)
  
  dat_585_change_df <- dat_585_change |>
    rasterToPoints() |>
    as_tibble() |>
    pivot_longer(3:10, names_to = "model", values_to = "preds") |>
    mutate(scenario = "ssp585") |>
    dplyr::select(x, y, scenario, model, preds)
  
  ensemble_mean_245_df <- ensemble_mean_245 |>
    rasterToPoints() |>
    as_tibble() |>
    mutate(model = "ensemble",
           scenario = "ssp245") |>
    dplyr::select(x, y, scenario, model, layer) |>
    rename(preds = layer)
  
  
  ensemble_mean_585_df <- ensemble_mean_585 |>
    rasterToPoints() |>
    as_tibble() |>
    mutate(model = "ensemble",
           scenario = "ssp585") |>
    dplyr::select(x, y, scenario, model, layer) |>
    rename(preds = layer)
  
  future_tibble <- ensemble_mean_245_df |> 
    rbind(dat_245_change_df,
          ensemble_mean_585_df,
          dat_585_change_df)
  
  return(future_tibble)
}

# ends
