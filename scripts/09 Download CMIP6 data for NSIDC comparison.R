############################################################
### Process CMIP6 sea ice data for comparison with NSIDC ###
############################################################

# 2024-06-25

# Need to select a subset of CMIP6 models to use for future climates

# Download a subset of 15 models mentioned in Roach et al. 2020
# constrained by those that have ssp245 and ssp585
# also constrained by those that have all required covariate data
# rank these against NSIDC using RMSE
# pick top 8 models to use in ensemble

# ensure first that this subset has data for each covariate and each experiment:
require(tidyverse)
require(rcmip6)

# set up CMIP6 query
query <- list(
  type           = "Dataset",
  source_id      = c("ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5",
                     "CESM2-WACCM", "CMCC-CM2-SR5", "EC-Earth3", "EC-Earth3-Veg",
                     "FGOALS-g3", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-LR",
                     "MRI-ESM2-0", "NESM3", "NorESM2-LM"),
  experiment_id  = c("historical", "ssp245", "ssp585"),
  variable_id    = c("tos", "sos", "zos", "siconc", "mlotst"),
  replica        = "false",
  latest         = "true",
  project        = "CMIP6",
  frequency      = "mon",
  member_id      = "r1i1p1f1",
  table_id       = c("Omon", "SImon"),
  grid_label     = "gn"
  )

# check for availability on portal
results <- cmip_search(query)

# check results
cmip_size(results) # 85Gb of data to download...

# output summary table to check results
df <- cmip_simplify(results) |> tibble() |> dplyr::select(source_id, experiment_id, member_id, variable_id, grid_label, nominal_resolution)
df <- df |> pivot_wider(names_from = variable_id, values_from = variable_id, values_fill = list(variable_id = NA))
df |> arrange(source_id, experiment_id) |> print(n = Inf)

# These 15 have all the required data...

# Now download the sea ice data for each model:
query <- list(
  type           = "Dataset",
  source_id      = c("ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5",
                     "CESM2-WACCM", "CMCC-CM2-SR5", "EC-Earth3", "EC-Earth3-Veg",
                     "FGOALS-g3", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-LR",
                     "MRI-ESM2-0", "NESM3", "NorESM2-LM"),
  experiment_id  = "historical",
  variable_id    = "siconc",
  replica        = "false",
  latest         = "true",
  project        = "CMIP6",
  frequency      = "mon",
  member_id      = "r1i1p1f1",
  table_id       = "SImon",
  grid_label     = "gn"
  )

# check for availability on portal
results <- cmip_search(query)

# check results
cmip_size(results)

# output summary table to check results
df <- cmip_simplify(results) |> tibble() |> dplyr::select(source_id, experiment_id, member_id, variable_id, grid_label, nominal_resolution)
df <- df |> pivot_wider(names_from = variable_id, values_from = variable_id, values_fill = list(variable_id = NA))
df |> arrange(source_id, experiment_id) |> print(n = Inf)

# download files required
cmip_root_set("/Users/home/ANTSIE data") # specify root directory
options(timeout = 1200)            # downloading from server so specify how long to keep connection open (will fail if not defined)
cmip_download(results)

# Pass the downloaded files to the cdo bash script to process from native grid to regular 1x1 degree grid

# ends


