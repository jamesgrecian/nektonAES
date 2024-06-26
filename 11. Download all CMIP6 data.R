
# 8 GCM subset is...
# ACCESS.CM2
# BCC.CSM2.MR
# CESM2.WACCM
# CMCC.CM2.SR5
# FGOALS.G3
# IPSL.CM6A.LR
# MRI_ESM2.0
# NESM3

# ensure first that this subset has data for each covariate and each experiment:
require(tidyverse)
require(rcmip6)

# set up CMIP6 query
query <- list(
  type           = "Dataset",
  source_id      = c("ACCESS-CM2", "BCC-CSM2-MR", "CESM2-WACCM", "CMCC-CM2-SR5",
                     "FGOALS-g3", "IPSL-CM6A-LR", "MRI-ESM2-0", "NESM3"),
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
