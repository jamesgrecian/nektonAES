#############################################################
### Fit four model ensembles to each species with weights ###
#############################################################

# 2024-06-24

# libraries
require(tidymodels)
require(sf)
sf_use_s2(FALSE)
require(spatialsample)
#require(future)

source("R/fit_ensemble.R")

# make the function a bit more straightforward to trouble shoot
# load the data into the environment
# pass data to function rather than loading within...
dat_all <- readRDS("data/presence_absence_data_10k_with_covariates_2024-06-17.rds")
dat_all <- dat_all |> mutate(PresAbs = factor(PresAbs))
dat_all <- dat_all |> mutate(cwts = hardhat::importance_weights(NA)) # create dummy weights column to pass to recipes later

# load pre-generated cv folds
folds_all <- readRDS("data/folds_weights_v2.rds")

# load results from species specific covariate selection
tidy_results_weights_all <- readRDS("data/tidy_results_weights_v2.rds")

#####################
### Fit ensembles ### 
#####################

# now for each species define which section of the data need to be used:
  
# 1. Bathylagus antarcticus
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[1]),
                          folds = folds_all[[1]],
                          term_table = tidy_results_weights_all[[1]])
saveRDS(model_out, "data/ensemble model outputs/model_df_01_v2.rds")
rm(model_out)

# 2. Notolepis coatsorum
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[2]),
                          folds = folds_all[[2]],
                          term_table = tidy_results_weights_all[[2]])
saveRDS(model_out, "data/ensemble model outputs/model_df_02_v2.rds")
rm(model_out)

# 3. Pleuragramma antarctica  
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[3]),
                          folds = folds_all[[3]],
                          term_table = tidy_results_weights_all[[3]])
saveRDS(model_out, "data/ensemble model outputs/model_df_03_v2.rds")
rm(model_out)

# 4. Electrona antarctica     
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[4]),
                          folds = folds_all[[4]],
                          term_table = tidy_results_weights_all[[4]])
saveRDS(model_out, "data/ensemble model outputs/model_df_04_v2.rds")
rm(model_out)

# 5. Electrona carlsbergi
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[5]),
                          folds = folds_all[[5]],
                          term_table = tidy_results_weights_all[[5]])
saveRDS(model_out, "data/ensemble model outputs/model_df_05_v2.rds")
rm(model_out)

# 6. Gymnoscopelus bolini
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[6]),
                          folds = folds_all[[6]],
                          term_table = tidy_results_weights_all[[6]])
saveRDS(model_out, "data/ensemble model outputs/model_df_06_v2.rds")
rm(model_out)

# 7.Gymnoscopelus braueri
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[7]),
                          folds = folds_all[[7]],
                          term_table = tidy_results_weights_all[[7]])
saveRDS(model_out, "data/ensemble model outputs/model_df_07_v2.rds")
rm(model_out)

# 8. Gymnoscopelus fraseri 
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[8]),
                          folds = folds_all[[8]],
                          term_table = tidy_results_weights_all[[8]])
saveRDS(model_out, "data/ensemble model outputs/model_df_08_v2.rds")
rm(model_out)

# 9. Gymnoscopelus nicholsi
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[9]),
                          folds = folds_all[[9]],
                          term_table = tidy_results_weights_all[[9]])
saveRDS(model_out, "data/ensemble model outputs/model_df_09_v2.rds")
rm(model_out)

#10. Gymnoscopelus opisthopterus
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[10]),
                          folds = folds_all[[10]],
                          term_table = tidy_results_weights_all[[10]])
saveRDS(model_out, "data/ensemble model outputs/model_df_10_v2.rds")
rm(model_out)

#11. Krefftichthys anderssoni
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[11]),
                          folds = folds_all[[11]],
                          term_table = tidy_results_weights_all[[11]])
saveRDS(model_out, "data/ensemble model outputs/model_df_11_v2.rds")
rm(model_out)

# 12. Protomyctophum bolini
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[12]),
                          folds = folds_all[[12]],
                          term_table = tidy_results_weights_all[[12]])
saveRDS(model_out, "data/ensemble model outputs/model_df_12_v2.rds")
rm(model_out)

# 13.Protomyctophum tenisoni
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[13]),
                          folds = folds_all[[13]],
                          term_table = tidy_results_weights_all[[13]])
saveRDS(model_out, "data/ensemble model outputs/model_df_13_v2.rds")
rm(model_out)

# 14. Euphausia superba
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[14]),
                          folds = folds_all[[14]],
                          term_table = tidy_results_weights_all[[14]])
saveRDS(model_out, "data/ensemble model outputs/model_df_14_v2.rds")
rm(model_out)

# 15. Alluroteuthis antarcticus
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[15]),
                          folds = folds_all[[15]],
                          term_table = tidy_results_weights_all[[15]])
saveRDS(model_out, "data/ensemble model outputs/model_df_15_v2.rds")
rm(model_out)

# 16. Bathyteuthis abyssicola
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[16]),
                          folds = folds_all[[16]],
                          term_table = tidy_results_weights_all[[16]])
saveRDS(model_out, "data/ensemble model outputs/model_df_16_v2.rds")
rm(model_out)

# 17. Galiteuthis glacialis
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[17]),
                          folds = folds_all[[17]],
                          term_table = tidy_results_weights_all[[17]])
saveRDS(model_out, "data/ensemble model outputs/model_df_17_v2.rds")
rm(model_out)

# 18. Gonatus antarcticus
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[18]),
                          folds = folds_all[[18]],
                          term_table = tidy_results_weights_all[[18]])
saveRDS(model_out, "data/ensemble model outputs/model_df_18_v2.rds")
rm(model_out)

# 19. Histioteuthis atlantica
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[19]),
                          folds = folds_all[[19]],
                          term_table = tidy_results_weights_all[[19]])
saveRDS(model_out, "data/ensemble model outputs/model_df_19_v2.rds")
rm(model_out)

# 20. Histioteuthis eltaninae
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[20]),
                          folds = folds_all[[20]],
                          term_table = tidy_results_weights_all[[20]])
saveRDS(model_out, "data/ensemble model outputs/model_df_20_v2.rds")
rm(model_out)

# 21. Kondakovia longimana
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[21]),
                          folds = folds_all[[21]],
                          term_table = tidy_results_weights_all[[21]])
saveRDS(model_out, "data/ensemble model outputs/model_df_21_v2.rds")
rm(model_out)

# 22. Martialia hyadesi
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[22]),
                          folds = folds_all[[22]],
                          term_table = tidy_results_weights_all[[22]])
saveRDS(model_out, "data/ensemble model outputs/model_df_22_v2.rds")
rm(model_out)

# 23. Mesonychoteuthis hamiltoni
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[23]),
                          folds = folds_all[[23]],
                          term_table = tidy_results_weights_all[[23]])
saveRDS(model_out, "data/ensemble model outputs/model_df_23_v2.rds")
rm(model_out)

# 24. Moroteuthis ingens 
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[24]),
                          folds = folds_all[[24]],
                          term_table = tidy_results_weights_all[[24]])
saveRDS(model_out, "data/ensemble model outputs/model_df_24_v2.rds")
rm(model_out)

# 25. Moroteuthis robsoni
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[25]),
                          folds = folds_all[[25]],
                          term_table = tidy_results_weights_all[[25]])
saveRDS(model_out, "data/ensemble model outputs/model_df_25_v2.rds")
rm(model_out)

# 26. Psychroteuthis glacialis
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[26]),
                          folds = folds_all[[26]],
                          term_table = tidy_results_weights_all[[26]])
saveRDS(model_out, "data/ensemble model outputs/model_df_26_v2.rds")
rm(model_out)

# 27. Slosarczykovia circumantarctica
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[27]),
                          folds = folds_all[[27]],
                          term_table = tidy_results_weights_all[[27]])
saveRDS(model_out, "data/ensemble model outputs/model_df_27_v2.rds")
rm(model_out)

# 28. Todarodes filippovae
model_out <- fit_ensemble(dat = dat_all |> filter(species == unique(species)[28]),
                          folds = folds_all[[28]],
                          term_table = tidy_results_weights_all[[28]])
saveRDS(model_out, "data/ensemble model outputs/model_df_28_v2.rds")
rm(model_out)

# ends
