############################
### Ensemble predictions ###
############################

# 2024-04-29

# a tidy version of the ensemble prediction script
# now wrapped as a function

# following on from ANTSIE script:
# new format for outputed models and AUCs
# using the saved model objects from the fit model ensemble script
# generate predictions for partial dependence plots
# generate predictions of spatial distributions across the southern ocean

# load libraries
require(tidymodels)
require(sf)
sf::sf_use_s2(FALSE)
require(spatialsample)
require(DALEXtra)
require(orsifronts)

# Mikes hack function
quantile.hardhat_importance_weights <- \(x, ...) rep(NA, length(x))

# load data
dat_all <- readRDS("data/presence_absence_data_10k_with_covariates_2024-06-17.rds")
folds_all <- readRDS("data/folds_weights_v2.rds")

################################
### Predict for each species ###
################################

source("R/predict_ensemble.R")

# 1. Bathylagus antarcticus
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_01_v2.rds"),
                 folds = folds_all[[1]],
                 dat = dat_all |> filter(species == unique(species)[1]),
                 plot_nm = "figures/species 01 ensemble plot_2024-06-21.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_01_v2.rds")

# 2. Notolepis coatsorum
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_02_v2.rds"),
                 folds = folds_all[[2]],
                 dat = dat_all |> filter(species == unique(species)[2]),
                 plot_nm = "figures/species 02 ensemble plot 2024-06-21.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_02_v2.rds")

# 3. Pleuragramma antarctica  
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_03_v2.rds"),
                 folds = folds_all[[3]],
                 dat = dat_all |> filter(species == unique(species)[3]),
                 plot_nm = "figures/species 03 ensemble plot 2024-06-21.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_03_v2.rds")

# 4. Electrona antarctica     
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_04_v2.rds"),
                 folds = folds_all[[4]],
                 dat = dat_all |> filter(species == unique(species)[4]),
                 plot_nm = "figures/species 04 ensemble plot 2024-06-21.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_04_v2.rds")

# 5. Electrona carlsbergi
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_05_v2.rds"),
                 folds = folds_all[[5]],
                 dat = dat_all |> filter(species == unique(species)[5]),
                 plot_nm = "figures/species 05 ensemble plot 2024-06-21.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_05_v2.rds")

# 6. Gymnoscopelus bolini
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_06_v2.rds"),
                 folds = folds_all[[6]],
                 dat = dat_all |> filter(species == unique(species)[6]),
                 plot_nm = "figures/species 06 ensemble plot 2024-06-21.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_06_v2.rds")




# 7.Gymnoscopelus braueri
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_07_v2.rds"),
                 folds = folds_all[[7]],
                 dat = dat_all |> filter(species == unique(species)[7]),
                 plot_nm = "figures/species 07 ensemble plot 2024-06-21.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_07_v2.rds")

# 8. Gymnoscopelus fraseri 
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_08_v2.rds"),
                 folds = folds_all[[8]],
                 dat = dat_all |> filter(species == unique(species)[8]),
                 plot_nm = "figures/species 08 ensemble plot 2024-06-21.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_08_v2.rds")

# 9. Gymnoscopelus nicholsi
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_09_v2.rds"),
                 folds = folds_all[[9]],
                 dat = dat_all |> filter(species == unique(species)[9]),
                 plot_nm = "figures/species 09 ensemble plot 2024-06-21.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_09_v2.rds")

# 10. Gymnoscopelus opisthopterus
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_10_v2.rds"),
                 folds = folds_all[[10]],
                 dat = dat_all |> filter(species == unique(species)[10]),
                 plot_nm = "figures/species 10 ensemble plot 2024-06-21.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_10_v2.rds")

# 11. Krefftichthys anderssoni
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_11_v2.rds"),
                 folds = folds_all[[11]],
                 dat = dat_all |> filter(species == unique(species)[11]),
                 plot_nm = "figures/species 11 ensemble plot 2024-06-21.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_11_v2.rds")

# 12. Protomyctophum bolini
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_12_v2.rds"),
                 folds = folds_all[[12]],
                 dat = dat_all |> filter(species == unique(species)[12]),
                 plot_nm = "figures/species 12 ensemble plot 2024-06-21.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_12_v2.rds")








# 13.Protomyctophum tenisoni
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_13.rds"),
                 folds = folds_all[[13]],
                 dat = dat_all |> filter(species == unique(species)[13]),
                 plot_nm = "figures/species 13 ensemble plot.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_13.rds")

# 14. Euphausia superba
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_14.rds"),
                 folds = folds_all[[14]],
                 dat = dat_all |> filter(species == unique(species)[14]),
                 plot_nm = "figures/species 14 ensemble plot.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_14.rds")

# 15. Alluroteuthis antarcticus
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_15.rds"),
                 folds = folds_all[[15]],
                 dat = dat_all |> filter(species == unique(species)[15]),
                 plot_nm = "figures/species 15 ensemble plot.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_15.rds")

# 16. Bathyteuthis abyssicola
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_16.rds"),
                 folds = folds_all[[16]],
                 dat = dat_all |> filter(species == unique(species)[16]),
                 plot_nm = "figures/species 16 ensemble plot.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_16.rds")

# 17. Galiteuthis glacialis
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_17.rds"),
                 folds = folds_all[[17]],
                 dat = dat_all |> filter(species == unique(species)[17]),
                 plot_nm = "figures/species 17 ensemble plot.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_17.rds")

# 18. Gonatus antarcticus
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_18.rds"),
                 folds = folds_all[[18]],
                 dat = dat_all |> filter(species == unique(species)[18]),
                 plot_nm = "figures/species 18 ensemble plot.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_18.rds")

# 19. Histioteuthis atlantica
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_19.rds"),
                 folds = folds_all[[19]],
                 dat = dat_all |> filter(species == unique(species)[19]),
                 plot_nm = "figures/species 19 ensemble plot.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_19.rds")

# 20. Histioteuthis eltaninae
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_20.rds"),
                 folds = folds_all[[20]],
                 dat = dat_all |> filter(species == unique(species)[20]),
                 plot_nm = "figures/species 20 ensemble plot.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_20.rds")

# 21. Kondakovia longimana
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_21.rds"),
                 folds = folds_all[[21]],
                 dat = dat_all |> filter(species == unique(species)[21]),
                 plot_nm = "figures/species 21 ensemble plot.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_21.rds")

# 22. Martialia hyadesi
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_22.rds"),
                 folds = folds_all[[22]],
                 dat = dat_all |> filter(species == unique(species)[22]),
                 plot_nm = "figures/species 22 ensemble plot.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_22.rds")

# 23. Mesonychoteuthis hamiltoni
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_23.rds"),
                 folds = folds_all[[23]],
                 dat = dat_all |> filter(species == unique(species)[23]),
                 plot_nm = "figures/species 23 ensemble plot.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_23.rds")

# 24. Moroteuthis ingens 
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_24.rds"),
                 folds = folds_all[[24]],
                 dat = dat_all |> filter(species == unique(species)[24]),
                 plot_nm = "figures/species 24 ensemble plot.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_24.rds")

# 25. Moroteuthis robsoni
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_25.rds"),
                 folds = folds_all[[25]],
                 dat = dat_all |> filter(species == unique(species)[25]),
                 plot_nm = "figures/species 25 ensemble plot.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_25.rds")

# 26. Psychroteuthis glacialis
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_26.rds"),
                 folds = folds_all[[26]],
                 dat = dat_all |> filter(species == unique(species)[26]),
                 plot_nm = "figures/species 26 ensemble plot.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_26.rds")

# 27. Slosarczykovia circumantarctica
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_27.rds"),
                 folds = folds_all[[27]],
                 dat = dat_all |> filter(species == unique(species)[27]),
                 plot_nm = "figures/species 27 ensemble plot.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_27.rds")

# 28. Todarodes filippovae
predict_ensemble(model_df = readRDS("data/ensemble model outputs/model_df_28.rds"),
                 folds = folds_all[[28]],
                 dat = dat_all |> filter(species == unique(species)[28]),
                 plot_nm = "figures/species 28 ensemble plot.jpg",
                 file_nm = "data/ensemble model outputs/ensemble_raster_28.rds")

# ends

# old themes
# incl. setting everything black...
#
# en_black <- en + theme(plot.background = element_rect(fill = "black"),
#                        panel.background = element_rect(fill = 'black'),
#                        legend.title = element_text(color = "white"),
#                        legend.text = element_text(color = "white"))
# 
# enst_black <- enst + theme(plot.background = element_rect(fill = "black"),
#                            panel.background = element_rect(fill = 'black'),
#                            legend.title = element_text(color = "white"),
#                            legend.text = element_text(color = "white"),
#                            strip.text = element_text(colour = 'white'))
# 
# p_grid_black <- p_grid + 
#   theme(plot.background = element_rect(fill = "black"),
#         panel.background = element_rect(fill = 'black'),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.title = element_text(color = "white"),
#         legend.text = element_text(color = "white"))
# 
# quartz(width = 11, height = 6.5)
# p_grid_black
# quartz.save(file = "PDP grid black.jpg",
#             type = "jpeg",
#             dev = dev.cur(),
#             dpi = 500)
# 
# # a wrapper function so you don't need 10 chunks of explainer scripts
# explain_wrapper <- function(input_model, input_data){
#   pred <- input_data |> analysis() |> st_drop_geometry() |> dplyr::select(c(PresAbs, strsplit(tms, "\\s\\+\\s", perl = TRUE)[[1]])) |>
#     mutate(cwts = hardhat::importance_weights(ifelse(PresAbs == 1, 1, sum(PresAbs == 1) / sum(PresAbs == 0))))
#   PresAbs <- input_data |> analysis() |> st_drop_geometry() |> dplyr::select(PresAbs)
#   pred_out <- explain_tidymodels(input_model, data = pred, y = PresAbs) |> model_profile(N = 1000, type = "partial")
#   pred_out <- pred_out[[2]] |> as_tibble()
#   return(pred_out)
# }
# 
