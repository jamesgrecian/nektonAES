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

# function operates the same as the general prediction function
# but no plotting
# input is:
# the model_df to scrape which covariates are needed for each species
# folds subset for focal species
# file path and name to output the two SSP raster prediction stacks

predict_futures <- function(model_df, folds, file_path_1, file_path_2){
  
  # identify which covariates required using fitted model
  model_covs <- names(model_df$workflows[[1]]$pre$mold$predictors)
  
  # custom helper function for change factor protocol
  source("R/change_factor.R")
  
  # load covariate stack of observations
  covariates <- readRDS("data/covariate_stack.rds")
  
  # GCMs
  model_list <- c("ACCESS-CM2", "BCC-CSM2-MR", "CESM2-WACCM", "CMCC-CM2-SR5",
                  "FGOALS-g3", "IPSL-CM6A-LR", "MRI-ESM2-0", "NESM3")
  
  ### at this point loop through model list? ##
  
  output_raster_245 <- stack()
  output_raster_585 <- stack()
  
  for(i in 1:length(model_list)){
    cat(paste0("\nprocessing ", model_list[i], "..."))
    
    # change factor protocol function
    covs_245 <- change_factor(covariate_stack = covariates,
                              model = model_list[i],
                              scenario = "ssp245")
    covs_585 <- change_factor(covariate_stack = covariates,
                              model = model_list[i],
                              scenario = "ssp585")
    covs_245 <- subset(covs_245, model_covs)
    covs_585 <- subset(covs_585, model_covs)
    
    ####################################
    ### Generate spatial predictions ###
    ####################################
    
    # generate predictions
    # https://github.com/tidymodels/planning/issues/26
    
    ### now need to generate 10 spatial predictions
    ### could simply pass the covariate data to the predict function
    ### or should this be some other weird explainer thing?
    explainer_generator <- function(input_model, input_data) {
      pred <-
        input_data |> analysis() |> st_drop_geometry() |> dplyr::select(all_of(model_covs))
      resp <-
        input_data |> analysis() |> st_drop_geometry() |> pull(PresAbs)
      explainer <-
        explain_tidymodels(input_model, data = pred, y = resp, verbose = F)
      return(explainer)
    }
    
    # create explainers for each model type
    cat("\n...processing glm explainers")
    glm_explainers <- map2(model_df$workflows[1:10], folds$splits, explainer_generator)
    cat("\n...processing gam explainers")
    gam_explainers <- map2(model_df$workflows[11:20], folds$splits, explainer_generator)
    cat("\n...processing rf explainers")
    rf_explainers <- map2(model_df$workflows[21:30], folds$splits, explainer_generator)
    cat("\n...processing brt explainers")
    brt_explainers <- map2(model_df$workflows[31:40], folds$splits, explainer_generator)
    
    # generate predictions for covariate stack
    cat("\n...generating ensemble predictions")
    glm_preds_stack_245 <- lapply(glm_explainers, function(x) terra::predict(covs_245, x)) |> raster::stack()
    gam_preds_stack_245 <- lapply(gam_explainers, function(x) terra::predict(covs_245, x)) |> raster::stack()
    rf_preds_stack_245 <- lapply(rf_explainers, function(x) terra::predict(covs_245, x)) |> raster::stack()
    brt_preds_stack_245 <- lapply(brt_explainers, function(x) terra::predict(covs_245, x)) |> raster::stack()
    
    glm_preds_stack_585 <- lapply(glm_explainers, function(x) terra::predict(covs_585, x)) |> raster::stack()
    gam_preds_stack_585 <- lapply(gam_explainers, function(x) terra::predict(covs_585, x)) |> raster::stack()
    rf_preds_stack_585 <- lapply(rf_explainers, function(x) terra::predict(covs_585, x)) |> raster::stack()
    brt_preds_stack_585 <- lapply(brt_explainers, function(x) terra::predict(covs_585, x)) |> raster::stack()
    
    # need final workflows from each model saved so that weights can be calculated...
    glm_auc_weights <- model_df |> filter(model == "GLM") |> pull(AUC)
    gam_auc_weights <- model_df |> filter(model == "GAM") |> pull(AUC)
    rf_auc_weights <- model_df |> filter(model == "RF") |> pull(AUC)
    brt_auc_weights <- model_df |> filter(model == "BRT") |> pull(AUC)
    
    # Create AUC weighted average of each model for the ensemble
    glm_preds_raster_245 <- weighted.mean(glm_preds_stack_245, w = glm_auc_weights)
    gam_preds_raster_245 <- weighted.mean(gam_preds_stack_245, w = gam_auc_weights)
    rf_preds_raster_245 <- weighted.mean(rf_preds_stack_245, w = rf_auc_weights)
    brt_preds_raster_245 <- weighted.mean(brt_preds_stack_245, w = brt_auc_weights)
    
    glm_preds_raster_585 <- weighted.mean(glm_preds_stack_585, w = glm_auc_weights)
    gam_preds_raster_585 <- weighted.mean(gam_preds_stack_585, w = gam_auc_weights)
    rf_preds_raster_585 <- weighted.mean(rf_preds_stack_585, w = rf_auc_weights)
    brt_preds_raster_585 <- weighted.mean(brt_preds_stack_585, w = brt_auc_weights)
    
    # ensemble mean
    ensemble_stack_245 <- stack(glm_preds_raster_245,
                                gam_preds_raster_245,
                                rf_preds_raster_245,
                                brt_preds_raster_245)
    
    ensemble_stack_585 <- stack(glm_preds_raster_585,
                                gam_preds_raster_585,
                                rf_preds_raster_585,
                                brt_preds_raster_585)
    
    auc_w <- c(mean(glm_auc_weights),
               mean(gam_auc_weights),
               mean(rf_auc_weights),
               mean(brt_auc_weights))
    
    ensemble_raster_245 <- weighted.mean(ensemble_stack_245, w = auc_w)
    ensemble_raster_585 <- weighted.mean(ensemble_stack_585, w = auc_w)
    
    names(ensemble_raster_245) <- model_list[i]
    names(ensemble_raster_585) <- model_list[i]
    
    output_raster_245 <- stack(output_raster_245, ensemble_raster_245)
    output_raster_585 <- stack(output_raster_585, ensemble_raster_585)
  }
  
  saveRDS(output_raster_245, file_path_1)
  saveRDS(output_raster_585, file_path_2)
  
}
