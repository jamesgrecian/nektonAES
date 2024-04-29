############################
### Ensemble predictions ###
############################

# 2024-04-29

# custom function to generate predictions from fitted models
# input arguments are:
# model_df - the fitted ensemble model workflow dataframe
# folds - the spatialsample folds for the focal species
# dat - the data frame for the focal specie - for plotting only
# plot_nm - filename for output plot
# file_nm - filename for output ensemble raster

predict_ensemble <- function(model_df, folds, dat, plot_nm, file_nm) {
  
  #######################
  ### Prep input data ###
  #######################
  
  # identify covariates from fitted model
  model_covs <- names(model_df$workflows[[1]]$pre$mold$predictors)
  
  # load in original data set for tick marks
  points <- dat |> dplyr::select(PresAbs, x, y) |> filter(PresAbs == 1)
  dat <- dat |> dplyr::select(c(PresAbs, model_covs))  # covariates of interest
  
  # don't need to plot all 10k points to illustrate distribution!
  # should be all the presence points and same number of absence points...
  dat_1 <- dat |> filter(PresAbs == 1)
  dat_0 <- dat |> filter(PresAbs == 0) |> slice_sample(n = nrow(dat_1))
  dat <- rbind(dat_1, dat_0)
  
  ################################
  ### Partial Dependence Plots ###
  ################################
  
  # a wrapper function so you don't need 10 chunks of explainer scripts
  explain_wrapper <- function(input_model, input_data) {
    pred <-
      input_data |> analysis() |> st_drop_geometry() |> dplyr::select(c(PresAbs, model_covs)) |>
      mutate(cwts = hardhat::importance_weights(ifelse(
        PresAbs == 1, 1, sum(PresAbs == 1) / sum(PresAbs == 0)
      )))
    PresAbs <-
      input_data |> analysis() |> st_drop_geometry() |> dplyr::select(PresAbs)
    pred_out <-
      explain_tidymodels(input_model, data = pred, y = PresAbs) |> model_profile(N = 1000, type = "partial")
    pred_out <- pred_out[[2]] |> as_tibble()
    return(pred_out)
  }
  
  glm_pdp_preds <- map2(model_df$workflows[1:10], folds$splits, explain_wrapper)
  gam_pdp_preds <- map2(model_df$workflows[11:20], folds$splits, explain_wrapper)
  rf_pdp_preds <-  map2(model_df$workflows[21:30], folds$splits, explain_wrapper)
  brt_pdp_preds <- map2(model_df$workflows[31:40], folds$splits, explain_wrapper)
  
  # format dataframe
  glm_pdp_preds <- glm_pdp_preds |> bind_rows(.id = "id")
  gam_pdp_preds <- gam_pdp_preds |> bind_rows(.id = "id")
  rf_pdp_preds <- rf_pdp_preds |> bind_rows(.id = "id")
  brt_pdp_preds <- brt_pdp_preds |> bind_rows(.id = "id")
  names(glm_pdp_preds) <- c("id", "var_name", "label", "x", "yhat", "ids")
  names(rf_pdp_preds) <- c("id", "var_name", "label", "x", "yhat", "ids")
  names(brt_pdp_preds) <- c("id", "var_name", "label", "x", "yhat", "ids")
  names(gam_pdp_preds) <- c("id", "var_name", "label", "x", "yhat", "ids")
  glm_pdp_preds$model <- "GLM"
  rf_pdp_preds$model <- "RF"
  brt_pdp_preds$model <- "BRT"
  gam_pdp_preds$model <- "GAM"
  
  # combine
  pdp_preds <- rbind(glm_pdp_preds, rf_pdp_preds, brt_pdp_preds, gam_pdp_preds)
  
  # refactor so order is intuitive in plot
  pdp_preds <- pdp_preds |>
    mutate(model = factor(model, levels = c("GLM", "GAM", "RF", "BRT")))
  
  # relabel to aid interpretation in plot
  pdp_preds <- pdp_preds |>
    mutate(
      var_name = recode_factor(
        var_name,
        sst = "sea surface temperature (ºC)",
        sst_grad = "sea surface temperature gradient",
        sal = "salinity (psu)",
        ssh = "sea surface height (m)",
        ssh_grad = "sea surface height gradient",
        mld = "mixed layer depth (m)",
        sic = "sea ice concentration (%)",
        bat = "bathymetry (m)"
      )
    )
  
  # format dat for plot tick marks
  dat <- dat |>
    pivot_longer(2:(length(model_covs) + 1),
                 names_to = "var_name",
                 values_to = "x") |>
    rename("yhat" = PresAbs)
  
  dat <- dat |>
    mutate(
      var_name = recode_factor(
        var_name,
        sst = "sea surface temperature (ºC)",
        sst_grad = "sea surface temperature gradient",
        sal = "salinity (psu)",
        ssh = "sea surface height (m)",
        ssh_grad = "sea surface height gradient",
        mld = "mixed layer depth (m)",
        sic = "sea ice concentration (%)",
        bat = "bathymetry (m)"
      )
    )
  
  facet_dat <- rbind(
    dat |> mutate(model = "GLM"),
    dat |> mutate(model = "GAM"),
    dat |> mutate(model = "RF"),
    dat |> mutate(model = "BRT")
  ) |>
    mutate(model = factor(model, levels = c("GLM", "GAM", "RF", "BRT")))
  
  # Partial dependence profile plots as a facet grid
  p_grid <- ggplot() +
    theme_bw(base_size = 10,
             base_family = "Helvetica Neue") +
    geom_line(
      aes(x = x, y = yhat, group = id),
      alpha = .5,
      linewidth = .5,
      data = pdp_preds
    ) +
    geom_point(
      aes(x = x, y = yhat),
      shape = "|",
      size = .5,
      alpha = .5,
      data = facet_dat
    ) +
    facet_grid(
      model ~ var_name,
      scales = "free_x",
      labeller = labeller(var_name = label_wrap_gen(
        width = 14, multi_line = TRUE
      ))
    ) + #wrap titles if they are too wide
    ylim(0, 1) +
    scale_x_continuous(breaks = scales::pretty_breaks(3),
                       labels = scales::comma) + # break the scales for each column so they don't overlap
    xlab(NULL) + ylab("Average prediction") +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(size = 8, family = "Helvetica Neue"),
      strip.text.y = element_text(size = 8, family = "Helvetica Neue"),
      axis.text = element_text(size = 8, family = "Helvetica Neue"),
      axis.title = element_text(size = 8, family = "Helvetica Neue"),
      axis.text.x = element_text(angle = 45, vjust = 0.5)
    ) # rotate axis labels to fit
  
  
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
      explain_tidymodels(input_model, data = pred, y = resp)
    return(explainer)
  }
  
  require(raster)
  # load covariate stack
  covs <- readRDS("data/covariate_stack.rds")
  # subset to only those covariates we are interested in
  covs <- subset(covs, model_covs)
  
  # create explainers for each model type
  glm_explainers <-
    map2(model_df$workflows[1:10], folds$splits, explainer_generator)
  gam_explainers <-
    map2(model_df$workflows[11:20], folds$splits, explainer_generator)
  rf_explainers <-
    map2(model_df$workflows[21:30], folds$splits, explainer_generator)
  brt_explainers <-
    map2(model_df$workflows[31:40], folds$splits, explainer_generator)
  # generate predictions for covariate stack
  glm_preds_stack <-
    lapply(glm_explainers, function(x)
      terra::predict(covs, x)) |> raster::stack()
  gam_preds_stack <-
    lapply(gam_explainers, function(x)
      terra::predict(covs, x)) |> raster::stack()
  rf_preds_stack <-
    lapply(rf_explainers, function(x)
      terra::predict(covs, x)) |> raster::stack()
  brt_preds_stack <-
    lapply(brt_explainers, function(x)
      terra::predict(covs, x)) |> raster::stack()
  
  # need final workflows from each model saved so that weights can be calculated...
  glm_auc_weights <- model_df |> filter(model == "GLM") |> pull(AUC)
  gam_auc_weights <- model_df |> filter(model == "GAM") |> pull(AUC)
  rf_auc_weights <- model_df |> filter(model == "RF") |> pull(AUC)
  brt_auc_weights <- model_df |> filter(model == "BRT") |> pull(AUC)
  
  # Create AUC weighted average of each model for the ensemble
  glm_preds_raster <- weighted.mean(glm_preds_stack, w = glm_auc_weights)
  gam_preds_raster <- weighted.mean(gam_preds_stack, w = gam_auc_weights)
  rf_preds_raster <- weighted.mean(rf_preds_stack, w = rf_auc_weights)
  brt_preds_raster <- weighted.mean(brt_preds_stack, w = brt_auc_weights)
  
  # ensemble mean
  ensemble_stack <-
    stack(glm_preds_raster,
          gam_preds_raster,
          rf_preds_raster,
          brt_preds_raster)
  auc_w <-
    c(
      mean(glm_auc_weights),
      mean(gam_auc_weights),
      mean(rf_auc_weights),
      mean(brt_auc_weights)
    )
  ensemble_raster <- weighted.mean(ensemble_stack, w = auc_w)
  plot(ensemble_stack)
  
  # Maps need to look pretty
  # could try this https://timogrossenbacher.ch/2016/12/beautiful-thematic-maps-with-ggplot2-only/
  
  # pretty gradient and polar mask functions
  source("R/discrete_gradient.R")
  source("R/polar_mask.R")
  polar_buffer <- polar_mask(radius_size = 5750000)
  
  # define projection
  crs_polar <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  
  # load in shapefile for background mask, clip and project
  world_shp <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
  CP <-
    sf::st_bbox(c(
      xmin = -180,
      xmax = 180,
      ymin = -90,
      ymax = 0
    ), crs = 4326) %>% sf::st_as_sfc()

  world_shp <- world_shp %>% sf::st_crop(CP)
  world_shp <- world_shp %>% sf::st_transform(crs_polar)
  world_shp <- world_shp %>% st_union()
  
  # load front shapefile
  fronts <- parkfronts |> 
    st_as_sf() |> 
    filter(front %in% c("SAF", "PF", "SACCF")) 
  
  # convert ensemble stack to a data frame for use with geom_raster
  ensemble_stack_df <-
    ensemble_stack |>
    rasterToPoints() |>
    as_tibble() |>
    rename(
      GLM = layer.1,
      GAM = layer.2,
      RF = layer.3,
      BRT = layer.4
    ) |>
    pivot_longer(3:6, names_to = "model", values_to = "preds")
  
  ensemble_stack_df <- ensemble_stack_df |>
    mutate(model = factor(model, levels = c("GLM", "GAM", "RF", "BRT")))
  
  # pretty facet plot for each model in the stack
  enst <- ggplot() +
    theme_void(base_size = 10,
               base_family = "Helvetica Neue") +
    geom_raster(aes(x = x, y = y, fill = preds),
                data = ensemble_stack_df) +
    geom_sf(
      aes(),
      data = world_shp,
      colour = "grey60",
      fill = "grey60"
    ) +
    geom_sf(
      data = polar_buffer$mask,
      fill = "white",
      color = "grey40",
      size = 0.5 / .pt
    ) +
    coord_sf(
      xlim = c(-6400000, 6400000),
      ylim = c(-6400000, 6400000),
      expand = FALSE,
      crs = crs_polar,
      ndiscr = 1000
    ) +
    xlab(NULL) + ylab(NULL) +
    facet_wrap( ~ model) +
    scale_fill_discrete_gradient(
      "Habitat Suitability",
      colours = viridis::viridis(10),
      bins = 10,
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = seq(0, 1, 0.2),
      guide = guide_colourbar(
        nbin = 500,
        raster = T,
        frame.colour = "grey40",
        ticks.colour = "grey40",
        frame.linewidth = .1,
        barwidth = 20,
        barheight = .5,
        direction = "horizontal",
        title.position = "top",
        #or "right"
        title.theme = element_text(hjust = 0.5,
                                   size = 10)
      )
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 10,
        family = "Helvetica Neue",
        color = "white"
      ),
      legend.text = element_text(size = 10, family = "Helvetica Neue"),
      legend.title = element_text(size = 10, family = "Helvetica Neue"),
      panel.background = element_rect(fill = "black"),
      legend.position = "none"
    )
  
  # pretty plot of the ensemble
  en <- ggplot() +
    theme_void(base_size = 10,
               base_family = "Helvetica Neue") +
    geom_raster(aes(x = x, y = y, fill = layer),
                data = rasterToPoints(ensemble_raster) |> as_tibble()) +
    geom_sf(
      aes(),
      data = world_shp,
      colour = "grey60",
      fill = "grey60"
    ) +
    geom_sf(
      aes(),
      data = fronts, 
      colour = "light grey", 
      linewidth = .1) +
    geom_text(aes(x = 2700000, y = -2100000, label = "SACCF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
    geom_text(aes(x = 3050000, y = -2550000, label = "PF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
    geom_text(aes(x = 3450000, y = -2950000, label = "SAF"), size = 2.5, family = "Helvetica Neue", angle = 40, colour = "grey") +
    geom_point(
      aes(x = x, y = y),
      data = points,
      colour = "black",
      fill = "white",
      shape = 21,
      size = .75) +
    geom_sf(
      data = polar_buffer$mask,
      fill = "white",
      color = "grey40",
      size = 0.5 / .pt
    ) +
    coord_sf(
      xlim = c(-6400000, 6400000),
      ylim = c(-6400000, 6400000),
      expand = FALSE,
      crs = crs_polar,
      ndiscr = 1000
    ) +
    xlab(NULL) + ylab(NULL) +
    #  ggtitle("Ensemble mean weighted by AUC") +
    scale_fill_discrete_gradient(
      "Habitat Suitability",
      colours = viridis::viridis(10),
      bins = 10,
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = seq(0, 1, 0.2),
      guide = guide_colourbar(
        nbin = 500,
        raster = T,
        frame.colour = "grey40",
        ticks.colour = "grey40",
        frame.linewidth = .1,
        barwidth = 20,
        barheight = .5,
        direction = "horizontal",
        title.position = "top",
        #or "right"
        title.theme = element_text(
          hjust = 0.5,
          size = 10,
          colour = "black"
        )
      )
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 10,
        family = "Helvetica Neue"
      ),
      legend.text = element_text(size = 10, family = "Helvetica Neue"),
      legend.title = element_text(size = 10, family = "Helvetica Neue"),
      panel.background = element_rect(fill = "black"),
      legend.position = "bottom"
    )
  
  # best way to output these for suppl material?
  require(patchwork)
  
  quartz(width = 7, height = 8.5)
  
  
  print(
    
    en / (enst + p_grid + plot_layout(widths = c(2, 3))) +
      plot_layout(heights = c(3, 2)) +
      plot_annotation(tag_levels = 'a') &
      theme(
        plot.tag = element_text(
          size = 12,
          family = "Helvetica Neue",
          face = "bold"
        ),
        plot.tag.position = c(.05, .95)
      )
  )
  
  quartz.save(
    file = plot_nm,
    type = "jpeg",
    dev = dev.cur(),
    dpi = 500
  )
  dev.off()
  
  saveRDS(ensemble_raster, file_nm) 
  
}
