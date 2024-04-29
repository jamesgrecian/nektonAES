####################################################################
### Function to fit four model ensemble to presence absence data ###
####################################################################

# define a function to fit the ensemble models to each species in turn
# this is the same as the previous version, just loaded remotely rather than inline
# all species are now in the folds object so can run without tweaking function

# Similar to previous attempts but this time with weights

# weights calculation following
# https://community.rstudio.com/t/error-when-generating-predictions-from-a-model-with-case-weights/167166

# updated 2024-02-29
# weights in the RF don't really work
# adapt to a down-sampled RF that subsamples the background point at each split
# new rules following Valavi et al. 2021  https://doi.org/10.1111/ecog.05615
# recommend a "down-sample RF" to address the class imbalance between presence and background points

fit_ensemble <- function(dat, folds, term_table) {
  
  # define projection
  prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" 
  
  # which covariates selected for species...
  print(term_table |> arrange(desc(mean)))
  
  tms <- term_table |> arrange(desc(mean)) |> dplyr::slice(1) |> pull(terms)
  
  # how many covariates selected
  n_covs <- length(strsplit(tms, "\\s\\+\\s", perl = TRUE)[[1]])
  
  ####################
  ### Ensembles... ###
  ####################
  
  # plan future
  #plan(multisession)
  
  ###############
  ### Fit GLM ###
  ###############
  
  # For GLMs optimise the combination of polynomial model terms to maximize model performance in terms of AUC
  # Fit up to and including third-order polynomials for the five predictor variables = 243 candidate models
  # Cross-validate across 10 spatial folds to calculate AUC
  # Select formula based on combination of polynomial terms that maximized AUC
  source("R/glm_poly.R") # generate all polynomial combinations of covariates
  
  # calculate possible polynomial combinations for the covariates
  poly_f <- glm_poly(cov_terms = strsplit(tms, "\\s\\+\\s", perl = TRUE)[[1]],
                     poly_degree = 3)
  
  writeLines(paste("number of polynomial combinations to try is:", length(poly_f)))
  
  # following on from discussion with Mike
  # use one recipe for fold CV
  # then another recipe for fitting the final model
  
  # define the formula
  glm_form <- as.formula(paste("PresAbs ~", tms))
  
  # define the first recipe
  # annoyingly need to pass the dataframe to the recipe even if though its not used...
  first_recipe <- recipes::recipe(
    dat,
    vars = c(all.vars(glm_form), "cwts"),  # Save the case weights column alongside our formula
    roles = c("outcome",
              rep("predictor", length(
                all.vars(glm_form)
              ) - 1),
              "case_weights")
  )
  
  cv_recipe <- first_recipe |>
    recipes::step_mutate(
      cwts = hardhat::importance_weights(
        ifelse(PresAbs == 1, 1, sum(PresAbs == 1) / sum(PresAbs == 0))
      ),
      role = "case_weight")  # Need to set the "case_weights" role explicitly:
  
  # define the glm model
  glm_model <- logistic_reg() |>
    set_engine("glm") |>
    set_mode("classification")
  
  # seems impossible to pass both the case weights and a list of formulas
  # instead loop through the formulas in turn
  # rank them and then re-run the best model
  results_summary <- list()
  
  # set up progress bar
  pb1 <- txtProgressBar(min = 1,
                        max = length(poly_f),
                        style = 3)
  
  for (i in 1:length(poly_f)) {
    # Set up a workflow - varying the formula with each iteration
    glm_wflow_wts <- workflow() |>
      add_model(glm_model,
                formula = poly_f[[i]]) |>
      add_case_weights(cwts) |>
      add_recipe(cv_recipe) |> # When doing anything relating to cross-validation use the sub-recipe
      fit_resamples(folds)
    
    setTxtProgressBar(pb1, i) # update progress bar
    
    results_summary[[i]] <-
      glm_wflow_wts |> collect_metrics() |> mutate(model = i)
    
  }
  
  # identify best fitting polynomial combination
  results_summary <- results_summary |> bind_rows()
  print(results_summary |> filter(.metric == "roc_auc") |> arrange(desc(mean)))
  
  mid <- results_summary |> 
    filter(.metric == "roc_auc") |>
    arrange(desc(mean)) |>
    dplyr::slice(1) |>
    pull(model)
  
  writeLines("Fitting final glm...")
  # define new recipe based on best model
  final_glm_wflow_wts <- workflow() |>
    add_model(glm_model,
              formula = poly_f[[mid]]) |>
    add_case_weights(cwts) |>
    add_recipe(first_recipe)
  
  # apply new recipe to each fold separately
  # calculate weights for each fold on the fly
  glm_final_model_fits_list <-
    lapply(folds$splits,
           FUN = function(x) fit(final_glm_wflow_wts, # the primary workflow
                                 analysis(x) |> # function to apply to each fold and manipulate case weights
                                   mutate(cwts = hardhat::importance_weights(
                                     ifelse(PresAbs == 1, 1, sum(PresAbs == 1) / sum(PresAbs == 0))
                                   ))
           )
    )
  
  ###############
  ### Fit GAM ###
  ###############
  
  # define GAM formula
  gam_form <-
    as.formula(paste(
      "PresAbs ~",
      paste0(
        "s(",
        strsplit(tms, "\\s\\+\\s", perl = TRUE)[[1]],
        ", bs = 'ts', k = 5)",
        collapse = " + "
      )
    ))
  
  # define gam model
  gam_spec <-
    gen_additive_mod() |>
    set_engine("mgcv", method = "REML") |>
    set_mode("classification")
  
  gam_wf_wts <-
    workflow() |>
    add_model(gam_spec,
              formula = gam_form) |>
    add_case_weights(cwts) |>
    add_recipe(first_recipe)
  
  writeLines("Fitting final gam...")
  
  # fit model to each fold
  gam_final_model_fits_list <-
    lapply(
      folds$splits,
      FUN = function(x)
        fit(
          gam_wf_wts,
          # the primary workflow
          analysis(x) |> # function to apply to each fold and manipulate case weights
            mutate(cwts = hardhat::importance_weights(ifelse(
              PresAbs == 1, 1, sum(PresAbs == 1) / sum(PresAbs == 0)
            )))
        )
    )
  
  ##############
  ### Fit RF ###
  ##############
  
  # set up an RF model using randomForest
  # mtry -> vary the number of parameters randomly tested and sampled at each branch between 1 and 3
  # fix the minimum number of data points allowed in a node (prevent branch splitting further) to 5
  # allow the number of trees to vary between 500 and 2500 in increments of 500?
  # rank all these models by AUC and select the best
  # both Reisinger and Tickly use randomForest package
  #plan(multisession)
  
  prNum <- as.numeric(table(dat$PresAbs)["1"])
  bgNum <- as.numeric(table(dat$PresAbs)["0"])
  
  rf_spec <- 
    rand_forest(
      mtry = tune(),       # vary number of predictors to randomly sample at each split
      min_n = 10,       # min number of data points in a node for node to be split further
      trees = tune(),      # vary number of trees to contain in the ensemble
      mode = "classification") |>
    set_engine("ranger",
               num.threads = 8,
               probability = TRUE, # probability tree
               sample.fraction = prNum/bgNum) # down-sample the background points
  
  # create workflow and specify model formula
  rf_wf <-
    workflow() |>
    add_model(rf_spec) |>
    add_case_weights(cwts) |>
    add_recipe(cv_recipe) # When doing anything relating to cross-validation use the sub-recipe
  
  writeLines("Tuning Random Forest...")
  
  # tune the randomForest model manually
  # allow mtry and number of trees to vary...
  tune_rf <- rf_wf |>
    tune_grid(resamples = folds,
              grid = expand.grid(mtry = c(2:(n_covs - 1)),
                                 # allow no. of covariates to randomly sample at each split to vary between 2 and n-1
                                 trees = seq(1000, 5000, by = 500))) # allow no. of trees to vary between 1000 and 2500
  
  # rank models by AUC
  tune_rf |>
    collect_metrics() |>
    filter(.metric == "roc_auc") |>
    arrange(desc(mean))
  
  # plot to check
  tune_rf |> autoplot(type = "marginals")
  
  # select best model...
  best_model_rf <-
    tune_rf |> select_best(metric = "roc_auc")
  
  # finalise workflow ready to for final fit
  # this creates a new workflow using the first recipe
  # not the cv recipe
  final_model_rf <-
    rf_wf |>
    update_recipe(first_recipe) |> # update to use the first recipe
    finalize_workflow(best_model_rf)
  
  writeLines("Fitting Random Forest...")
  
  # fit model to each fold
  rf_final_model_fits_list <-
    lapply(
      folds$splits,
      FUN = function(x)
        fit(
          final_model_rf, # the primary workflow
          analysis(x) |> # function to apply to each fold and manipulate case weights
            mutate(cwts = hardhat::importance_weights(
              ifelse(PresAbs == 1, 1, sum(PresAbs == 1) / sum(PresAbs == 0)))
            )))
  
  ###############
  ### Fit BRT ###
  ###############
  
  # similar cross-validation approach was used to parameterise the BRT model
  # shrinkage parameter set to 0.001
  # number of trees varied between 1000 and 10000
  # tree complexity allowed to vary between 1 and 4
  # mtry vary between 2 and 4
  #plan(multisession)
  
  brt_spec <-
    boost_tree(
      mtry = tune(),
      # number of covariates to randomly sample at each split
      trees = tune(),
      # vary trees
      min_n = 10,
      # min number of data points in a node for it to be split further Tickley 10, Reisinger 20
      tree_depth = tune(),
      # vary branches
      learn_rate = .001,
      # learning rate is also known as the shrinkage parameter Tickley 0.001
      sample_size = .8
    ) |>         # proportion of data to sample at each step - smaller = more randomness in trees
    set_engine("xgboost",
               counts = T) |>
    set_mode("classification")
  
  brt_wf <-
    workflow() |>
    add_model(brt_spec) |>
    add_case_weights(cwts) |>
    add_recipe(cv_recipe) # When doing anything relating to cross-validation use the sub-recipe
  
  #dials_random <- grid_random(
  #  mtry(c(2, 4)),         # keep mtry the same as before
  #  trees(c(1000, 10000)), # tickly 5000 Ryan 10000
  #  tree_depth(c(1, 4)),   # number of branches allowed in tree
  #  size = 100) # only 10 models to try - should increase this but computationally v expensive...
  
  dials_grid <- expand.grid(
    mtry = c(2:(n_covs - 1)),
    trees = seq(1000, 10000, by = 1000),
    tree_depth = c(1:4)
  )
  
  writeLines("Tuning Boosted Regression Tree...")
  
  # manually tune the brt model
  tune_brt <- brt_wf |>
    tune_grid(resamples = folds,
              grid = dials_grid)
  
  # rank all these models by AUC
  tune_brt |>
    collect_metrics() |>
    filter(.metric == "roc_auc") |>
    arrange(desc(mean))
  
  tune_brt |> autoplot(type = "marginals")
  
  best_model_brt <-
    tune_brt |> select_best(metric = "roc_auc")
  
  # finalise workflow ready to for final fit
  # this creates a new workflow using the first recipe
  # not the cv recipe
  final_model_brt <-
    brt_wf |>
    update_recipe(first_recipe) |> # When doing anything relating to cross-validation use the sub-recipe
    finalize_workflow(best_model_brt)
  
  writeLines("Fitting Boosted Regression Tree...")
  
  # fit model to each fold
  brt_final_model_fits_list <-
    lapply(
      folds$splits,
      FUN = function(x)
        fit(
          final_model_brt,
          # the primary workflow
          analysis(x) |> # function to apply to each fold and manipulate case weights
            mutate(cwts = hardhat::importance_weights(ifelse(
              PresAbs == 1, 1, sum(PresAbs == 1) / sum(PresAbs == 0)
            ))))
    )
  
  ################################
  ### Save models and AUC info ###
  ################################
  
  writeLines("Writing output...")
  
  # calculate AUC for each fold
  glm_auc_weights <-
    final_glm_wflow_wts |>
    update_recipe(cv_recipe) |>
    fit_resamples(folds) |>
    unnest(.metrics) |>
    filter(.metric == "roc_auc") |>
    pull(.estimate)
  
  gam_auc_weights <- gam_wf_wts |>
    update_recipe(cv_recipe) |>
    fit_resamples(folds) |>
    unnest(.metrics) |>
    filter(.metric == "roc_auc") |>
    pull(.estimate)
  
  rf_auc_weights <- final_model_rf |>
    update_recipe(cv_recipe) |>
    fit_resamples(folds) |>
    unnest(.metrics) |>
    filter(.metric == "roc_auc") |>
    pull(.estimate)
  
  brt_auc_weights <-
    final_model_brt |>
    update_recipe(cv_recipe) |>
    fit_resamples(folds) |>
    unnest(.metrics) |>
    filter(.metric == "roc_auc") |>
    pull(.estimate)
  
  model_out <- bind_rows(
    tibble(
      species = unique(dat$species)[1],
      model = "GLM",
      fold = 1:10,
      AUC = glm_auc_weights,
      workflows = glm_final_model_fits_list
    ),
    tibble(
      species = unique(dat$species)[1],
      model = "GAM",
      fold = 1:10,
      AUC = gam_auc_weights,
      workflows = gam_final_model_fits_list
    ),
    tibble(
      species = unique(dat$species)[1],
      model = "RF",
      fold = 1:10,
      AUC = rf_auc_weights,
      workflows = rf_final_model_fits_list
    ),
    tibble(
      species = unique(dat$species)[1],
      model = "BRT",
      fold = 1:10,
      AUC = brt_auc_weights,
      workflows = brt_final_model_fits_list
    )
  )
  return(model_out)
}
