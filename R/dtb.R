#' Summarize the number of picks with specific data types, including 'Shotholes'
#' as boreholes and 'Contour' as 'Pseudo-observation'
#'
#' @param picks tibble of picks data
#'
#' @return tibble of summarized data with the number of observations per data
#'   type.
#' @export
summarize_data_type <- function(picks) {
  cnts <- picks |>
    mutate(data_type = fct_recode(
      data_type,
      "Pseudo-observation" = "Contour",
    )) |>
    group_by(data_type) |>
    tally()

  return(cnts)
}

#' Summarize the proportion of borehole sources in the picks data
#'
#' @param picks tibble of picks data.
#'
#' @return tibble of summarized data with the number and proportion of borehole
#'   data sources.
#' @export
summarize_borehole_sources <- function(picks) {
  source_perc <- picks |>
    filter(data_type == "Borehole") |>
    group_by(source_type) |>
    tally() |>
    mutate(pct = round(n / sum(n), 3) * 100)
  return(source_perc)
}

#' Basic preprocessing recipe for bedrock topography
#'
#' A tidymodels `recipes::recipe` object with the basic preprocessing operations
#' used for DTB prediction. This includes decomposing aspect (as radians) into
#' the strength of the northerly and easterly directions, and imputing missing
#' values.
#'
#' @param data tibble of training data
#'
#' @return a `recipe` specification
#' @export
recipe_dtb <- function(data) {
  data |>
    select(-id) |>
    recipes::recipe(bedrock_dep ~ .) |>
    recipes::step_mutate(
      aspect = dplyr::if_else(is.na(aspect), 0, aspect),
      easternness = cos(aspect * (pi / 180)),
      northerness = sin(aspect * (pi / 180))
    ) |>
    recipes::step_rm(dplyr::all_of("aspect")) |>
    recipes::step_impute_mean(recipes::all_numeric_predictors())
}

#' Preprocessing recipe that uses only terrain-related predictors
#'
#' @param data tibble of training data
#'
#' @return a `recipe` specification
#' @export
recipe_terrain <- function(data) {
  data |>
    recipe_dtb() |>
    recipes::step_rm(dplyr::any_of(c("xcoords", "ycoords")))
}

#' Preprocessing recipe that uses Euclidean distance fields as a spatial
#' feature engineering technique
#'
#' @param data tibble of training data
#'
#' @return a `recipe` specification
#' @export
recipe_edfs <- function(data) {
  tl <- tibble(xcoords = min(data$xcoords), ycoords = max(data$ycoords))
  tr <- tibble(xcoords = max(data$xcoords), ycoords = max(data$ycoords))
  bl <- tibble(xcoords = min(data$xcoords), ycoords = min(data$ycoords))
  br <- tibble(xcoords = max(data$xcoords), ycoords = min(data$ycoords))
  centre <- tibble(xcoords = mean(data$xcoords), ycoords = mean(data$ycoords))

  data |>
    recipe_dtb() |>
    spatialrecipes::step_spatial_dist2d(
      lat = "ycoords", lon = "xcoords", ref_lat = tl$ycoords, ref_lon = tl$xcoords,
      name = "edf1"
    ) |>
    spatialrecipes::step_spatial_dist2d(
      lat = "ycoords", lon = "xcoords", ref_lat = tr$ycoords, ref_lon = tr$xcoords,
      name = "edf2"
    ) |>
    spatialrecipes::step_spatial_dist2d(
      lat = "ycoords", lon = "xcoords", ref_lat = bl$ycoords, ref_lon = bl$xcoords,
      name = "edf3"
    ) |>
    spatialrecipes::step_spatial_dist2d(
      lat = "ycoords", lon = "xcoords", ref_lat = br$ycoords, ref_lon = br$xcoords,
      name = "edf4"
    ) |>
    spatialrecipes::step_spatial_dist2d(
      lat = "ycoords", lon = "xcoords", ref_lat = centre$ycoords, ref_lon = centre$xcoords,
      name = "edf5"
    ) |>
    recipes::step_rm(dplyr::all_of(c("xcoords", "ycoords")))
}

#' Preprocessing recipe that uses the values and distances to neighbouring
#' spatial locations as features
#'
#' @param data tibble of training data
#'
#' @return a `recipe` specification
#' @export
recipe_neighbours <- function(data) {
  recipe_dtb(data) |>
    recipes::step_normalize(recipes::all_numeric_predictors()) |>
    spatialrecipes::step_spatial_neighbors(
      dplyr::all_of(c("xcoords", "ycoords")),
      outcome = "bedrock_dep",
      neighbors = 16
    ) |>
    recipes::step_rm(dplyr::all_of(c("xcoords", "ycoords")))
}

#' Preprocessing recipe that uses spatial lag variables derived from the
#' Gaussian weighted mean of nearest neighbouring points as features
#'
#' @param data tibble of training data
#'
#' @return a `recipe` specification
#' @export
recipe_lags <- function(data) {
  recipe_dtb(data) |>
    recipes::step_normalize(recipes::all_numeric_predictors()) |>
    spatialrecipes::step_spatial_lag(
      dplyr::all_of(c("xcoords", "ycoords")),
      outcome = "bedrock_dep",
      weight_func = "gaussian",
      neighbors = 5
    ) |>
    spatialrecipes::step_spatial_lag(
      dplyr::all_of(c("xcoords", "ycoords")),
      outcome = "bedrock_dep",
      weight_func = "gaussian",
      neighbors = 10
    ) |>
    spatialrecipes::step_spatial_lag(
      dplyr::all_of(c("xcoords", "ycoords")),
      outcome = "bedrock_dep",
      weight_func = "gaussian",
      neighbors = 15
    ) |>
    recipes::step_rm(dplyr::all_of(c("xcoords", "ycoords")))
}

#' A `parsnip` model specification for random forests
#'
#' @return a `model_spec` object
#' @export
spec_rf <- function(r) {
  clf <- rand_forest(mtry = tune(), min_n = tune(), trees = 500L) |>
    set_mode("regression") |>
    set_engine(
      "ranger",
      seed = 1234,
      splitrule = tune(),
      num.threads = !!future::availableCores()
    )

  return(clf)
}

#' A `parsnip` model specification for XGBoost
#'
#' @return a `model_spec` object
#' @export
spec_xgb <- function(recipe) {
  clf <- boost_tree(
    trees = 500L,
    tree_depth = tune(),
    sample_size = tune(),
    learn_rate = tune(),
    mtry = tune(),
    stop_iter = 20L
  ) |>
    set_mode("regression") |>
    set_engine(
      "xgboost",
      tree_method = "gpu_hist",
      nthread = !!future::availableCores()
    )

  return(clf)
}

#' A `parsnip` model specification for Cubist regression trees
#'
#' @return a `model_spec` object
#' @export
spec_cubist <- function(recipe) {
  clf <- cubist_rules(
    committees = tune(),
    neighbors = tune()
  ) |>
    set_engine("Cubist")

  return(clf)
}

#' A `tune` hyperparameter tuning grid for use with the `rand_forest` model
#' specification
#'
#' @return a `tibble` with the tuning combinations.
#' @export
grid_rf <- function() {
  p <- parameters(
    mtry(c(5, 38)),
    min_n(c(1, 10)),
    splitrule = splitting_rule(c("extratrees", "variance"))
  )
  grid_regular(p, levels = c(3, 3, 2))
}

#' A `tune` hyperparameter tuning grid for use with the `cubist_rules` model
#' specification
#'
#' @return a `tibble` with the tuning combinations.
#' @export
grid_cubist <- function() {
  p <- parameters(
    committees(c(1, 50)),
    neighbors(c(1, 9))
  )
  grid_regular(p, levels = 3)
}

#' A `tune` hyperparameter tuning grid for use with the `boost_tree` model
#' specification
#'
#' @return a `tibble` with the tuning combinations.
#' @export
grid_xgb <- function() {
  p <- parameters(
    tree_depth(c(3, 21)),
    sample_prop(c(0.6, 1.0)),
    learn_rate(c(-10, -1)),
    mtry(c(5, 38))
  )
  set.seed(42)
  grid_latin_hypercube(p, size = 25)
}

#' Helper function to perform nested cross-validation using a recipe and
#' model combination.
#'
#' @param data tibble of training data with 'bedrock_dep' and predictor columns.
#' @param id character with the name of the modellign task.
#' @param recipe a `recipe` object used for preprocessing.
#' @param learner a list with two elements - a `model_spec` as the first
#'   element, and a tibble of tuning combinations as the second element.
#'
#' @return tibble of cross-validation results
#' @export
cross_validate_dtb <- function(data, id, recipe, learner) {
  data <- data |>
    select(-id)

  ctrl <- control_grid(
    save_pred = TRUE,
    allow_par = FALSE,
    save_workflow = TRUE,
    verbose = FALSE
  )

  model_spec <- learner[[1]]
  grid <- learner[[2]]

  wflow <- workflow() |>
    add_recipe(recipe) |>
    add_model(model_spec)

  # cross-validate
  p(message = sprintf("Cross-validating %s", id))

  set.seed(1)
  resamples <- nested_cv(
    data,
    outside = vfold_cv(v = 10L, strata = "bedrock_dep"),
    inside = vfold_cv(v = 3L, strata = "bedrock_dep")
  )

  cv_res <- cross_validate(wflow, resamples = resamples, grid = grid, control = ctrl)
  cv_res$id <- id
  return(cv_res)
}

#' Helper function to tune and fit a model specification
#'
#' @param data tibble of training data with 'bedrock_dep' and predictor columns.
#' @param id character with the name of the modellign task.
#' @param recipe a `recipe` object used for preprocessing.
#' @param learner a list with two elements - a `model_spec` as the first
#'   element, and a tibble of tuning combinations as the second element.
#'
#' @return a fitted model object
#' @export
train_dtb <- function(data, id, recipe, learner) {
  data <- data |>
    select(-id)

  ctrl <- control_grid(
    save_pred = TRUE,
    allow_par = FALSE,
    save_workflow = TRUE
  )

  model_spec <- learner[[1]]
  grid <- learner[[2]]

  wflow <- workflow() |>
    add_recipe(recipe) |>
    add_model(model_spec)

  # train model
  p(message = sprintf("Fitting %s", id))

  set.seed(44)
  resamples <- vfold_cv(data, v = 3L, strata = "bedrock_dep")
  model_metrics <- metric_set(rmse)

  tune_res <- tune_grid(
    wflow,
    resamples = resamples,
    grid = grid,
    metrics = model_metrics,
    control = ctrl
  )

  best_model <- finalize_workflow(wflow, select_best(tune_res, metric = "rmse"))
  final_model <- fit(best_model, data)

  return(final_model)
}

#' Helper function to perform DTB spatial prediction
#'
#' @param model a fitted `workflow` object.
#' @param predictors a wrapped `SpatRaster` object.
#' @param bbox optional `sf` bbox used to define the extent of the prediction.
#'
#' @return a wrapped `SpatRaster` object with the DTB predictions.
#' @export
predict_dtb <- function(id, model, predictors, bbox = NULL) {
  predictors <- terra::unwrap(predictors)

  if (!is.null(bbox)) {
    bbox <- as.numeric(bbox)
    terra::window(predictors) <- ext(bbox[1], bbox[3], bbox[2], bbox[4])
  }

  p(message = sprintf("Fitting %s", id))

  preds <- predict(predictors, model = model) |>
    setNames(id)

  preds <- terra::wrap(preds)
  terra::window(predictors) <- NULL
  return(preds)
}

#' Helper function to extract feature importances from a fitted workflow object
#'
#' Feature importances use permutation importances from the `vip` package.
#'
#' @param data tibble of training data with 'bedrock_dep' and predictor columns.
#' @param id character with the name of the modelling task.
#' @param workflow a fitted workflow object.
#'
#' @return tibble
#' @export
importances_dtb <- function(id, data, workflow) {
  prepped_recipe <- workflow |>
    extract_preprocessor() |>
    prep(data)

  p(message = sprintf("Cross-validating %s", id))

  fimp <- vip::vi_permute(
    object = extract_fit_parsnip(workflow),
    train = bake(prepped_recipe, new_data = NULL),
    target = "bedrock_dep",
    metric = "rsq",
    pred_wrapper = function(object, newdata)
      predict(object, newdata)[[1]]
  )

  fimp$id <- id
  return(fimp)
}
