library(here)
library(glue)
library(tidyverse)
library(terra)
library(sf)
library(tidymodels)
library(rules)
library(spatialrecipes)
library(nngeo)
source(here("R/dtb.R"))
source(here("R/resampling.R"))

tidymodels_prefer()
conflicted::conflicts_prefer(recipes::update)
conflicted::conflicts_prefer(terra::extract)
conflicted::conflicts_prefer(dplyr::filter)

# flags ----
kwargs <-
  optparse::OptionParser() |>
  optparse::add_option("--config", default = "edmonton500", type = "character") |>
  optparse::parse_args()

config <- kwargs$config

# setup ----
conf <- config::get(config = config)

# read data ----
predictors <- rast(conf$predictors)

picks_comp <- read_csv("data/bedrock-depth-picks.csv") |>
  mutate(pick_date = as_date(pick_date)) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(3402)

picks_nlp <- read_csv("outputs/bedrock-depth-picks-nlp.csv") |>
  mutate(pick_date = as_date(pick_date)) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(3402)

picks_nlp <- picks_nlp |>
  mutate(
    pick_date = as.Date(pick_date),
    gr_elev_source = "DEM",
    source_table = "Autopicked"
  )

# extract training data ----
picks <-
  combine_picks(picks_comp, picks_nlp) |>
  st_crop(unlist(conf$region))

data <- extract_training(picks, predictors)

# machine learning ----
## split into training and test sets ----
set.seed(23)
splits <- initial_split(data, prop = 0.8, strata = "bedrock_dep")

## define model ----
regr <- rand_forest(trees = 300L, min_n = tune(), mtry = tune()) |>
  set_mode("regression") |>
  set_engine(
    "ranger",
    seed = 1234,
    splitrule = "extratrees",
    num.threads = !!future::availableCores()
  )

rec <- head(data) |>
  recipe_dtb() |>
  step_spatial_lag(
    all_of(c("xcoords", "ycoords")),
    outcome = "bedrock_dep",
    weight_func = "gaussian",
    neighbors = 5
  ) |>
  step_spatial_lag(
    all_of(c("xcoords", "ycoords")),
    outcome = "bedrock_dep",
    weight_func = "gaussian",
    neighbors = 10
  ) |>
  step_spatial_lag(
    all_of(c("xcoords", "ycoords")),
    outcome = "bedrock_dep",
    weight_func = "gaussian",
    neighbors = 15
  ) |>
  step_spatial_clusterdist(
    xcoords, ycoords, bedrock_dep,
    ref = c("xcoords", "ycoords"),
    num_comp = tune()
  )

wflow <- workflow() |>
  add_recipe(rec) |>
  add_model(regr)

grid <- wflow |>
  extract_parameter_set_dials() |>
  update(
    mtry = mtry(c(5, 38)),
    min_n = min_n(c(1, 10)),
    num_comp = num_comp(c(5, 50))
  ) |>
  grid_regular(levels = 3)

## fit model ----
# tuning
set.seed(44)
resamples <- vfold_cv(data, v = 10L, strata = "bedrock_dep")

tune_res <- tune_grid(
  wflow,
  resamples = resamples,
  grid = grid,
  metrics = metric_set(rmse),
  control = control_grid(verbose = TRUE, save_pred = TRUE, allow_par = FALSE),
)

dir.create(here("models"), showWarnings = FALSE)
write_rds(tune_res, conf$tuning)

# finalization
wflow_best <- finalize_workflow(wflow, select_best(tune_res, metric = "rmse"))
model <- fit(wflow_best, data)
write_rds(model, conf$model)

## test set predictions ----
ypreds <- last_fit(wflow_best, splits, metrics = metric_set(rmse))

rmse <- ypreds |>
  collect_metrics(summarize = TRUE) |>
  filter(.metric == "rmse") |>
  pull(.estimate)

cat("rmse:", rmse, "\n")
write_rds(ypreds, here(conf$testset)

# raster prediction ----
window(predictors) <- ext(unlist(conf$region))

preds <- predict(predictors, model) |>
  setNames("DTB")
writeRaster(preds, conf$dtb, overwrite = TRUE)

btopo <- predictors$dem - preds
writeRaster(btopo, conf$btopo, overwrite = TRUE)

# quantile prediction intervals ----
rec_prepped <- extract_recipe(model)
prepped_data <- bake(rec_prepped, new_data = data)

rf_raw <-
  ranger::ranger(
    bedrock_dep ~ .,
    data = prepped_data,
    mtry = 6,
    min.node.size = 1,
    keep.inbag = TRUE,
    num.threads = future::availableCores(),
    quantreg = TRUE
  )

pred_int <- predict(
  predictors,
  rf_raw,
  fun = function(model, newdata, recipe, type, quantiles) {
    newdata <- bake(recipe, new_data = newdata)
    preds <- predict(model, newdata, type = type, quantiles = quantiles)$predictions
    preds <- as.data.frame(preds)
    preds <- setNames(preds, paste0("q", quantiles))
    preds$.pred_int <- preds[[2]] - preds[[1]]
    preds
  },
  type = "quantiles",
  quantiles = c(0.1, 0.9),
  recipe = rec_prepped,
  filename = conf$predint,
  overwrite = TRUE
)
