library(here)
library(tidyverse)
library(terra)
library(sf)
library(progressr)
library(tidymodels)
library(rules)
library(spatialrecipes)
source(here("R/functions/dtb.R"))
source(here("R/functions/resampling.R"))

# read data ----
data <- readRDS(here("data/processed/training-data.rds"))

bnd_ab <- st_read(here("projdata/bnd-ab.gpkg")) |>
  st_transform(3402)

predictors <- rast(here("data/processed/predictors.tif"))

data <- data |>
  st_as_sf(coords = c("xcoords", "ycoords"), crs = 3402, remove = FALSE) |>
  st_crop(bnd_ab) |>
  st_drop_geometry()

# machine learning ----
# define model
regr <- rand_forest(trees = 500L, min_n = tune(), mtry = tune()) |>
  set_mode("regression") |>
  set_engine(
    "ranger",
    seed = 1234,
    splitrule = "extratrees",
    num.threads = !!future::availableCores()
  )

grid <- regr |>
  extract_parameter_set_dials() |>
  update(
    mtry = mtry(c(5, 38)),
    min_n = min_n(c(1, 10))
  ) |>
  grid_regular(levels = 3)

rec <- recipe_lags(head(data))

wflow <- workflow() |>
  add_recipe(rec) |>
  add_model(regr)

## fit model ----
set.seed(44)
resamples <- vfold_cv(data, v = 3L, strata = "bedrock_dep")

ctrl <- control_grid()

tune_res <- tune_grid(
  wflow,
  resamples = resamples,
  grid = grid,
  metrics = metric_set(rmse),
  control = ctrl
)

autoplot(tune_res)

wflow_best <- finalize_workflow(wflow, select_best(tune_res, metric = "rmse"))
model <- fit(wflow_best, data)
write_rds(model, here("models/model-dtb-prov.rds"))

## nested cross-validation ----
set.seed(42)
resamples <- nested_cv(
  data,
  outside = vfold_cv(v = 10, strata = "bedrock_dep"),
  inside = vfold_cv(v = 3L, strata = "bedrock_dep")
)

cv_res <- cross_validate(
  wflow,
  resamples,
  grid = grid,
  metrics = metric_set(rmse, mae),
  metric = "rmse",
  control = ctrl
)

write_rds(cv_res, here("models/cross-validation-dtb-prov.rds"))

cv_res |>
  unnest(metrics) |>
  filter(.metric == "rmse") |>
  summarize(mean = mean(.estimate), std_err = sd(.estimate) / sqrt(n()))

# raster prediction ----
window(predictors) <- ext(vect(bnd_ab))

preds <-
  predict(predictors, model) |>
  setNames("DTB")

preds <- mask(preds, bnd_ab)
plot(preds, col = colorRampPalette(viridis::magma(9, direction = -1))(100))

writeRaster(preds, here("outputs/dtb-rf-prov.tif"), overwrite = TRUE)

# quantile prediction intervals ----
rec <- recipe_lags(data)

rec_prepped <- prep(
  rec,
  training = data,
  refresh = TRUE
)
prepped_data <- bake(rec_prepped, new_data = NULL)

rf_raw <-
  ranger::ranger(
    bedrock_dep ~ .,
    data = prepped_data,
    mtry = 5,
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
  recipe = rec_prepped
)

pred_int <- mask(pred_int, bnd_ab)
plot(pred_int, col = colorRampPalette(viridis::magma(9, direction = -1))(100))

writeRaster(pred_int, here("outputs/dtb-rf-prov-pred-int.tif"), overwrite = TRUE)
