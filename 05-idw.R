library(here)
library(glue)
library(sf)
library(terra)
library(tidyverse)
library(tidymodels)
library(gstat)
library(future)
library(doFuture)

# setup ----
plan("multisession")

# read datasets ----
data <- readRDS(here("data/processed/training-data.rds"))

bnd_ab <- st_read(here("projdata/bnd-ab.gpkg")) |>
  st_transform(3402)

predictors <- rast(here("data/processed/predictors.tif"))

data <- data |>
  select(-id)

data <- data |>
  mutate(bedrock_elev = dem - bedrock_dep) |>
  st_as_sf(coords = c("xcoords", "ycoords"), crs = 3402) |>
  st_crop(bnd_ab)

# cross-validation ----
set.seed(42)
folds <- vfold_cv(data = data, v = 10, strata = "bedrock_elev")

folds$predictions <- map(folds$splits, function(splits) {
  train <- training(splits)
  test <- testing(splits)

  idw(
    formula = bedrock_elev ~ 1,
    locations = as(train, "Spatial"),
    newdata = as(test, "Spatial"),
    nmax = 16,
    idp = 2.0
  ) |>
    as.data.frame() |>
    as_tibble() |>
    select(var1.pred)
})

folds$metrics <- map2(
  folds$splits,
  folds$predictions,
  function(rsplit, preds) {
    test <- testing(rsplit)
    test$.pred <- preds$var1.pred
    metric_set(rmse, mae, rsq)(test, bedrock_elev, .pred)
  }
)

saveRDS(folds, here("models/cross-validation-idw-prov.rds"))

# raster interpolation ----
interpolate_gstat <- function(model, x, crs, ...) {
  v <- st_as_sf(x, coords = c("x", "y"), crs = crs)
  p <- predict(model, v, ...)
  as.data.frame(p)[, 1:2]
}

model <- gstat(formula = bedrock_elev ~ 1, data = data, nmax = 16)

window(predictors) <- ext(vect(bnd_ab))

result <- interpolate(
  predictors[[c("xcoords", "ycoords")]],
  model = model,
  fun = interpolate_gstat,
  crs = crs(predictors)
)
result <- result[[1]]

window(predictors) <- NULL

result <- mask(result, bnd_ab)
writeRaster(result, here("outputs/idw-prov.tif"), overwrite = TRUE)
