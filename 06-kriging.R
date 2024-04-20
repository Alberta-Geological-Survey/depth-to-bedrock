library(here)
library(glue)
library(sf)
library(terra)
library(tidyverse)
library(tidymodels)
library(gstat)
library(future.apply)
library(doFuture)

# setup ----
# parallel processing setup
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

# fit ordinary kriging on DTB ----
v <- variogram(bedrock_dep ~ 1, data, cutoff = 60000)
m <- fit.variogram(v, model = vgm(model = "Exp", nugget = 10))
plot(v, m)

# cross-validation ----
set.seed(42)
folds <- vfold_cv(data, v = 10, strata = "bedrock_dep")

folds$predictions <- future_lapply(folds$splits, function(splits) {
  train <- training(splits)
  test <- testing(splits)
  v_train <- variogram(bedrock_dep ~ 1, train, cutoff = 60000)
  m_train <- fit.variogram(v_train, model = vgm(model = "Exp", nugget = 10))
  krige(
    formula = bedrock_dep ~ 1,
    locations = train,
    newdata = test,
    model = m_train,
    nmax = 16
  )
}, future.packages = c("gstat", "sf"))

folds$metrics <-
  map2(folds$splits, folds$predictions, function(rsplit, preds) {
    test <- testing(rsplit)
    test <- st_drop_geometry(test)
    test$.pred <- preds$var1.pred
    rmse(test, bedrock_dep, .pred)
  })
write_rds(folds, here("models/cross-validation-kriging-prov.rds"))

folds |>
  unnest(metrics) |>
  summarize(mean = mean(.estimate), std_err = sd(.estimate) / sqrt(n()))

# raster prediction ----
interpolate_gstat <- function(model, x, ...) {
  v <- st_as_sf(x, coords = c("xcoords", "ycoords"), crs = 3402)
  p <- predict(model, v, ...)
  as.data.frame(p)[, 1:2]
}

model <- gstat(
  formula = bedrock_dep ~ 1,
  model = m,
  data = data,
  nmax = 16
)
window(predictors) <- ext(vect(bnd_ab))
result <- interpolate(predictors, model = model, fun = interpolate_gstat)
result[result < 0] <- 0

window(predictors) <- NULL

result <- mask(result, bnd_ab)
writeRaster(result, here("outputs/kriging-prov.tif"), overwrite = TRUE)
