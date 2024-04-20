library(here)
library(tidyverse)
library(terra)
library(sf)
library(progressr)
library(tidymodels)
library(rules)
library(spatialrecipes)
source(here("R/dtb.R"))
source(here("R/resampling.R"))

# read data ----
predictors <- rast(here("data/processed/predictors.tif"))
data <- readRDS(here("data/processed/training-data.rds"))

# create experiment design ----
subdataset_rois <- list(
  mnts = st_bbox(c(xmin = 459473, ymin = 5621257, xmax = 530280, ymax = 5688716)),
  nwab = st_bbox(c(xmin = 220554, ymin = 6427403, xmax = 478118, ymax = 6656935)),
  saos = st_bbox(c(xmin = 633875, ymin = 6109081, xmax = 821422, ymax = 6295954)),
  wcab = st_bbox(c(xmin = 265324, ymin = 5923836, xmax = 463220, ymax = 6093134))
)

subdatasets <- map(
  subdataset_rois,
  ~ data |>
    st_as_sf(
      coords = c("xcoords", "ycoords"),
      crs = 3402,
      remove = FALSE
    ) |>
    st_crop(.x) |>
    st_drop_geometry()
)

models <- list(
  RF = spec_rf(),
  Cubist = spec_cubist(),
  XGBoost = spec_xgb()
)

grids <- list(
  RF = grid_rf(),
  Cubist = grid_cubist(),
  XGBoost = grid_xgb()
)

preprocs <- list(
  terrain = recipe_terrain(head(data)),
  coordinates = recipe_dtb(head(data)),
  geodists = recipe_edfs(head(data)),
  neighbors = recipe_neighbours(head(data)),
  lags = recipe_lags(head(data))
)

# create a design matrix for the experiments
design <- expand_grid(
  data = map2(subdatasets, subdataset_rois, ~ list(.x, .y)),
  learner = map2(models, grids, ~ list(.x, .y)),
  preproc = preprocs
)

design <- design |>
  mutate(id = paste(names(learner), names(preproc), names(data), sep = "_")) |>
  relocate(id, .before = "data")

design <- design |>
  unnest_wider(data, names_sep = "_") |>
  rename("data" = "data_1", "bbox" = "data_2")

# machine learning ----
## fit models ----
model_res <- with_progress({
  p <- progressor(along = design$id)
  mapply(
    train_dtb,
    id = design$id,
    data = design$data,
    recipe = design$preproc,
    learner = design$learner,
    SIMPLIFY = FALSE,
    USE.NAMES = TRUE
  )
})

write_rds(model_res, here("models/experiments-models.rds"))

# cross-validation ----
cv_res <- with_progress({
  p <- progressor(along = design$id)
  mapply(
    cross_validate_dtb,
    id = design$id,
    data = design$data,
    recipe = design$preproc,
    learner = design$learner,
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
})

cv_res <- list_rbind(cv_res)
write_rds(cv_res, here("models/experiments-cross-validation.rds"))

# raster prediction ----
preds <- with_progress({
  p <- progressor(along = names(model_res))
  mapply(
    predict_dtb,
    id = names(model_res),
    model = model_res,
    bbox = design$bbox,
    MoreArgs = list(predictors = wrap(predictors)),
    SIMPLIFY = FALSE,
    USE.NAMES = TRUE
  )
})

write_rds(preds, here("models/experiments-dtb.rds"))

# feature importances ----
fimp_res <- with_progress({
  p <- progressor(along = names(model_res))
  mapply(
    importances_dtb,
    id = names(model_res),
    data = design$data,
    workflow = model_res,
    SIMPLIFY = FALSE,
    USE.NAMES = TRUE
  )
})
fimp_res <- list_rbind(fimp_res)

write_rds(fimp_res, here("models/experiments-importances.rds"))
