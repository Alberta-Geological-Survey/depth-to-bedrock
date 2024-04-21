library(here)
library(glue)
library(tidyverse)
library(terra)
library(sf)
library(tidymodels)
library(rules)
library(spatialrecipes)
library(stacks)
library(nngeo)
source(here("R/dtb.R"))
source(here("R/resampling.R"))

tidymodels_prefer()
conflicted::conflicts_prefer(recipes::update)
conflicted::conflicts_prefer(terra::extract)
conflicted::conflicts_prefer(dplyr::filter)

dir.create(here("models"), showWarnings = FALSE)
dir.create(here("outputs"), showWarnings = FALSE)
ncores <- future::availableCores()

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
rec <- head(data) |>
  recipe_dtb() |>
  step_spatial_clusterdist(
    xcoords, ycoords, bedrock_dep,
    ref = c("xcoords", "ycoords"),
    num_comp = 10L
  )

rec_norm <- rec |> 
  step_normalize(all_numeric_predictors())

rec_crds <- rec |> 
  step_select(c(xcoords, ycoords, bedrock_dep))

rec_dem <- rec |> 
  step_select(c(xcoords, ycoords, dem, bedrock_dep)) |> 
  step_normalize(all_numeric_predictors())

rf <- rand_forest(trees = 500L, min_n = tune(), mtry = tune()) |>
  set_mode("regression") |>
  set_engine("ranger", seed = 1234, splitrule = "extratrees", num.threads = !!ncores)

xgb <- boost_tree(
    trees = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    sample_size = tune()
  ) |>
  set_mode("regression") |>
  set_engine("xgboost", nthread = !!ncores)

knn <- nearest_neighbor(weight_func = tune(), neighbors = tune()) |>
  set_mode("regression") |>
  set_engine("kknn")

## fit models ----
# tuning
set.seed(44)
resamples <- vfold_cv(data, v = 5L, strata = "bedrock_dep")

wflowset <- workflow_set(
  preproc = list(base = rec, base = rec, norm = rec_norm, crds = rec_crds,
                 crds_dem = rec_dem),
  models = list(rf, xgb, knn, knn, knn),
  cross = FALSE
)

wflowset <- wflowset |> 
  option_add(
    id = "base_rand_forest", 
    grid = grid_random(list(min_n = min_n(c(1, 10)), mtry = mtry(c(5, 41))))
  )

res <- workflow_map(
  wflowset,
  resamples = resamples,
  metrics = metric_set(rmse),
  control = control_stack_grid(),
  grid = 20L,
  seed = 44
)
autoplot(res)
write_rds(res, conf$tuning)

# stacking
model <- stacks() |>
  add_candidates(res) |>
  blend_predictions() |>
  fit_members()

write_rds(model, conf$model)

## test set predictions ----
ypreds <- last_fit(model, splits, metrics = metric_set(rmse))

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
