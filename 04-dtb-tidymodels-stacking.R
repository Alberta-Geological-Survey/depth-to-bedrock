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
config <- "edmonton250"

# setup ----
conf <- config::get(config = config)

# read data ----
predictors <- rast(here(conf$predictors))

picks_comp <-
  read_csv(here("data/bedrock-depth-picks.csv")) |>
  mutate(pick_date = as_date(pick_date)) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(3402)

picks_nlp <-
  read_csv(here("outputs/bedrock-depth-picks-nlp.csv")) |>
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
  step_rm(-c(xcoords, ycoords, all_outcomes()))

rec_dem <- rec |>
  step_rm(-c(xcoords, ycoords, dem, all_outcomes())) |>
  step_normalize(all_numeric_predictors())

rf <- rand_forest(trees = 500L, min_n = tune(), mtry = tune()) |>
  set_mode("regression") |>
  set_engine("ranger", seed = 1234, splitrule = "extratrees",
             num.threads = !!ncores)

xgb <- boost_tree(
    trees = tune(),
    tree_depth = tune(),
    learn_rate = 0.1,
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
  grid = 10L,
  seed = 44
)
score_plot <- autoplot(res)
ggsave(here(conf$scoreplot), score_plot, width = 6, height = 4)
write_rds(res, conf$tuning)

# stacking
model <- stacks() |>
  add_candidates(res) |>
  blend_predictions() |>
  fit_members()

write_rds(model, here(conf$model))

## test set predictions ----
ypreds <- augment(model, testing(splits)) |>
  select(c(bedrock_dep, .pred, .resid))

rmse <- rmse(ypreds, bedrock_dep, .pred)
cat("rmse:", rmse$.estimate, "\n")

write_rds(ypreds, here(conf$testset))

# raster prediction ----
preds <- predict(predictors, model) |>
  setNames("DTB")
btopo <- predictors$dem - preds

writeRaster(preds, here(conf$dtb), overwrite = TRUE)
writeRaster(btopo, here(conf$btopo), overwrite = TRUE)
