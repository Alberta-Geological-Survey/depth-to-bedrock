library(here)
library(data.table)
library(terra)
library(sf)
library(mlr3verse)
library(mlr3pipelines)
library(mlr3hyperband)
library(mlr3spatialops)
library(ggplot2)
library(nngeo)
source(here("R/dtb.R"))

# setup ----
ncores = future::availableCores()
dir.create(here("models"), showWarnings = FALSE)
dir.create(here("outputs"), showWarnings = FALSE)

# setup ----
conf = config::get(config = "edmonton250")

# read data ----
predictors = rast(conf$predictors)

picks_comp = fread(here("data/bedrock-depth-picks.csv"))
picks_comp[, pick_date := as.Date(pick_date)]
picks_comp = picks_comp |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(3402)

picks_nlp = fread(here("outputs/bedrock-depth-picks-nlp.csv"))
picks_nlp[, pick_date := as.Date(pick_date)]
picks_nlp[, gr_elev_source := "DEM"]
picks_nlp[, source_table := "Autopicked"]

picks_nlp = picks_nlp |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(3402)

# extract training data ----
aoi = ext(unlist(conf$region))

picks = combine_picks(picks_comp, picks_nlp) |>
  vect() |>
  crop(aoi) |>
  as.data.frame(geom = "XY") |>
  st_as_sf(coords = c("x", "y"), crs = 3400)

data = as.data.table(extract_training(picks, predictors))
data = na.omit(data)

# machine learning ----
# define task
task = TaskRegr$new(id = "dtb", backend = data, target = "bedrock_dep")
task$add_strata("bedrock_dep")
task$set_col_roles("id", roles = "name")

# define preprocessing
pop_function = po("mutate", mutation = list(
  aspect = ~ ifelse(is.na(aspect), 0, aspect),
  easternness = ~ cos(aspect * (pi / 180)),
  northerness = ~ sin(aspect * (pi / 180))
))

pop_cluster = PipeOpSpatialDist$new(
  param_vals = list(
    k = 25,
    xcolname = "xcoords",
    ycolname = "ycoords",
    affect_columns = selector_name(c("xcoords", "ycoords", "bedrock_dep"))
  ))

# define model
rf = lrn("regr.ranger", num.threads = ncores, splitrule = "extratrees")
xgb = lrn("regr.xgboost", nthread = ncores, nrounds = 500, eta = 0.1,
          max_depth = 15, subsample = 0.67)

knn1 =
  po("select", selector = selector_name(c("xcoords", "ycoords", "dem"))) %>>%
  po("scale") %>>%
  lrn("regr.kknn", kernel = "gaussian", k = 12, id = "knn.dem")

knn2 = po("select", selector = selector_name(c("xcoords", "ycoords"))) %>>%
  lrn("regr.kknn", kernel = "gaussian", k = 12, id = "knn.crds12")

knn3 = po("select", selector = selector_name(c("xcoords", "ycoords"))) %>>%
  lrn("regr.kknn", kernel = "inv", k = 25, id = "knn.crds25")

library(mlr3extralearners)
cubist = lrn("regr.cubist", neighbors = 5, committees = 25)

base_learners = lapply(list(rf, xgb, knn1, knn2, knn3, cubist), as_learner)

stack = pipeline_stacking(
  base_learners = base_learners,
  super_learner = lrn("regr.cv_glmnet"),
  folds = 3,
  use_features = FALSE
)

stack = pop_function %>>%
  pop_cluster %>>%
  stack |>
  as_learner()

# fit model ----
stack$train(task)
saveRDS(stack, here(conf$model))

# performance estimation
set.seed(42)
train_ids = sample(task$row_ids, 0.8 * length(task$row_ids))
test_ids = setdiff(task$row_ids, train_ids)

stack$clone(deep = TRUE)$train(task, row_ids = train_ids)
yhat = stack$predict(task, row_ids = test_ids)

score_rmse = yhat$score(msr("regr.rmse"))
score_mae = yhat$score(msr("regr.mae"))

autoplot(yhat) +
  ggtitle(glue::glue("RMSE: round({score_rmse}, 2)"))

saveRDS(yhat, here(conf$testset))

# raster prediction ----
predfun = function(model, newdata, ...) {
  x = as.data.table(model$predict_newdata(newdata))
  return(x$response)
}

preds = setNames(predict(predictors, stack, fun = predfun, na.rm = TRUE), "DTB")
btopo = predictors$dem - preds

writeRaster(preds, here(conf$dtb), overwrite = TRUE)
writeRaster(btopo, here(conf$btopo), overwrite = TRUE)
