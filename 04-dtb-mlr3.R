library(here)
library(data.table)
library(terra)
library(sf)
library(mlr3verse)
library(mlr3pipelines)
library(mlr3hyperband)
library(mlr3spatialops)
library(nngeo)
source(here("R/dtb.R"))

future::plan("multisession")
ncores = future::availableCores()

dir.create(here("models"), showWarnings = FALSE)
dir.create(here("outputs"), showWarnings = FALSE)

# setup ----
conf = config::get(config = "edmonton100")

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

# machine learning ----
# define task
task = TaskRegr$new(id = "dtb", backend = data, target = "bedrock_dep")
task$set_col_roles("id", roles = "name")

# define preprocessing
pop_function = po("mutate", mutation = list(
  aspect = ~ ifelse(is.na(aspect), 0, aspect),
  easternness = ~ cos(aspect * (pi / 180)),
  northerness = ~ sin(aspect * (pi / 180))
))

pop_impute = po("imputemean", affect_columns = selector_type("numeric"))

pop_cluster = PipeOpSpatialDist$new()
pop_cluster$param_set$values$k = 25L
pop_cluster$param_set$values$xcolname = "xcoords"
pop_cluster$param_set$values$ycolname = "ycoords"
pop_cluster$param_set$values$affect_columns = selector_name(c("xcoords", "ycoords", "bedrock_dep"))

preproc = pop_impute %>>%
  pop_function %>>%
  pop_cluster

# define model
rf = lrn("regr.ranger", num.threads = ncores)
xgb  = lrn("regr.xgboost", nthread = ncores)
knn = as_learner(po("scale") %>>% lrn("regr.kknn", kernel = "gaussian", k = 12))
knn_sp = as_learner(po("select", selector = selector_name(c("xcoords", "ycoords"))) %>>% knn)
regr = pipeline_stacking(
  base_learners = list(rf, xgb, knn, knn_sp),
  super_learner = lrn("regr.cv_glmnet"),
  folds = 3
)

pipeline = as_learner(preproc %>>% regr)

pipeline$param_set$set_values(
  regr.ranger.splitrule = to_tune(p_fct(c("variance", "extratrees"))),
  regr.ranger.mtry.ratio = to_tune(p_dbl(0.1, 1)),
  regr.ranger.min.node.size = to_tune(p_int(1, 10)),
  regr.xgboost.eta = to_tune(p_dbl(0.01, 0.3)),
  regr.xgboost.max_depth = to_tune(p_int(5, 20)),
  regr.xgboost.nrounds = to_tune(p_int(100, 500)),
  regr.xgboost.subsample = to_tune(p_dbl(0.5, 1)),
  geodist.k  = to_tune(p_int(5, 50))
)

instance = auto_tuner(
  tuner = tnr("random_search"),
  learner = pipeline,
  resampling = rsmp("cv", folds = 3L),
  measure = msr("regr.rmse"),
  terminator = trm("evals", n_evals = 20L)
)

# fit model ----
instance$train(task)
autoplot(instance$tuning_instance)
saveRDS(instance, here(conf$model))

# raster prediction ----
predfun = function(model, newdata, ...) {
  x = as.data.table(model$predict_newdata(newdata))
  return(x$response)
}

window(predictors) = aoi
preds = setNames(predict(predictors, instance, fun = predfun), "DTB")
btopo = predictors$dem - preds
window(predictors) = NULL

writeRaster(preds, here(conf$dtb), overwrite = TRUE)
writeRaster(btopo, here(conf$btopo), overwrite = TRUE)
