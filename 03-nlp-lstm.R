library(keras)
library(abind)
library(tidyverse)
library(rsample)
library(data.table)
library(here)
library(colorspace)
library(ggplot2)
library(patchwork)
library(furrr)
library(yardstick)
source(here("R/lstm.R"))

sequence_length <- 25

config <- config::get()
random_seed <- config$seed

# read prepared lithologs ----
# splits <- board |>
#   pin_read("nlp-data-splits")
#
# unlabelled <- board |>
#   pin_read("nlp-newdata")
#
# picks <- board |>
#   pin_read("picks-bedrock-top")

picks_ww <- picks |>
  filter(
    source_type == "Water well",
    project != "Moose Lake"
  ) |>
  mutate(id = as.integer(id)) |>
  drop_na(id) |>
  rename(gicwellid = "id")

training_data <- training(splits)
testing_data <- testing(splits)

train <- training_data |>
  unite(col = "material", material, colour, na.rm = TRUE, sep = " ") |>
  select(gicwellid, int_top_dep, longitude, latitude, material, unit) |>
  mutate(unit = if_else(unit == "bedrock", 1, 0))

test <- testing_data |>
  unite(col = "material", material, colour, na.rm = TRUE, sep = " ") |>
  select(gicwellid, int_top_dep, longitude, latitude, material, unit) |>
  mutate(unit = if_else(unit == "bedrock", 1, 0))

new_data <- unlabelled |>
  unite(col = "material", material, colour, na.rm = TRUE, sep = " ") |>
  select(gicwellid, int_top_dep, longitude, latitude, material)

# encode material descriptions ----
text_vectorization <-
  layer_text_vectorization(output_mode = "multi_hot", max_tokens = 100)

# fit the preprocessing layer
adapt(text_vectorization, train$material)

# vector encode the train data
train_vec <- bind_cols(
  train,
  text_vectorization(train$material) |>
    as.matrix() |>
    as_tibble(.name_repair = "minimal") |>
    set_names(get_vocabulary(text_vectorization))
) |>
  select(-material)

test_vec <- bind_cols(
  test,
  text_vectorization(test$material) |>
    as.matrix() |>
    as_tibble(.name_repair = "minimal") |>
    set_names(get_vocabulary(text_vectorization))
) |>
  select(-material)

new_data_vec <- bind_cols(
  new_data,
  text_vectorization(new_data$material) |>
    as.matrix() |>
    as_tibble(.name_repair = "minimal") |>
    set_names(get_vocabulary(text_vectorization))
) |>
  select(-material)

# normalize numeric features
numeric_features <- c("int_top_dep", "longitude", "latitude")
numeric_cols <- which(names(train_vec) %in% numeric_features)

transform <- standardscaler(train_vec[, numeric_cols])

train_vec[, numeric_cols] <- predict(transform, train_vec[, numeric_cols])
test_vec[, numeric_cols] <- predict(transform, test_vec[, numeric_cols])
new_data_vec[, numeric_cols] <- predict(transform, new_data_vec[, numeric_cols])

# split data into each sequence ----
train_groups <- group_by(train_vec, gicwellid) |>
  group_split() |>
  set_names(unique(train$gicwellid))

test_groups <- group_by(test_vec, gicwellid) |>
  group_split() |>
  set_names(unique(test$gicwellid))

new_data_groups <- group_by(new_data_vec, gicwellid) |>
  group_split() |>
  set_names(unique(new_data$gicwellid))

# pad each sequence ----
plan(multisession)
train_padded <- future_map(train_groups, padding, sequence_length = sequence_length)
test_padded <- future_map(train_groups, padding, sequence_length = sequence_length)
new_data_padded <- future_map(new_data_groups, padding, sequence_length = sequence_length)
plan(sequential)

# split X, y ----
x_train <- map(train_padded, function(x) x[, colnames(x) != "unit"])
y_train <- map(train_padded, function(x) as.integer(x[, colnames(x) == "unit"]))

x_test <- map(test_padded, function(x) x[, colnames(x) != "unit"])
y_test <- map(test_padded, function(x) as.integer(x[, colnames(x) == "unit"]))

new_data <- map(new_data_padded, function(x) x[, colnames(x) != "unit"])

# data reshaping ----
# the data must be in the shape of (batch_size, seq_len, vocab)
n_features <- ncol(x_train[[1]])

x_train <- abind(x_train, along = 0)
y_train <- abind(y_train, along = 0)
input_shape <- dim(x_train)[2:3]

dim(y_train) <- c(dim(y_train), 1)

samples <- dim(x_train)[1]
timesteps <- dim(x_train)[2]
features <- dim(x_train)[3]

print(c(samples, timesteps, features))

# model ----
model <-
  keras_model_sequential() |>
  layer_masking(mask_value = -99999, input_shape = list(timesteps, features)) |>
  bidirectional(layer_lstm(units = 128, return_sequences = TRUE)) |>
  bidirectional(layer_lstm(units = 64, return_sequences = TRUE)) |>
  bidirectional(layer_lstm(units = 32, return_sequences = TRUE)) |>
  layer_dense(units = 1, activation = "sigmoid")

model |> compile(
  optimizer = "rmsprop",
  loss = "binary_crossentropy",
  metrics = c("accuracy")
)

summary(model)

history <- model |> fit(
  x = x_train,
  y = y_train,
  batch_size = 32,
  epochs = 5,
  validation_split = 0.2,
)

# predict test set ----
x_test <- abind(x_test, along = 0)
y_test <- abind(y_test, along = 0)

# convert to classification
y_hat <- predict(model, x_test)
y_hat <- y_hat > 0.5

# reshape predictions
dim(y_hat) <- dim(y_hat)[1:2]

# unpad (slice top n intervals)
logs_test <- testing(splits)
logs_test <- inner_join(logs_test, select(picks_ww, c(gicwellid, bedrock_dep)))

logs_test <- logs_test |>
  unite(col = "material", material, colour, na.rm = TRUE, sep = " ") |>
  select(gicwellid, int_top_dep, longitude, latitude, material, unit, bedrock_dep) |>
  mutate(unit = if_else(unit == "bedrock", 1, 0)) |>
  group_by(gicwellid) |>
  slice_head(n = sequence_length) |>
  ungroup()

logs_test <- logs_test |>
  mutate(
    unit = if_else(int_top_dep >= bedrock_dep, "bedrock", "surficial"),
    unit = factor(unit, levels = c("bedrock", "surficial"))
  )

y_hat <- unpad_sequences(y_hat, logs_test)
logs_test <- bind_cols(logs_test, y_hat)

logs_test <- logs_test |>
  mutate(.pred_class = str_to_sentence(.pred_class))

test_tops <- logs_test |> pick_bedrock()
logs_test <- inner_join(logs_test, test_tops)

logs_test |>
  ggplot(aes(bedrock_dep - .bedrock_dep)) +
  geom_histogram()

# predict unlabelled logs ----
# reshape
new_data <- abind(new_data, along = 0)

# predict
preds <- predict(model, new_data)
preds <- preds >= 0.5
dim(preds) <- dim(preds)[1:2]

preds_unpadded <- preds |>
  unpad_sequences(unlabelled)

predicted_logs <- unlabelled |>
  bind_cols(preds_unpadded)

preds_top <- predicted_logs |> pick_bedrock()
predicted_logs <- left_join(predicted_logs, preds_top)

View(predicted_logs |> select(id, int_top_dep, material, .pred_class, .bedrock_dep))
