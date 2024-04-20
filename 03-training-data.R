library(here)
library(tidyverse)
library(terra)
library(sf)
library(nngeo)
library(conflicted)
source(here("R/dtb.R"))

conflicts_prefer(terra::extract)
conflicts_prefer(dplyr::filter)

# read data ----
predictors <- rast(here("data/processed/predictors.tif"))

picks_comp <-
  readRDS(here("projdata/picks.rds")) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(3402)

picks_nlp <-
  readRDS(here("data/processed/picks-nlp.rds")) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(3402)

# combine pick sets ----
## only use predicted picks where others are sparse
k <- ifelse(nrow(picks_nlp) < 100, nrow(picks_nlp), 100)
idx_within_dist <- st_nn(picks_comp, picks_nlp, maxdist = 3000, k = k)
idx_within_dist <- unique(unlist(idx_within_dist))

picks_nlp_thinned <- picks_nlp |>
  filter(!row_number() %in% idx_within_dist)

picks <- bind_rows(picks_comp, picks_nlp_thinned)
st_write(picks, here("outputs/picked-combined.gpkg"), delete_dsn = TRUE)

# extract training data ----
data <- predictors |>
  extract(vect(picks), ID = FALSE) |>
  as_tibble()

data <- bind_cols(
  data,
  picks |> as.data.frame() |> select(id, bedrock_dep)
)

# average dtb in case of multiple picks per cell
data <- data |>
  group_by(xcoords, ycoords) |>
  summarise(
    across(where(is.numeric), mean),
    across(where(is.character), list),
    across(where(is.logical), all)
  ) |>
  ungroup()

glimpse(data)

write_rds(data, here("data/processed/training-data.rds"))
