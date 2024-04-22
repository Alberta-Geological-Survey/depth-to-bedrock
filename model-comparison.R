library(tidyverse)
library(tidymodels)

res100 <- readRDS("models/dtb-testset-edmonton100.rds") |>
  mutate(resolution = 100)
res250 <- readRDS("models/dtb-testset-edmonton250.rds") |>
  mutate(resolution = 250)
res500 <- readRDS("models/dtb-testset-edmonton500.rds") |>
  mutate(resolution = 500)

res <- bind_rows(res100, res250, res500)

res |>
  group_by(resolution) |>
  rmse(bedrock_dep, .pred) |>
  ggplot(aes(resolution, .estimate)) +
  geom_point() +
  geom_line() +
  xlab("Horizontal resolution [m]") +
  ylab("RMSE [m]")

res |>
  mutate(resolution = paste(resolution, "m")) |>
  ggplot(aes(bedrock_dep, .pred)) +
  geom_hex(binwidth = 1) +
  scale_fill_viridis_c(option = "inferno") +
  facet_wrap(vars(resolution)) +
  xlab("Observed depth [m]") +
  ylab("Predicted depth [m]")

