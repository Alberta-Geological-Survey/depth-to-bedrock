#' Simplify variable names for plotting
#'
#' @param df tibble containing a 'Variable' columns
#'
#' @return recoded data
#' @export
recode_variables <- function(df) {
  df |>
    mutate(
      Variable = str_remove(Variable, "^euclidean_"),
      Variable = str_remove(Variable, "sur_refl_"),
      Variable = str_replace(
        Variable,
        "bedrock_dep_lag_[0-9]+_gaussian",
        paste0("lag", str_extract(Variable, "[0-9]+"))
      ),
      Variable = str_replace(
        Variable,
        "nn_bedrock_dep_dist_bedrock_dep_[0-9]+",
        paste0("nndist", str_extract(Variable, "[0-9]+"))
      ),
      Variable = str_replace(
        Variable,
        "nn_bedrock_dep_[0-9]+",
        paste0("nn", str_extract(Variable, "[0-9]+"))
      )
    )
}

#' Plot raster grids of terrain-related variables
#'
#' @param object RasterStack object.
#' @return ggplot object.
#' @export
plot_terrain_predictors <- function(object) {
  subset_grids <- list(
    object$slope,
    object$c_long,
    object$c_cros,
    object$tri,
    object$vrm,
    object$mrvbf,
    object$mrrtf,
    object$openness_pos,
    object$ho,
    object$hu,
    object$nh,
    object$texture,
    object$vdchn5,
    object$vdepth,
    object$swi,
    object$mtpi
  )

  subset_grids_names <- c(
    "Slope [radians]",
    "Curvature (longitudinal)",
    "Curvature (cross-sectional)",
    "Terrain ruggedness index [m]",
    "Vector ruggedness measure",
    "Multires idx valley bot. flatness",
    "Multires idx ridge top flatness",
    "Topographic openness",
    "Height over drainage [m]",
    "Height under summits [m]",
    "Normalized height",
    "Texture",
    "Height above channel network [m]",
    "Valley depth [m]",
    "SAGA wetness index",
    "Multiscale topographic position idx"
  )

  subset_grids_pal <- list(
    slope = scale_fill_gradientn(
      colours = cpt("grass_slope"),
      na.value = "transparent"
    ),
    c_long = scale_fill_gradientn(
      colours = cpt("grass_curvature"),
      na.value = "transparent"
    ),
    c_cros = scale_fill_distiller(
      palette = "Spectral",
      na.value = "transparent"
    ),
    tri1 = scale_fill_continuous_sequential(
      palette = "Viridis",
      na.value = "transparent"
    ),
    vrm1 = scale_fill_continuous_sequential(
      palette = "Sunset",
      na.value = "transparent",
      trans = "sqrt"
    ),
    mrvbf = scale_fill_viridis_c(
      option = "inferno",
      na.value = "transparent",
      direction = -1
    ),
    mrrtf = scale_fill_viridis_c(
      option = "magma",
      na.value = "transparent",
      direction = -1
    ),
    openness_pos = scale_fill_continuous_sequential(
      palette = "PuBuGn",
      na.value = "transparent",
      rev = FALSE,
      trans = "sqrt"
    ),
    ho = scale_fill_continuous_sequential(
      palette = "Reds 3",
      na.value = "transparent"
    ),
    hu = scale_fill_continuous_sequential(
      palette = "Blues 3",
      na.value = "transparent"
    ),
    nh = scale_fill_continuous_diverging(
      palette = "Blue-Red",
      mid = 0.5,
      na.value = "transparent"
    ),
    texture = scale_fill_viridis_c(
      na.value = "transparent",
      option = "plasma"
    ),
    vdchn5 = scale_fill_gradientn(
      colours = cpt("grass_bcyr"),
      na.value = "transparent"
    ),
    vdepth = scale_fill_continuous_sequential(
      palette = "YlGnBu",
      na.value = "transparent"
    ),
    swi = scale_fill_continuous_diverging(
      palette = "Blue-Red 2",
      mid = 13,
      na.value = "transparent",
      rev = TRUE
    ),
    mtpi = scale_fill_gradientn(
      colours = cpt("grass_haxby"),
      na.value = "transparent"
    )
  )

  raster_plot <- function(grd, title, pal) {
    ggplot() +
      geom_spatraster(data = grd) +
      pal +
      ggtitle(title) +
      theme(
        axis.title = element_blank(),
        legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(ifelse(is.null(levels(layer)), 0.5, 0.2), "cm"),
        legend.box.spacing = unit(0, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        axis.text = element_blank(),
        plot.margin = ggplot2::margin(t = 0, b = 6, r = 0, l = 12),
        plot.title = element_text(size = 10)
      )
  }

  p <- pmap(list(subset_grids, subset_grids_names, subset_grids_pal), raster_plot)
  wrap_plots(p)
}

#' Plots cross sections
#'
#' @param profile_line The `sf` object containing a line geometry to use for the
#'   profile.
#' @param grids A `RasterStack` or `RasterBrick` containing the bedrock
#'   topography grids.
#' @param dem A `RasterLayer` object containing the dem to use to show the
#'   ground surface elevation on the cross-section.
#'
#' @return ggplot.
#' @export
plot_profiles <- function(profile_line, grids, dem, label) {
  # create points along line
  profile_pts <- profile_line |>
    st_line_sample(density = 1 / 500, type = "regular") |>
    st_cast("POINT") |>
    st_sf(crs = 3402) |>
    st_set_geometry("geometry")

  profile_pts$dist <- st_distance(profile_pts[1, ], profile_pts)[1, ]
  profile_pts$dist <- as.numeric(profile_pts$dist)

  # create points at each end of the line for labels
  line_ends <- profile_line |>
    st_cast("POINT") |>
    slice(1, n())

  # make a square bbox around the line
  profile_bbox <- profile_line |>
    st_buffer(10000) |>
    st_bbox()

  max_dimension <- which.max(c(
    abs(profile_bbox[1] - profile_bbox[3]),
    abs(profile_bbox[2] - profile_bbox[4])
  ))

  if (max_dimension == 1) {
    max_idx <- c(1, 3)
    min_idx <- c(2, 4)
  } else {
    max_idx <- c(2, 4)
    min_idx <- c(1, 3)
  }
  max_length <- abs(profile_bbox[max_idx[1]] - profile_bbox[max_idx[2]])
  min_length <- abs(profile_bbox[min_idx[1]] - profile_bbox[min_idx[2]])
  diff_length <- max_length - min_length

  profile_bbox[min_idx[1]] <- profile_bbox[min_idx[1]] - (diff_length / 2)
  profile_bbox[min_idx[2]] <- profile_bbox[min_idx[2]] + (diff_length / 2)
  profile_extent <- ext(as.numeric(profile_bbox)[c(1, 3, 2, 4)])

  # sample the grids
  dem_resampled <- terra::resample(dem, grids)
  extracted_df <- terra::extract(
    c(grids, dem_resampled),
    st_coordinates(profile_pts)
  )
  names(extracted_df)[ncol(extracted_df)] <- "Land surface"
  profile_pts <- bind_cols(profile_pts, extracted_df)

  # summarize min, max elevations to limit y axis
  min_elev <- profile_pts |>
    st_drop_geometry() |>
    select(-dist) |>
    min()

  max_elev <- profile_pts |>
    st_drop_geometry() |>
    select(-dist) |>
    max()

  # plot the profile
  n_colors <- nlyr(grids)
  pal <- qualitative_hcl(n = n_colors, palette = "Dark 2")
  pal[4] <- "black"

  p_profile <- profile_pts |>
    st_drop_geometry() |>
    pivot_longer(-c(dist, picks)) |>
    mutate(control = "Geological Pick") |>
    ggplot(aes(x = dist, y = value, colour = name)) +
    geom_line() +
    scale_colour_manual(values = pal) +
    geom_point(aes(y = picks, fill = control), shape = 21, stroke = 0, size = 2) +
    scale_fill_manual(values = "black") +
    ylim(c(min_elev, max_elev)) +
    xlab("Distance (m)") +
    ylab("EL (masl)") +
    ggtitle("Cross-section Profile") +
    labs(colour = "Model", shape = "Control") +
    labs(fill = "Control", colour = "Surface") +
    annotate("text", x = 0, y = max(profile_pts$dem), label = label) +
    annotate("text",
      x = max(profile_pts$dist), y = max(profile_pts$dem),
      label = label
    ) +
    theme(plot.title = element_text(size = 11))

  p_map <- grids |>
    select(-picks) |>
    crop(profile_extent) |>
    as.data.frame(xy = TRUE) |>
    pivot_longer(-c(x, y)) |>
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = value)) +
    facet_wrap(vars(name)) +
    tidyterra::scale_fill_whitebox_c(na.value = "transparent") +
    labs(fill = "EL (masl)") +
    geom_sf(data = profile_line, linetype = "dashed") +
    geom_sf(data = line_ends) +
    geom_sf_label(data = line_ends, aes(label = label)) +
    coord_sf(datum = sf::st_crs(3402)) +
    ggtitle("Bedrock Topography") +
    theme(
      axis.title = element_blank(),
      plot.title = element_text(size = 11),
      axis.text = element_text(size = 6),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 8)
    )

  profile_bounding_box <- profile_extent |>
    st_bbox() |>
    st_as_sfc() |>
    st_set_crs(3400)

  p_map / p_profile
}

#' Plot truth vs predicted
#'
#' @param preds tibble with columns 'bedrock_dep' and '.pred'
#'
#' @return ggplot
#' @export
plot_truth_predicted <- function(preds) {
  p1 <- ggplot(preds, aes(x = bedrock_dep, y = .pred)) +
    geom_bin_2d(binwidth = 5) +
    geom_abline(slope = 1, intecept = 0, linetype = 3) +
    geom_smooth(method = "lm", colour = "black") +
    scale_fill_viridis_c(option = "magma", direction = -1, trans = "sqrt") +
    labs(fill = "Count") +
    xlab("Bedrock Depth [m]") +
    ylab("Predicted [m]")

  p2 <- ggplot(preds, aes(x = bedrock_dep, y = .pred - bedrock_dep)) +
    geom_bin_2d(binwidth = 5, show.legend = FALSE) +
    scale_fill_viridis_c(option = "magma", direction = -1, trans = "sqrt") +
    xlab("Bedrock Depth [m]") +
    ylab("Residuals [m]")

  rmse_result <- round(rmse(preds, bedrock_dep, .pred)$.estimate, 1)

  p1 / p2 + plot_annotation(
    title = "Model performance from 10-fold cross-validation",
    subtitle = glue("Model metrics are RMSE = {rmse_result}")
  )
}

#' @export
create_btopos <- function(dtbs, dem) {
  dem_area <- crop(dem, dtbs)
  btopo <- dem_area - dtbs
  setNames(btopo, names(dtbs))
}

#' Extract LISA statistcs from cross-validation results
#'
#' @param res resample results tibble
#'
#' @return sf object
#' @export
extract_lisa <- function(res) {
  test_data <- map_dfr(res$splits, testing, .id = "Fold")
  preds <- map_dfr(res$predictions, function(x) x, .id = "Fold")
  preds$residual <- preds$bedrock_dep - preds$.pred

  member_points <-
    bind_cols(select(test_data, xcoords, ycoords), preds) |>
    st_as_sf(coords = c("xcoords", "ycoords"), crs = 3402)

  lisa <- member_points |>
    mutate(
      nb = st_knn(geometry, k = 7),
      wt = st_weights(nb),
      moran = local_moran(residual, nb, wt)
    ) |>
    unnest(moran) |>
    mutate(
      mean = as.character(mean),
      mean = if_else(p_ii_sim > 0.005, "Insig", mean),
      mean = factor(mean, levels = c(
        "Insig",
        "Low-High",
        "High-Low",
        "Low-Low",
        "High-High")
      )
    )
  lisa
}
