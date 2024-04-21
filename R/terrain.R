#' Terrain analysis
#'
#' @param dem_input SpatRaster DEM.
#' @param res resolution (in units of the DEM) used for the terrain analysis.
#'   Default is 500.
#'
#' @return SpatRaster object of terrain analysis grid.
#' @export
terrain_analysis <- function(dem_input, res = 500) {
  # create bridge to saga-gis
  saga <- Rsagacmd::saga_gis(raster_format = "SAGA", all_outputs = FALSE)

  # resample to target resolution
  dem <- saga$grid_tools$resampling(
    input = dem_input,
    target_definition = 0,
    target_user_size = res,
    output = tempfile(fileext = ".sgrd")
  )
  dem <- setNames(dem, "dem")

  # florinsky curvature
  lsp_local <- saga$ta_morphometry$slope_aspect_curvature(
    elevation = dem,
    slope = tempfile(fileext = ".sgrd"),
    aspect = tempfile(fileext = ".sgrd"),
    c_long = tempfile(fileext = ".sgrd"),
    c_cros = tempfile(fileext = ".sgrd"),
    c_tota = tempfile(fileext = ".sgrd"),
    c_roto = tempfile(fileext = ".sgrd"),
    method = 8,
    unit_slope = 1,
    unit_aspect = 1
  )
  lsp_local <- terra::rast(lsp_local) |>
    setNames(c("slope", "aspect", "c_long", "c_cros", "c_tota", "c_roto"))

  tri <- saga$ta_morphometry$terrain_ruggedness_index_tri(
    dem = dem,
    tri = tempfile(fileext = ".sgrd"),
    radius = 1,
    dw_weighting = 1
  ) |>
    setNames("tri")

  # vector ruggedness index
  vrm <- saga$ta_morphometry$vector_ruggedness_measure_vrm(
    dem = dem,
    vrm = tempfile(fileext = ".sgrd"),
    radius = 1,
    dw_weighting = 1
  ) |>
    setNames("vrm")

  # mrvbf
  mrvbf <-
    saga$ta_morphometry$multiresolution_index_of_valley_bottom_flatness_mrvbf(
      dem = dem,
      t_slope = as.numeric(Rsagacmd::mrvbf_threshold(terra::res(dem)[1])),
      p_slope = 2,
      p_pctl = 2,
      mrvbf = tempfile(fileext = ".sgrd"),
      mrrtf = tempfile(fileext = ".sgrd")
    )

  mrvbf <- terra::rast(mrvbf) |>
    setNames(c("mrvbf", "mrrtf"))

  # topographic openness
  openness <- saga$ta_lighting$topographic_openness(
    dem = dem,
    pos = tempfile(fileext = ".sgrd"),
    neg = tempfile(fileext = ".sgrd"),
    radius = terra::res(dem)[1] * 30,
    method = 2
  )

  openness <- terra::rast(openness) |>
    setNames(c("openness_pos", "openness_neg"))

  # topographic wetness index
  dem_filled = tempfile(fileext = ".sgrd")
  filled <- saga$ta_preprocessor$fill_sinks_xxl_wang_liu(
    elev = dem,
    filled = dem_filled,
    minslope = 0.01
  )

  swi <- saga$ta_hydrology$saga_wetness_index(
    dem = filled,
    twi = tempfile(fileext = ".sgrd")
  ) |>
    setNames("swi")

  # height above channels
  strahler_order = tempfile(fileext = ".sgrd")
  strahler <- saga$ta_channels$strahler_order(
    dem = filled,
    strahler = strahler_order
  )

  vdchn_func <- function(threshold, strahler, dem) {
    channels <- saga$grid_calculus$grid_calculator(
      grids = strahler,
      formula = glue::glue("ifelse(a < {threshold}, (-99999), a)"),
      result = tempfile(fileext = ".sgrd")
    )
    saga$ta_channels$vertical_distance_to_channel_network(
      elevation = dem,
      channels = channels,
      distance = tempfile(fileext = ".sgrd")
    )
  }

  strahler_thresholds <- 5:7
  vdchns <- lapply(strahler_thresholds, vdchn_func, strahler = strahler,
                   dem = dem)
  vdchns <- terra::rast(vdchns) |>
    setNames(glue::glue("vdchn{strahler_thresholds}"))

  # proximity to channels
  hdchn_func <- function(threshold, strahler, dem) {
    channels <- saga$grid_calculus$grid_calculator(
      grids = strahler,
      formula = glue::glue("ifelse(a < {threshold}, (-99999), a)"),
      result = tempfile(fileext = ".sgrd")
    )
    tryCatch({
      saga$imagery_vigra$distance_vigra(
        input = channels,
        output = tempfile(fileext = ".sgrd"),
        norm = 2
      )
    }, error = function(e) {
      saga$grid_tools$proximity_grid(
        features = channels,
        distance = tempfile(fileext = ".sgrd")
      )
    })
  }

  hdchns <- lapply(strahler_thresholds, hdchn_func, strahler = strahler, dem = dem)
  hdchns <- terra::rast(hdchns) |>
    setNames(glue::glue("hdchn{strahler_thresholds}"))

  # valley depth
  vall_depth <- saga$ta_channels$valley_depth(
    elevation = dem,
    valley_depth = tempfile(fileext = ".sgrd")
  ) |>
    setNames("vdepth")

  # terrain surface texture
  texture <- saga$ta_morphometry$terrain_surface_texture(
    dem = dem,
    texture = tempfile(fileext = ".sgrd"),
    epsilon = 0.2
  ) |>
    setNames("texture")

  # relative heights
  rsps <- lapply(vdchns, function(vh, vd) {
    saga$grid_calculus$grid_calculator(
      grids = list(vh, vd),
      formula = "g1 / (g1 + g2)",
      result = tempfile(fileext = ".sgrd")
    )
  }, vd = vall_depth)

  rsps <- terra::rast(rsps) |>
    setNames(glue::glue("rsp{strahler_thresholds}"))

  # relative heights and slope positions
  rsp2 <- saga$ta_morphometry$relative_heights_and_slope_positions(
    dem = dem,
    ho = tempfile(fileext = ".sgrd"),
    hu = tempfile(fileext = ".sgrd"),
    nh = tempfile(fileext = ".sgrd"),
    sh = tempfile(fileext = ".sgrd")
  ) |>
    setNames(c("ho", "hu", "nh", "sh"))

  # multiscale tpi
  mtpi <- saga$ta_morphometry$multi_scale_topographic_position_index_tpi(
    dem = dem,
    tpi = tempfile(fileext = ".sgrd"),
    scale_max = floor(1500 / res),
    scale_num = 5
  ) |>
    setNames("mtpi")

  # create grids stack
  stack <- c(
    dem,
    lsp_local,
    tri,
    vrm,
    mrvbf,
    openness,
    swi,
    vdchns,
    hdchns,
    vall_depth,
    texture,
    rsps,
    mtpi,
    rsp2
  )

  return(stack)
}
