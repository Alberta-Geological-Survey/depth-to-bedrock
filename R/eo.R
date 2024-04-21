#' Download ALOS DEM from Planetary Computer
#'
#' @return SpatRaster
#' @export
download_alos = function() {
  s_obj = stac("https://planetarycomputer.microsoft.com/api/stac/v1/")
  it_obj = s_obj |>
    stac_search(
      collections = "alos-dem",
      bbox = c(-121, 49, -108, 61.5)
    ) |>
    get_request()

  print(paste("Number of tiles:", length(it_obj$features)))

  # merge tiles
  urls = sapply(it_obj$features, function(x) paste0("/vsicurl/", x$assets$data$href))
  dem_tiles = lapply(urls, rast)
  dem = sprc(dem_tiles)
  dem_merged = merge(dem)

  # reproject to 10tm
  dem_merged = project(dem_merged, "EPSG:3402", threads = availableCores())
  writeRaster(dem_merged, here("data/alos-dem.tif"), overwrite = TRUE)
  return(dem_merged)
}

#' Download MODIS reflectance from EarthData
#'
#' @param start character for starting date to obtain MODIS scenes. In the
#'   format "YYYY-MM-DD".
#' @param end character for ending date.
#'
#' @return SpatRaster
#' @export
download_modis = function(start = "2020-07-01", end = "2020-08-30") {
  product = "MOD09A1"

  bnd = ext(c(xmin = -121, xmax = -108, ymin = 49, ymax = 61.5))
  bnd = vect(bnd)
  crs(bnd) = "EPSG:4326"

  mf = getNASA(
    product,
    start,
    end,
    aoi = bnd,
    download = TRUE,
    path = tempdir(),
    username = Sys.getenv("EARTHDATA_USER"),
    password = Sys.getenv("EARTHDATA_KEY"),
    overwrite = TRUE
  )

  from = c(1, 3, 11, 14)
  to = c(2, 3, 11, 14)
  reject = c("01,10", "1", "1", "1")
  qa_bits = cbind(from, to, reject)

  rasters = sapply(mf, function(f) {
    try({
      x = rast(f)
      quality_mask = luna::modis_mask(x[[12]], 16, qa_bits)
      x = x[[1:7]]
      mask(x, quality_mask)
    })
  })

  modis = sprc(rasters)
  modis = merge(modis)

  modis_tm = project(
    modis,
    "epsg:3402",
    threads = availableCores()
  )

  bnd_dst = project(bnd, "epsg:3402")
  modis_tm_cropped = crop(modis_tm, bnd_dst)
  writeRaster(modis_tm_cropped, here("data/modis.tif"), overwrite = TRUE)
  return(modis_tm_cropped)
}
