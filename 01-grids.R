library(here)
library(glue)
library(terra)
library(Rsagacmd)
library(future)
library(httr2)
source(here("R/eo.R"))
source(here("R/terrain.R"))

# setup environment ----
config = "edmonton500"
conf = config::get(config = config)

# prepare terrain grids ----
dem = rast(conf$dem)
reflectance = rast(here(conf$reflectance))
aoi = ext(unlist(conf$region))

# terrain analysis
dem = crop(dem, aoi)
terrain_grids = terrain_analysis(dem, conf$resolution)
xcoords = setNames(init(terrain_grids[[1]], "x"), "xcoords")
ycoords = setNames(init(terrain_grids[[1]], "y"), "ycoords")

# align reflectance ----
reflectance_aligned = resample(reflectance, terrain_grids[[1]])
reflectance_aligned = setNames(
  reflectance_aligned,
  c("l8.blue", "l8.green", "l8.red", "l8.nir", "l8.swir1", "l8.swir2")
)

# store predictors ----
stack = c(terrain_grids, reflectance_aligned, xcoords, ycoords)

writeRaster(
  stack,
  conf$predictors,
  gdal = c("BIGTIFF=IF_SAFER", "COMPRESS=DEFLATE", "COPY_SRC_OVERVIEWS=YES", "TILED=YES"),
  overwrite = TRUE
)

# store alberta boundary ----
ab_bnd = request("https://geospatial.alberta.ca/titan/rest/services") |>
  req_url_path_append("boundary/goa_administrative_area/MapServer/0/query") |>
  req_url_query(f = "geojson", where = "1=1", outFields = "*", returnGeometry = "true")

ab_bnd = vect(ab_bnd$url) |>
  project("EPSG:3402")

writeVector(ab_bnd, here("projdata/bnd-ab.gpkg"), overwrite = TRUE)
