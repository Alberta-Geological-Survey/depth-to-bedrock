library(here)
library(glue)
library(terra)
library(Rsagacmd)
library(future)
library(httr2)
source(here("R/eo.R"))
source(here("R/terrain.R"))

# setup environment ----
config = "edmonton100"
conf = config::get(config = config)

# prepare terrain grids ----
dem = rast(conf$dem)
modis = rast(conf$reflectance)

# terrain analysis
terrain_grids = terrain_analysis(dem, conf$resolution)
xcoords = setNames(init(terrain_grids[[1]], "x"))
ycoords = setNames(init(terrain_grids[[1]], "y"))

# align modis ----
modis_aligned = resample(modis, terrain_grids[[1]])
modis_aligned = setNames(modis_aligned, make.names(names(modis)))
names(modis_aligned) = gsub("MOD_Grid_500m_Surface_Reflectance.", "",
                            names(modis_aligned))

# store predictors ----
stack = c(terrain_grids, modis_aligned, xcoords, ycoords)

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
