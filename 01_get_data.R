library(tidyverse)
library(here)
library(sf)
library(terra)
library(stars)
library(askpass)
library(DATRASextra)
library(CopernicusMarine)

# read a credential from the environment, otherwise ask for it in a popup
ask_credential <- function(env_var, prompt) {
  value <- Sys.getenv(env_var)
  if (!nzchar(value)) value <- askpass(paste0(prompt, ": "))
  if (is.null(value) || !nzchar(value)) {
    stop("no ", prompt, " given, set ", env_var, " in ~/.Renviron or run interactively to be prompted")
  }
  value
}

# copernicus ships the grid as cell boundaries and the time axis as intervals counted
# in days since 1900, neither of which terra will read, so both are rewritten as a
# plain offset and step, and as plain dates
to_regular_grid <- function(s) {
  # a single depth level arrives as a length 1 axis that terra cannot take
  s <- adrop(s)
  dims <- st_dimensions(s)
  # the first two axes are the spatial ones, whatever the product calls them
  axes <- names(dims)[1:2]
  box <- st_bbox(s)
  s <- st_set_dimensions(
    s, axes[1],
    offset = box[["xmin"]], delta = (box[["xmax"]] - box[["xmin"]]) / dims[[axes[1]]]$to
  )
  s <- st_set_dimensions(
    s, axes[2],
    offset = box[["ymin"]], delta = (box[["ymax"]] - box[["ymin"]]) / dims[[axes[2]]]$to
  )
  for (ax in setdiff(names(st_dimensions(s)), axes)) {
    v <- st_get_dimension_values(s, ax)
    v <- if (inherits(v, c("Date", "POSIXct"))) as.Date(v) else as.Date(as.numeric(v), origin = "1900-01-01")
    s <- st_set_dimensions(s, ax, values = v)
  }
  s
}

# pad a bounding box
pad_bbox <- function(bbox, pad) as.numeric(bbox) + c(-pad, -pad, pad, pad)

# one year of a copernicus layer over the survey extent, a year at a time so no single
# subset request gets large enough for the service to refuse it
cms_year <- function(product, layer, variable, year, region) {
  cms_download_subset(
    product = product,
    layer = layer,
    variable = variable,
    region = region,
    timerange = c(paste0(year, "-01-01"), paste0(year, "-12-31")),
    progress = FALSE
  )
}

# 01 Set up folders ----
dir.create(here("data", "survey"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("data", "extra"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("data", "covariates"), recursive = TRUE, showWarnings = FALSE)

years_survey <- 1993:2025

# 02 Survey data ----
# one zipped exchange file per survey-year, in survey-specific subfolders
# BTS-GSA17 is the adriatic solemon survey, outside the north east atlantic
download_datras(
  path = here("data", "survey"),
  surveys = setdiff(list_surveys()$survey, c("BTS-GSA17", "Can-Mar")),
  download_ca = FALSE,
  years = years_survey,
  return_data = FALSE
)

# 03 Survey extent ----
# the climate layers are cut to the hauls rather than to a political boundary, so every
# haul sits inside the data. only HH carries positions
hh <- read_datras(here("data", "survey"), drop_hl = TRUE, drop_ca = TRUE)$HH

# padded so hauls near the edge are surrounded by data rather than sitting against it
region <- pad_bbox(
  c(
    min(hh$lon, na.rm = TRUE), min(hh$lat, na.rm = TRUE),
    max(hh$lon, na.rm = TRUE), max(hh$lat, na.rm = TRUE)
  ),
  0.5
)

# terra takes an extent as xmin, xmax, ymin, ymax
region_ext <- ext(region[c(1, 3, 2, 4)])

# 05 Depth ----
# gebco ships one 7.5 gb global netcdf inside the zip, and gdal cannot read a deflated
# member in place, so it is unpacked to the session temp folder, cut down to the survey
# extent, and thrown away rather than left sitting in a synced project folder
depth_zip <- list.files(here("data", "depth"), pattern = "\\.zip$", full.names = TRUE)

depth_nc <- unzip(
  depth_zip,
  files = grep("\\.nc$", unzip(depth_zip, list = TRUE)$Name, value = TRUE),
  exdir = tempdir(),
  junkpaths = TRUE
)

dtm <- rast(depth_nc) |> crop(region_ext)

# the dtm is a terrain model that carries land as positive elevation, so drop
# everything at or above sea level and flip the sign, leaving depth positive
depth <- ifel(dtm >= 0, NA, -dtm)
names(depth) <- "depth_m"

writeRaster(
  depth,
  here("data", "covariates", "depth.tif"),
  overwrite = TRUE,
  gdal = c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER")
)

unlink(depth_nc)

# 07 Copernicus Marine login ----
# taken from ~/.Renviron when set, otherwise asked for in a popup
cms_set_username(ask_credential("COPERNICUS_UID", "Copernicus Marine username"))
cms_set_password(ask_credential("COPERNICUS_PWD", "Copernicus Marine password"))

# cms_login returns the account details on success and NULL on failure
if (is.null(cms_login())) {
  stop("Copernicus login failed, the username or password was not accepted")
}

# the product is monthly and complete over the survey years
months <- format(
  seq(
    as.Date(paste0(min(years_survey), "-01-01")),
    as.Date(paste0(max(years_survey), "-12-01")),
    by = "month"
  ),
  "%Y-%m"
)

# 08 Bottom temperature ----
# glorys12v1 publishes potential temperature at the sea floor as a 2d monthly field,
# so there is no vertical axis to select from
temperature <- years_survey |>
  map(\(y) {
    cms_year(
      product = "GLOBAL_MULTIYEAR_PHY_001_030",
      layer = "cmems_mod_glo_phy_my_0.083deg_P1M-m",
      variable = "bottomT",
      year = y,
      region = region
    ) |>
      to_regular_grid() |>
      rast()
  }) |>
  rast()

# the subsets arrive without a crs, which would make every later crop complain
crs(temperature) <- "EPSG:4326"

stopifnot(nlyr(temperature) == length(months))
names(temperature) <- months
time(temperature) <- as.Date(paste0(months, "-15"))

writeRaster(
  temperature,
  here("data", "covariates", "temp.tif"),
  overwrite = TRUE,
  gdal = c("COMPRESS=DEFLATE")
)
