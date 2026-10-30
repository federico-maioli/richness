library(tidyverse)
library(here)
library(sf)
library(terra)
library(DATRASextra)

# 01 Read the survey data ----
# age data is not needed for presence absence, dropping it keeps the object small
survey <- read_datras(here("data", "survey"), drop_ca = TRUE)

# 02 Clean ----
survey <- clean_datras(survey)

# 03 Swept area ----
survey <- add_swept_area(survey)

# 04 Species information ----
survey <- add_species_info(survey)

# 05 Haul table ----
hauls <- as_tibble(survey$HH) |>
  transmute(
    survey = as.character(Survey),
    gear = as.character(Gear),
    haul_id = as.character(haul.id),
    year = as.integer(as.character(Year)),
    month = as.integer(Month),
    quarter = as.integer(as.character(Quarter)),
    day = as.integer(Day),
    swept_area = SweptArea,
    lat = lat,
    lon = lon,
    depth = as.numeric(Depth),
    bot_temp = as.numeric(BotTemp)
  )

# 06 Drop the baltic, the barents sea and greenland waters ----
# bits is the only survey that samples the baltic proper. its low salinity holds a
# distinct community, and pooling it with the atlantic hauls forces species that live
# on both sides - like gasterosteus_aculeatus, two thirds baltic by presence - into one
# shared depth response that fits neither population well.
# ns-idps and is-idps reach 60-76n and as far west as -46e - the barents sea and the
# denmark strait toward greenland, despite their names. together they are 339 of 95193
# hauls, too far from the rest of the shelf for the spatial field to connect them to
# anything, and etrs89-laea itself distorts badly that far from continental europe
hauls <- hauls |>
  filter(!survey %in% c("BITS", "NS-IDPS", "IS-IDPS"))

# 07 Projected coordinates in km ----
# etrs89-laea, the standard pan european equal area grid, rather than a single utm
# zone - hauls now span the nordic seas to iberia, and a utm zone badly distorts that
# far from its own meridian
hauls <- hauls |>
  filter(!is.na(lon), !is.na(lat)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  st_transform(3035) |>
  mutate(
    X = st_coordinates(geometry)[, 1] / 1000,
    Y = st_coordinates(geometry)[, 2] / 1000
  ) |>
  st_drop_geometry() |>
  relocate(X, Y, .after = lon)

# 08 Depth from the seabed map where the survey has none ----
# the survey measurement comes first, the emodnet dtm only fills the gaps
depth_map <- rast(here("data", "covariates", "depth.tif"))

hauls <- hauls |>
  mutate(depth = coalesce(depth, terra::extract(depth_map, cbind(lon, lat))[, 1]))


# 10 Temperature for the month the haul was fished ----
# the survey's own bottom temperature comes first, most hauls never recorded one
# though, so the map fills the gap - same pattern as depth above
temperature_map <- rast(here("data", "covariates", "temp.tif"))

monthly_temperature <- hauls |>
  mutate(layer = sprintf("%d-%02d", year, month)) |>
  filter(layer %in% names(temperature_map)) |>
  group_by(layer) |>
  group_modify(\(g, key) {
    tibble(
      haul_id = g$haul_id,
      temp_map = terra::extract(temperature_map[[key$layer]], cbind(g$lon, g$lat))[, 1]
    )
  }) |>
  ungroup() |>
  select(haul_id, temp_map)

hauls <- hauls |>
  left_join(monthly_temperature, by = "haul_id") |>
  mutate(temperature = coalesce(bot_temp, temp_map)) |>
  select(-bot_temp, -temp_map)

# 11 One row per haul and species ----
# the six classes covering bony fish, sharks and rays, chimaeras, sturgeons,
# lampreys and hagfish, everything else in the catch is invertebrate or algae
fish_classes <- c(
  "Teleostei", "Elasmobranchii", "Holocephali",
  "Chondrostei", "Petromyzonti", "Myxini"
)

# records identified only to genus or family would inflate richness, so they go
caught <- as_tibble(survey$HL) |>
  filter(
    class %in% fish_classes,
    rank == "Species",
    !is.na(ScientificName_WoRMS)
  ) |>
  distinct(
    haul_id = as.character(haul.id),
    species = str_replace_all(str_to_lower(ScientificName_WoRMS), "[^a-z0-9]+", "_")
  ) |>
  filter(haul_id %in% hauls$haul_id)

# 12 Presence absence in wide format ----
presence <- caught |>
  mutate(present = 1L) |>
  pivot_wider(names_from = species, values_from = present, values_fill = 0L)

# every haul keeps its row, so a species missing from a haul is a real absence
survey_wide <- hauls |>
  left_join(presence, by = "haul_id") |>
  mutate(across(-all_of(names(hauls)), \(x) replace_na(x, 0L))) |>
  arrange(survey, year, month, haul_id)

write_csv(survey_wide, here("data", "survey_presence_absence.csv"))

# 13 Spatiotemporal prediction grid ----
# a regular 4 km lattice in the same etrs89-laea coordinates as the hauls, kept only
# within max_haul_km of an actual trawl so it never extrapolates into unsurveyed water
grid_resolution_km <- 4
max_haul_km <- 20

grid <- make_survey_grid(
  x = hauls$X,
  y = hauls$Y,
  resolution = grid_resolution_km,
  max_dist = max_haul_km
)

# X and Y are kilometres while the projection counts metres, so the geometry is built
# from metre copies and the km columns are kept for the model
grid <- grid |>
  mutate(easting = X * 1000, northing = Y * 1000) |>
  st_as_sf(coords = c("easting", "northing"), crs = 3035) |>
  st_transform(4326) |>
  mutate(
    lon = st_coordinates(geometry)[, 1],
    lat = st_coordinates(geometry)[, 2],
    depth = terra::extract(depth_map, st_coordinates(geometry))[, 1]
  ) |>
  st_drop_geometry() |>
  select(X, Y, lon, lat, depth)

# 14 Drop the nodes that fall on land ----
# the seabed map already carries land as a gap, and at 15 arc seconds it draws the
# coast far more closely than a coastline polygon would
grid <- grid |>
  filter(!is.na(depth))

# 15 Drop cells shallower than the survey reaches ----
# the shallowest one percent of hauls sits at this depth, and a trawl survey barely
# reaches above it
min_depth_m <- as.numeric(quantile(hauls$depth, 0.01, na.rm = TRUE))

grid <- grid |>
  filter(depth >= min_depth_m)

# 15b Drop cells deeper than any haul ever sampled ----
# gebc's bathymetry runs far deeper than the surveys do, and a cell 20km from a shelf
# haul can still sit on the far side of the shelf break in water the model never saw -
# poly(log(depth), 2) fit on the training range diverges fast past its edge, producing
# a wedge of runaway predictions along the whole shelf break rather than a gentle
# extrapolation. the cutoff is the hard max here, not a percentile to match the shallow
# trim above - depth is heavily right skewed (p99 is 528m against a max near 1900m) so
# a percentile would cut away a real, well sampled depth range instead of just the edge
max_depth_m <- max(hauls$depth, na.rm = TRUE)

grid <- grid |>
  filter(depth <= max_depth_m)

# 16 Repeat the grid across the survey years ----
# each cell gets that year's mean bottom temperature across the twelve months, rather
# than a single month's snapshot
grid_years <- 2015:2025

grid_st <- map(grid_years, \(yr) {
  year_layers <- grep(sprintf("^%d-", yr), names(temperature_map), value = TRUE)
  annual_temperature <- mean(temperature_map[[year_layers]], na.rm = TRUE)

  # the copernicus grid is coarser than the lattice and stops short of the coast, so
  # the handful of gaps left at the edges are filled from the neighbouring cells.
  # na.policy only means cells that already hold a value are never touched
  for (pass in 1:5) {
    annual_temperature <- focal(annual_temperature, w = 5, fun = mean, na.policy = "only", na.rm = TRUE)
  }

  grid |>
    mutate(
      year = yr,
      quarter = 1L,
      temperature = terra::extract(annual_temperature, cbind(lon, lat))[, 1]
    )
}) |>
  list_rbind()

# 17 Tidy the column order ----
# year stays a real year, one the model was fitted on, so the fitted year effect is
# reused rather than drawn from the prior
grid_st <- grid_st |>
  relocate(X, Y, lon, lat, year, quarter)

# 18 Drop the handful of cells outside the training temperature range too ----
# same reasoning as the depth cutoff above, just a much smaller effect here - the
# copernicus grid barely strays past what the hauls themselves recorded
grid_st <- grid_st |>
  filter(
    temperature >= min(hauls$temperature, na.rm = TRUE),
    temperature <= max(hauls$temperature, na.rm = TRUE)
  )

write_csv(grid_st, here("data", "prediction_grid.csv"))

# 20 Cod and sole ----
# TODO
