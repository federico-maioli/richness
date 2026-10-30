library(tidyverse)
library(here)
library(sf)

# which grid nodes already sit inside a marine protected area - the lock in constraint
# for the prioritization exercise in 19_prioritization.R. the wdpa download is the
# full global database (4.7gb of nested zips), so it is read straight out of the zip
# with a spatial filter rather than ever extracted to disk - gdal chains /vsizip/
# through a zip inside a zip just as readily as a single one

# 01 Read the three wdpa shapefile chunks, filtered to the study area ----
bbox_wkt <- "POLYGON((-18 34, 15 34, 15 63, -18 63, -18 34))"

mpa_raw <- map(0:2, \(i) {
  path <- sprintf("/vsizip/vsizip/data/mpa/WDPA_Oct2026_Public_shp.zip/WDPA_Oct2026_Public_shp_%d.zip", i)
  st_read(path, layer = "WDPA_Oct2026_Public_shp-polygons", wkt_filter = bbox_wkt, quiet = TRUE)
}) |>
  list_rbind() |>
  st_as_sf()

cat("wdpa features in the study area:", nrow(mpa_raw), "\n")

# 02 Keep sites that are both marine and actually in place ----
# gis_m_area is the area of the site's marine component - a purely terrestrial
# protected area inside the bbox (there are many) has this at zero. "proposed" sites
# are not yet legally designated, so they are not a constraint on where new protection
# can go
mpa <- mpa_raw |>
  filter(GIS_M_AREA > 0, STATUS != "Proposed") |>
  st_transform(3035) |>
  st_make_valid()

cat("marine and designated:", nrow(mpa), "\n")

# 03 Flag the grid nodes that fall inside one ----
# x and y are kilometres while epsg:3035 counts metres, same as everywhere else this
# grid gets turned into geometry
nodes <- read_csv(here("results", "richness_rarity_by_node.csv"), show_col_types = FALSE) |>
  mutate(easting = X * 1000, northing = Y * 1000) |>
  st_as_sf(coords = c("easting", "northing"), crs = 3035, remove = FALSE)

protected <- st_intersects(nodes, mpa, sparse = FALSE) |> apply(1, any)

nodes_mpa <- nodes |>
  st_drop_geometry() |>
  transmute(X, Y, already_protected = protected)

cat(
  "\nnodes already inside an mpa:", sum(nodes_mpa$already_protected), "of", nrow(nodes_mpa),
  sprintf("(%.1f%%)", 100 * mean(nodes_mpa$already_protected)), "\n"
)

write_csv(nodes_mpa, here("results", "mpa_by_node.csv"))

cat("written to results/mpa_by_node.csv\n")
