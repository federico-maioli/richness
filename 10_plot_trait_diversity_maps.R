library(tidyverse)
library(here)
library(sf)
library(rnaturalearth)
library(patchwork)

trait_diversity <- read_csv(here("results", "trait_diversity_grid.csv"), show_col_types = FALSE)

# 01 Average across years ----
node_diversity <- trait_diversity |>
  summarise(rao_q = mean(rao_q), fdis = mean(fdis), .by = c(X, Y))

# 02 Land, for context, in the same km etrs89-laea units as the grid ----
land <- ne_countries(scale = "medium", returnclass = "sf") |>
  st_transform(3035) |>
  st_crop(st_bbox(c(
    xmin = min(node_diversity$X) * 1000, ymin = min(node_diversity$Y) * 1000,
    xmax = max(node_diversity$X) * 1000, ymax = max(node_diversity$Y) * 1000
  ), crs = 3035)) |>
  st_geometry() |>
  (\(g) g / 1000)() |>
  st_set_crs(NA)

# 03 One map per index ----
map_rao <- ggplot(node_diversity) +
  geom_tile(aes(X, Y, fill = rao_q), width = 4, height = 4) +
  geom_sf(data = land, fill = "grey85", colour = NA, inherit.aes = FALSE) +
  coord_sf(xlim = range(node_diversity$X), ylim = range(node_diversity$Y), expand = FALSE) +
  scale_fill_viridis_c(name = "Rao's Q") +
  labs(title = "Rao's quadratic entropy", x = NULL, y = NULL) +
  theme_light()

map_fdis <- ggplot(node_diversity) +
  geom_tile(aes(X, Y, fill = fdis), width = 4, height = 4) +
  geom_sf(data = land, fill = "grey85", colour = NA, inherit.aes = FALSE) +
  coord_sf(xlim = range(node_diversity$X), ylim = range(node_diversity$Y), expand = FALSE) +
  scale_fill_viridis_c(name = "FDis") +
  labs(title = "Functional dispersion", x = NULL, y = NULL) +
  theme_light()

combined <- map_rao / map_fdis

ggsave(here("results", "trait_diversity_maps.png"), combined, width = 8, height = 11, dpi = 200)

cat("maps written to results/trait_diversity_maps.png\n")
