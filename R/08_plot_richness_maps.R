library(tidyverse)
library(here)
library(sf)
library(rnaturalearth)
library(patchwork)

comparison <- read_csv(here("results", "richness_subset_comparison.csv"), show_col_types = FALSE)

# species counts for the panel titles - read fresh rather than hardcoded, so they can't
# go stale against results/richness_subset_comparison.csv the way a typed-in number would
norm_name <- function(x) str_replace_all(str_to_lower(x), "[^a-z0-9]+", "_")
fit_log <- read_csv(here("results", "fit_log.csv"), show_col_types = FALSE)
n_full <- sum(fit_log$converged)
n_dryad <- read_delim(here("data", "traits", "dryad_traits.csv"), delim = ";", show_col_types = FALSE) |>
  transmute(species = norm_name(Species)) |>
  distinct() |>
  filter(species %in% fit_log$species[fit_log$converged]) |>
  nrow()
n_sensitive <- read_csv(here("data", "sensitive_species_aphia_ids.csv"), show_col_types = FALSE) |>
  transmute(species = norm_name(species)) |>
  distinct() |>
  filter(species %in% fit_log$species[fit_log$converged]) |>
  nrow()

# 01 Average across years ----
# one spatial snapshot per node rather than eleven near identical ones - the point here
# is the spatial pattern, not its small year to year wobble
node_richness <- comparison |>
  summarise(
    richness_full = mean(richness_full),
    richness_dryad = mean(richness_dryad),
    richness_sensitive = mean(richness_sensitive),
    .by = c(X, Y)
  )

# 02 Land, for context, in the same km etrs89-laea units as the grid ----
land <- ne_countries(scale = "medium", returnclass = "sf") |>
  st_transform(3035) |>
  st_crop(st_bbox(c(
    xmin = min(node_richness$X) * 1000, ymin = min(node_richness$Y) * 1000,
    xmax = max(node_richness$X) * 1000, ymax = max(node_richness$Y) * 1000
  ), crs = 3035)) |>
  st_geometry() |>
  (\(g) g / 1000)() |>
  st_set_crs(NA)

# 03 One map per richness metric ----
map_full <- ggplot(node_richness) +
  geom_tile(aes(X, Y, fill = richness_full), width = 4, height = 4) +
  geom_sf(data = land, fill = "grey85", colour = NA, inherit.aes = FALSE) +
  coord_sf(
    xlim = range(node_richness$X), ylim = range(node_richness$Y), expand = FALSE
  ) +
  scale_fill_viridis_c(name = "richness") +
  labs(title = paste0("Full community (", n_full, " species)"), x = NULL, y = NULL) +
  theme_light()

map_dryad <- ggplot(node_richness) +
  geom_tile(aes(X, Y, fill = richness_dryad), width = 4, height = 4) +
  geom_sf(data = land, fill = "grey85", colour = NA, inherit.aes = FALSE) +
  coord_sf(
    xlim = range(node_richness$X), ylim = range(node_richness$Y), expand = FALSE
  ) +
  scale_fill_viridis_c(name = "richness") +
  labs(title = paste0("Dryad trait species (", n_dryad, " species)"), x = NULL, y = NULL) +
  theme_light()

map_sensitive <- ggplot(node_richness) +
  geom_tile(aes(X, Y, fill = richness_sensitive), width = 4, height = 4) +
  geom_sf(data = land, fill = "grey85", colour = NA, inherit.aes = FALSE) +
  coord_sf(
    xlim = range(node_richness$X), ylim = range(node_richness$Y), expand = FALSE
  ) +
  scale_fill_viridis_c(name = "richness") +
  labs(title = paste0("WKABSENS sensitive species (", n_sensitive, " species)"), x = NULL, y = NULL) +
  theme_light()

combined <- map_full / map_dryad / map_sensitive

ggsave(here("results", "richness_maps.png"), combined, width = 8, height = 16, dpi = 200)

cat("maps written to results/richness_maps.png\n")
