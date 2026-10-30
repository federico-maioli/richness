library(tidyverse)
library(here)
library(sf)
library(rnaturalearth)
library(patchwork)

ses <- read_csv(here("results", "trait_diversity_ses.csv"), show_col_types = FALSE)

# 01 Land, for context, in the same km etrs89-laea units as the grid ----
land <- ne_countries(scale = "medium", returnclass = "sf") |>
  st_transform(3035) |>
  st_crop(st_bbox(c(
    xmin = min(ses$X) * 1000, ymin = min(ses$Y) * 1000,
    xmax = max(ses$X) * 1000, ymax = max(ses$Y) * 1000
  ), crs = 3035)) |>
  st_geometry() |>
  (\(g) g / 1000)() |>
  st_set_crs(NA)

# 02 One map per index, diverging around zero since ses is a z score ----
# ses is one value per node already, not repeated across years, unlike the raw indices
limit_rao <- max(abs(ses$rao_q_ses))
limit_fdis <- max(abs(ses$fdis_ses))

map_rao <- ggplot(ses) +
  geom_tile(aes(X, Y, fill = rao_q_ses), width = 4, height = 4) +
  geom_sf(data = land, fill = "grey85", colour = NA, inherit.aes = FALSE) +
  coord_sf(xlim = range(ses$X), ylim = range(ses$Y), expand = FALSE) +
  scale_fill_distiller(palette = "RdBu", limits = c(-limit_rao, limit_rao), name = "SES") +
  labs(title = "Rao's Q, standardised effect size", x = NULL, y = NULL) +
  theme_light()

map_fdis <- ggplot(ses) +
  geom_tile(aes(X, Y, fill = fdis_ses), width = 4, height = 4) +
  geom_sf(data = land, fill = "grey85", colour = NA, inherit.aes = FALSE) +
  coord_sf(xlim = range(ses$X), ylim = range(ses$Y), expand = FALSE) +
  scale_fill_distiller(palette = "RdBu", limits = c(-limit_fdis, limit_fdis), name = "SES") +
  labs(title = "FDis, standardised effect size", x = NULL, y = NULL) +
  theme_light()

combined <- map_rao / map_fdis

ggsave(here("results", "trait_diversity_ses_maps.png"), combined, width = 8, height = 11, dpi = 200)

cat("maps written to results/trait_diversity_ses_maps.png\n")
