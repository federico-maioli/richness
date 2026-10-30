library(tidyverse)
library(here)
library(sf)
library(rnaturalearth)
library(patchwork)

# the main figure candidate - a single map that shows, in one glance, where richness
# and rarity weighted conservation value agree and where they don't, next to the
# capture efficiency curve that quantifies the mismatch. a richness map and a rarity
# map side by side ask the reader to mentally overlay two colour scales; a bivariate
# map does that overlay for them

# teal/magenta bivariate scheme (stevens), rows = richness tercile, cols = rarity
# tercile, low to high - not fully colourblind safe (few bivariate schemes are), so the
# 3x3 legend grid doubles as a key for identifying a cell by position, not colour alone
bivariate_colours <- c(
  "1-1" = "#e8e8e8", "2-1" = "#ace4e4", "3-1" = "#5ac8c8",
  "1-2" = "#dfb0d6", "2-2" = "#a5add3", "3-2" = "#5698b9",
  "1-3" = "#be64ac", "2-3" = "#8c62aa", "3-3" = "#3b4994"
)

# 01 Classify each node into a richness x rarity tercile ----
nodes <- read_csv(here("results", "richness_rarity_by_node.csv"), show_col_types = FALSE) |>
  mutate(
    richness_tier = ntile(richness, 3),
    rarity_tier = ntile(rarity_value, 3),
    bivariate_class = paste(richness_tier, rarity_tier, sep = "-")
  )

cat("nodes per bivariate class:\n")
print(nodes |> count(bivariate_class) |> arrange(bivariate_class))

# the cells the richness story misses - high rarity value but only mid or low richness
gap_cells <- nodes |> filter(rarity_tier == 3, richness_tier < 3)
cat("\nhigh rarity value but not top richness tercile:", nrow(gap_cells), "of", nrow(nodes), "nodes\n")

# 02 Land, for context ----
land <- ne_countries(scale = "medium", returnclass = "sf") |>
  st_transform(3035) |>
  st_crop(st_bbox(c(
    xmin = min(nodes$X) * 1000, ymin = min(nodes$Y) * 1000,
    xmax = max(nodes$X) * 1000, ymax = max(nodes$Y) * 1000
  ), crs = 3035)) |>
  st_geometry() |>
  (\(g) g / 1000)() |>
  st_set_crs(NA)

# 03 The map ----
map_panel <- ggplot(nodes) +
  geom_tile(aes(X, Y, fill = bivariate_class), width = 4, height = 4) +
  geom_sf(data = land, fill = "grey95", colour = NA, inherit.aes = FALSE) +
  coord_sf(xlim = range(nodes$X), ylim = range(nodes$Y), expand = FALSE) +
  scale_fill_manual(values = bivariate_colours, guide = "none") +
  labs(title = "Where richness and rarity weighted value diverge", x = NULL, y = NULL) +
  theme_light() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# 04 The 3x3 legend, built as its own small plot and inset into the map ----
legend_data <- expand_grid(richness_tier = 1:3, rarity_tier = 1:3) |>
  mutate(bivariate_class = paste(richness_tier, rarity_tier, sep = "-"))

legend_panel <- ggplot(legend_data, aes(richness_tier, rarity_tier, fill = bivariate_class)) +
  geom_tile() +
  scale_fill_manual(values = bivariate_colours, guide = "none") +
  labs(x = "Richness →", y = "Rarity value →") +
  theme_minimal() +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title = element_text(size = 8), plot.background = element_rect(fill = "white", colour = "grey40"),
    panel.grid = element_blank()
  )

map_with_legend <- map_panel +
  inset_element(legend_panel, left = 0.02, bottom = 0.02, right = 0.22, top = 0.22)

# 05 The capture efficiency curve, as the quantitative half of the figure ----
curve_data <- read_csv(here("results", "capture_efficiency.csv"), show_col_types = FALSE) |>
  pivot_longer(c(random, richness, optimal), names_to = "selection", values_to = "captured") |>
  mutate(selection = factor(selection, levels = c("random", "richness", "optimal")))

curve_panel <- ggplot(curve_data, aes(100 * q, 100 * captured, colour = selection)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(
    values = c(random = "grey60", richness = "#3b4994", optimal = "#be64ac"),
    labels = c(random = "Random", richness = "Ranked by richness", optimal = "Ranked by rarity value (optimal)")
  ) +
  labs(
    title = "What a richness based reserve misses",
    x = "Share of cells protected (%)", y = "Share of rarity weighted value captured (%)", colour = NULL
  ) +
  theme_light() +
  theme(legend.position = "bottom", legend.direction = "vertical")

# 06 Combine ----
figure <- map_with_legend + curve_panel + plot_layout(widths = c(1.3, 1))

ggsave(here("results", "main_figure_richness_rarity.png"), figure, width = 13, height = 7, dpi = 300)

cat("\nwritten to results/main_figure_richness_rarity.png\n")
