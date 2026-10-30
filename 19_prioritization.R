library(tidyverse)
library(here)
library(sf)
library(rnaturalearth)
library(patchwork)

# a 30% area prioritization, with the mpas already in place locked in, repeated once
# per candidate metric to see how much the choice of metric actually moves the answer.
# a simple greedy fill (take the mpas, then the highest ranked remaining cells up to
# the budget) rather than an ilp solver - with a single budget constraint and no
# connectivity or cost term, greedy and optimal coincide here, so the extra machinery
# would not change the answer
budget_share <- 0.30

# 01 One table of every candidate metric, per node ----
mpa <- read_csv(here("results", "mpa_by_node.csv"), show_col_types = FALSE)

ses <- read_csv(here("results", "trait_diversity_ses.csv"), show_col_types = FALSE) |>
  select(X, Y, richness, rao_q, rao_q_ses, fdis, fdis_ses)

rarity <- read_csv(here("results", "richness_rarity_by_node.csv"), show_col_types = FALSE) |>
  select(X, Y, rarity_value)

sensitive <- read_csv(here("results", "richness_subset_comparison.csv"), show_col_types = FALSE) |>
  summarise(richness_sensitive = mean(richness_sensitive), .by = c(X, Y))

nodes <- mpa |>
  left_join(ses, by = c("X", "Y")) |>
  left_join(rarity, by = c("X", "Y")) |>
  left_join(sensitive, by = c("X", "Y"))

cat("nodes:", nrow(nodes), "| already protected:", sum(nodes$already_protected), "\n")

metrics <- c(
  "richness", "rarity_value", "rao_q", "fdis", "rao_q_ses", "fdis_ses", "richness_sensitive"
)

# 02 Greedy fill: lock in the mpas, then take the best remaining cells to budget ----
budget <- ceiling(budget_share * nrow(nodes))
cat("budget:", budget, "of", nrow(nodes), "nodes (", 100 * budget_share, "%)\n")

prioritise <- function(score) {
  selected <- nodes$already_protected
  remaining <- budget - sum(selected)
  if (remaining > 0) {
    candidates <- which(!selected)
    top <- candidates[order(score[candidates], decreasing = TRUE)[seq_len(remaining)]]
    selected[top] <- TRUE
  }
  selected
}

selection <- map(metrics, \(m) prioritise(nodes[[m]])) |>
  set_names(metrics) |>
  as_tibble()

solutions <- bind_cols(nodes |> select(X, Y, already_protected), selection)
write_csv(solutions, here("results", "prioritization_solutions.csv"))

# 03 How much do the metrics agree, on the cells each one chose on top of the mpas ----
# jaccard on the full solution would be inflated by the locked in cells every solution
# shares, so this compares the discretionary choices only
discretionary <- selection |> mutate(across(everything(), \(x) x & !nodes$already_protected))

jaccard <- function(a, b) sum(a & b) / sum(a | b)

pairs <- combn(metrics, 2, simplify = FALSE)
agreement <- map_dfr(pairs, \(p) tibble(
  metric_a = p[1], metric_b = p[2],
  jaccard = jaccard(discretionary[[p[1]]], discretionary[[p[2]]])
))

cat("\npairwise jaccard overlap, discretionary cells only (0 = no overlap, 1 = identical):\n")
print(agreement |> arrange(jaccard))

write_csv(agreement, here("results", "prioritization_agreement.csv"))

cat(
  "\nmean pairwise jaccard:", round(mean(agreement$jaccard), 3),
  "| min:", round(min(agreement$jaccard), 3), "| max:", round(max(agreement$jaccard), 3), "\n"
)

# 04 Selection frequency map - how many of the 7 metrics choose to ADD each cell ----
# every metric trivially "agrees" on the mpas already in place, since those are locked
# in regardless of ranking - counting them as agreement would bury the real result
# (almost no overlap on what gets added) under the one thing every solution shares
nodes_freq <- nodes |>
  mutate(
    n_metrics = rowSums(discretionary),
    status = case_when(
      already_protected ~ "Already protected",
      n_metrics == length(metrics) ~ "Chosen by every metric",
      n_metrics == 0 ~ "Chosen by none",
      TRUE ~ "Metric dependent"
    )
  )

cat("\ncells added by every metric:", sum(nodes_freq$n_metrics == length(metrics)), "\n")
cat("cells added by at least one but not all:", sum(nodes_freq$n_metrics > 0 & nodes_freq$n_metrics < length(metrics)), "\n")

land <- ne_countries(scale = "medium", returnclass = "sf") |>
  st_transform(3035) |>
  st_crop(st_bbox(c(
    xmin = min(nodes$X) * 1000, ymin = min(nodes$Y) * 1000,
    xmax = max(nodes$X) * 1000, ymax = max(nodes$Y) * 1000
  ), crs = 3035)) |>
  st_geometry() |>
  (\(g) g / 1000)() |>
  st_set_crs(NA)

# already protected cells get a flat neutral tile first - they are not a metric
# agreeing on anything, they were never a choice - then candidate cells are coloured
# by how many metrics would have added them, so the two kinds of cell read as visually
# distinct rather than blending into one misleading scale
map_freq <- ggplot(nodes_freq) +
  geom_tile(
    data = nodes_freq |> filter(already_protected),
    aes(X, Y), fill = "grey65", width = 4, height = 4
  ) +
  geom_tile(
    data = nodes_freq |> filter(!already_protected),
    aes(X, Y, fill = n_metrics), width = 4, height = 4
  ) +
  geom_sf(data = land, fill = "grey85", colour = NA, inherit.aes = FALSE) +
  coord_sf(xlim = range(nodes$X), ylim = range(nodes$Y), expand = FALSE) +
  scale_fill_viridis_c(name = "metrics\nchoosing\nto add this\ncell (grey =\nalready\nprotected)", breaks = c(0, 7)) +
  labs(title = "How many of the 7 metrics select each cell", x = NULL, y = NULL) +
  theme_light()

ggsave(here("results", "prioritization_selection_frequency.png"), map_freq, width = 7, height = 9, dpi = 200)

# 05 Heatmap of the pairwise agreement ----
agreement_matrix <- bind_rows(
  agreement,
  agreement |> rename(metric_a = metric_b, metric_b = metric_a),
  tibble(metric_a = metrics, metric_b = metrics, jaccard = 1)
)

map_agreement <- ggplot(agreement_matrix, aes(metric_a, metric_b, fill = jaccard)) +
  geom_tile() +
  geom_text(aes(label = round(jaccard, 2)), size = 3) +
  scale_fill_viridis_c(name = "Jaccard", limits = c(0, 1)) +
  labs(title = "Agreement between metrics on which cells to add", x = NULL, y = NULL) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here("results", "prioritization_agreement.png"), map_agreement, width = 7, height = 6, dpi = 200)

cat("\nwritten to results/prioritization_solutions.csv, prioritization_agreement.csv,\n")
cat("prioritization_selection_frequency.png, prioritization_agreement.png\n")
