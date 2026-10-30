library(tidyverse)
library(here)

# schoener's d and bhattacharyya's coefficient between the spatial distributions of the
# different richness and diversity metrics - a threshold free complement to the 30%
# budget jaccard comparison in 19_prioritization.R. both treat a surface as a
# probability mass over the grid (hence the normalise to sum to 1 below), so they only
# apply to metrics that are non negative everywhere - the ses indices are mostly
# negative z scores, not a mass, so they are left out here
ses <- read_csv(here("results", "trait_diversity_ses.csv"), show_col_types = FALSE) |>
  select(X, Y, richness, rao_q, fdis)

rarity <- read_csv(here("results", "richness_rarity_by_node.csv"), show_col_types = FALSE) |>
  select(X, Y, rarity_value)

sensitive <- read_csv(here("results", "richness_subset_comparison.csv"), show_col_types = FALSE) |>
  summarise(richness_sensitive = mean(richness_sensitive), .by = c(X, Y))

nodes <- ses |>
  left_join(rarity, by = c("X", "Y")) |>
  left_join(sensitive, by = c("X", "Y"))

metrics <- c("richness", "rarity_value", "rao_q", "fdis", "richness_sensitive")

# 01 Normalise every surface to a spatial probability mass ----
p <- nodes |>
  select(all_of(metrics)) |>
  mutate(across(everything(), \(x) x / sum(x)))

# 02 Schoener's d and bhattacharyya's coefficient, every pair ----
schoeners_d <- function(a, b) 1 - 0.5 * sum(abs(a - b))
bhattacharyya <- function(a, b) sum(sqrt(a * b))

pairs <- combn(metrics, 2, simplify = FALSE)
overlap <- map_dfr(pairs, \(pr) tibble(
  metric_a = pr[1], metric_b = pr[2],
  schoeners_d = schoeners_d(p[[pr[1]]], p[[pr[2]]]),
  bhattacharyya = bhattacharyya(p[[pr[1]]], p[[pr[2]]])
))

cat("spatial distribution overlap between metrics (0 = no overlap, 1 = identical):\n")
print(overlap |> arrange(schoeners_d), n = 20)

write_csv(overlap, here("results", "distribution_overlap.csv"))

cat(
  "\nschoener's d - mean:", round(mean(overlap$schoeners_d), 3),
  "| min:", round(min(overlap$schoeners_d), 3), "| max:", round(max(overlap$schoeners_d), 3), "\n"
)
cat(
  "bhattacharyya - mean:", round(mean(overlap$bhattacharyya), 3),
  "| min:", round(min(overlap$bhattacharyya), 3), "| max:", round(max(overlap$bhattacharyya), 3), "\n"
)

# 03 Heatmaps, one per statistic ----
symmetrise <- function(d, stat) {
  bind_rows(
    d |> select(metric_a, metric_b, value = all_of(stat)),
    d |> select(metric_a = metric_b, metric_b = metric_a, value = all_of(stat)),
    tibble(metric_a = metrics, metric_b = metrics, value = 1)
  )
}

plot_heatmap <- function(mat, title) {
  ggplot(mat, aes(metric_a, metric_b, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 2)), size = 3) +
    scale_fill_viridis_c(name = title, limits = c(0, 1)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

map_d <- plot_heatmap(symmetrise(overlap, "schoeners_d"), "Schoener's D")
map_bc <- plot_heatmap(symmetrise(overlap, "bhattacharyya"), "Bhattacharyya's coefficient")

library(patchwork)
ggsave(here("results", "distribution_overlap.png"), map_d + map_bc, width = 13, height = 6, dpi = 200)

cat("\nwritten to results/distribution_overlap.csv and results/distribution_overlap.png\n")
