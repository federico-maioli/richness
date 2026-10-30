library(tidyverse)
library(here)

# same question as 13_capture_efficiency.R, but asking it of functional diversity
# instead of rarity: if a reserve is picked by species richness (the easy, commonly
# used proxy), how much of the community's functional diversity does it actually
# protect, compared to picking cells for functional diversity directly? richness and
# functional diversity are not the same thing - a richness hotspot can be full of
# species that all do the same ecological job, while a less rich cell can hold a
# genuinely distinct functional mix

# 01 Richness and the two functional indices, averaged to one value per node ----
# a reserve is a place, picked once, not re-picked every year - same reasoning as the
# rarity version
richness_by_node <- read_csv(here("results", "richness_grid.csv"), show_col_types = FALSE) |>
  summarise(richness = mean(richness_mean), .by = c(X, Y))

trait_by_node <- read_csv(here("results", "trait_diversity_grid.csv"), show_col_types = FALSE) |>
  summarise(rao_q = mean(rao_q), fdis = mean(fdis), .by = c(X, Y))

nodes <- richness_by_node |> inner_join(trait_by_node, by = c("X", "Y"))
cat("nodes:", nrow(nodes), "\n")

# 02 Capture curves ----
capture <- function(score, value, q) {
  top <- order(score, decreasing = TRUE)[seq_len(ceiling(q * length(score)))]
  sum(value[top]) / sum(value)
}

auc <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

qs <- seq(0.01, 0.5, by = 0.01)

report_at <- function(curve_rich, curve_opt, q_target, label) {
  i <- which.min(abs(qs - q_target))
  cat(
    label, "- at q =", qs[i], ": richness based selection captures",
    round(100 * curve_rich[i], 1), "% of total", label, ", efficiency",
    round(100 * curve_rich[i] / curve_opt[i], 1), "% of the optimal selection\n"
  )
}

build_curves <- function(target, label) {
  curve_random <- qs
  curve_rich <- sapply(qs, \(q) capture(nodes$richness, target, q))
  curve_opt <- sapply(qs, \(q) capture(target, target, q))

  auc_random <- auc(qs, curve_random)
  cat(
    "\n", label, "- area under curve, as a multiple of random: richness",
    round(auc(qs, curve_rich) / auc_random, 3), "| optimal",
    round(auc(qs, curve_opt) / auc_random, 3),
    "| overall efficiency", round(auc(qs, curve_rich) / auc(qs, curve_opt), 3), "\n"
  )
  report_at(curve_rich, curve_opt, 0.1, label)
  report_at(curve_rich, curve_opt, 0.3, label)

  tibble(metric = label, q = qs, random = curve_random, richness = curve_rich, optimal = curve_opt)
}

result <- bind_rows(
  build_curves(nodes$rao_q, "Rao's Q"),
  build_curves(nodes$fdis, "FDis")
)

write_csv(result, here("results", "capture_efficiency_functional.csv"))

# 03 Plot ----
result_long <- result |>
  pivot_longer(c(random, richness, optimal), names_to = "selection", values_to = "captured") |>
  mutate(selection = factor(selection, levels = c("random", "richness", "optimal")))

ggplot(result_long, aes(100 * q, 100 * captured, colour = selection)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~metric) +
  scale_colour_manual(
    values = c(random = "grey60", richness = "#21918c", optimal = "#440154"),
    labels = c(random = "Random", richness = "Ranked by richness", optimal = "Ranked by the index itself (optimal)")
  ) +
  labs(
    title = "Capture efficiency of richness based cell selection, for functional diversity",
    x = "Share of cells selected (%)", y = "Share of total index value captured (%)", colour = NULL
  ) +
  theme_light() +
  theme(strip.background = element_blank())

ggsave(here("results", "capture_efficiency_functional.png"), width = 11, height = 5, dpi = 200)

cat("\nwritten to results/capture_efficiency_functional.csv and results/capture_efficiency_functional.png\n")
