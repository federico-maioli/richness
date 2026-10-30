library(tidyverse)
library(here)

# does picking cells by plain species richness capture the species that are actually
# rare, or does it just find the cells with the most common species piled on top of
# each other? compares a richness-ranked selection against the best possible selection
# (ranked by the rarity weighted total itself) and against a random one, across
# reserve sizes from 1 to 50 percent of the area

# 01 Species by node probability matrix, averaged across years ----
# "cells" here are spatial nodes, not node-year rows - a reserve is a place, picked
# once, not re-picked every year
fit_log <- read_csv(here("results", "fit_log.csv"), show_col_types = FALSE)
species_group <- fit_log |> filter(converged) |> select(species, group)

grid <- read_csv(here("data", "prediction_grid.csv"), show_col_types = FALSE) |>
  mutate(node = paste(X, Y))

nodes <- grid |> distinct(node, X, Y) |> arrange(node)
node_idx <- match(grid$node, nodes$node)
n_nodes <- nrow(nodes)

# a node's depth is fixed but its temperature trim can differ by year (see
# 02_prepare_data.R section 18), so it is not always present in all eleven years -
# rowsum/count handles that properly rather than assuming a fixed divisor
P <- sapply(species_group$species, \(sp) {
  group <- species_group$group[species_group$species == sp]
  prob <- readRDS(here("results", "predictions", group, paste0(sp, ".rds")))$mean_prob
  node_sum <- rowsum(prob, node_idx)[, 1]
  node_count <- rowsum(rep(1, length(prob)), node_idx)[, 1]
  node_sum / node_count
})

cat("P:", nrow(P), "nodes x", ncol(P), "species\n")

# 02 Rarity weights ----
# the continuous equivalent of inverse range size - a species present with near
# certainty everywhere gets a small weight, one barely ever predicted present gets a
# large one, same logic as williams et al.'s range weighted richness
w <- 1 / colSums(P)

richness <- rowSums(P)
wr <- as.vector(P %*% w)

# saved per node rather than just folded into the capture curves below, so a map can
# be built from it without reloading all 221 species' predictions again
write_csv(
  nodes |> mutate(richness = richness, rarity_value = wr),
  here("results", "richness_rarity_by_node.csv")
)

# 03 Capture curves ----
capture <- function(score, wr, q) {
  top <- order(score, decreasing = TRUE)[seq_len(ceiling(q * length(score)))]
  sum(wr[top]) / sum(wr)
}

qs <- seq(0.01, 0.5, by = 0.01)
curve_random <- qs
curve_rich <- sapply(qs, \(q) capture(richness, wr, q))
curve_opt <- sapply(qs, \(q) capture(wr, wr, q))
efficiency <- curve_rich / curve_opt

# trapezoidal area under each curve, as a multiple of the random baseline's own area -
# qs sums to the same thing for every curve, so this is just each curve's average
# height over random's average height
auc <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
auc_random <- auc(qs, curve_random)

cat("\narea under curve (q = 0.01 to 0.5), as a multiple of random selection:\n")
cat("random:   ", round(auc(qs, curve_random) / auc_random, 3), "\n")
cat("richness: ", round(auc(qs, curve_rich) / auc_random, 3), "\n")
cat("optimal:  ", round(auc(qs, curve_opt) / auc_random, 3), "\n")
cat("overall efficiency (richness auc / optimal auc):", round(auc(qs, curve_rich) / auc(qs, curve_opt), 3), "\n")

# qs is a floating point sequence, so qs == 0.1 can silently miss - the nearest point
# is used instead
report_at <- function(q_target) {
  i <- which.min(abs(qs - q_target))
  cat(
    "at q =", qs[i], "(top", round(100 * qs[i]), "% of cells): richness based selection captures",
    round(100 * curve_rich[i], 1), "% of total rarity weighted value, efficiency",
    round(100 * efficiency[i], 1), "% of the optimal selection\n"
  )
}
cat("\n")
report_at(0.1)
report_at(0.3)

# 04 Save and plot ----
result <- tibble(q = qs, random = curve_random, richness = curve_rich, optimal = curve_opt, efficiency = efficiency)
write_csv(result, here("results", "capture_efficiency.csv"))

result_long <- result |>
  select(q, random, richness, optimal) |>
  pivot_longer(-q, names_to = "selection", values_to = "captured") |>
  mutate(selection = factor(selection, levels = c("random", "richness", "optimal")))

ggplot(result_long, aes(100 * q, 100 * captured, colour = selection)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(
    values = c(random = "grey60", richness = "#21918c", optimal = "#440154"),
    labels = c(random = "Random", richness = "Ranked by richness", optimal = "Ranked by rarity weighted value (optimal)")
  ) +
  labs(
    title = "Capture efficiency of richness based cell selection",
    x = "Share of cells selected (%)", y = "Share of rarity weighted value captured (%)", colour = NULL
  ) +
  theme_light()

ggsave(here("results", "capture_efficiency.png"), width = 7, height = 5, dpi = 200)

cat("\nwritten to results/capture_efficiency.csv and results/capture_efficiency.png\n")
