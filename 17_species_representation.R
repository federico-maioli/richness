library(tidyverse)
library(here)

# species representation - the standard systematic conservation planning question,
# more concrete than the abstract "rarity weighted value" the capture curves use. for
# a reserve covering q percent of the area, how many of the 221 species have at least
# half of their total predicted range inside it, under each selection strategy? and at
# a realistic reserve size, which specific species does a richness only strategy leave
# essentially unprotected

# 01 Species by node probability matrix, averaged across years ----
fit_log <- read_csv(here("results", "fit_log.csv"), show_col_types = FALSE)
species_group <- fit_log |> filter(converged) |> select(species, group)

grid <- read_csv(here("data", "prediction_grid.csv"), show_col_types = FALSE) |>
  mutate(node = paste(X, Y))

nodes <- grid |> distinct(node, X, Y) |> arrange(node)
node_idx <- match(grid$node, nodes$node)

P <- sapply(species_group$species, \(sp) {
  group <- species_group$group[species_group$species == sp]
  prob <- readRDS(here("results", "predictions", group, paste0(sp, ".rds")))$mean_prob
  node_sum <- rowsum(prob, node_idx)[, 1]
  node_count <- rowsum(rep(1, length(prob)), node_idx)[, 1]
  node_sum / node_count
})

richness <- rowSums(P)
w <- 1 / colSums(P)
wr <- as.vector(P %*% w)
total_range <- colSums(P)

cat("P:", nrow(P), "nodes x", ncol(P), "species\n")

# 02 Representation per species, for a reserve picked by a given ranking ----
# representation = the share of a species' total predicted occurrence that falls
# inside the selected cells - 1 would mean the whole predicted range is protected
represented <- function(score, q) {
  top <- order(score, decreasing = TRUE)[seq_len(ceiling(q * nrow(P)))]
  colSums(P[top, , drop = FALSE]) / total_range
}

target <- 0.5
qs <- seq(0.01, 0.5, by = 0.01)

n_represented <- function(score) {
  sapply(qs, \(q) sum(represented(score, q) >= target))
}

curve_richness <- n_represented(richness)
curve_rarity <- n_represented(wr)
# representation under random selection is q in expectation for every species, same
# reasoning as the capture curves, so this is exact rather than simulated
curve_random <- sapply(qs, \(q) sum(q >= target))

cat(
  "\nof", ncol(P), "species, number with at least", 100 * target,
  "% of their range protected, by reserve size:\n"
)
report <- tibble(q = qs, random = curve_random, richness = curve_richness, rarity = curve_rarity)
print(report |> filter(q %in% c(0.1, 0.2, 0.3)))

write_csv(report, here("results", "species_representation.csv"))

# 03 At a realistic reserve size, which species does richness based selection leave behind ----
q_focus <- 0.2
rep_richness <- represented(richness, q_focus)
rep_rarity <- represented(wr, q_focus)

left_behind <- tibble(
  species = colnames(P),
  group = species_group$group[match(colnames(P), species_group$species)],
  total_range = total_range,
  representation_richness = rep_richness,
  representation_rarity = rep_rarity
) |>
  filter(representation_richness < 0.1, representation_rarity >= target) |>
  arrange(representation_richness)

cat(
  "\nat a", 100 * q_focus, "% reserve: species with under 10% of their range protected by",
  "richness based selection, despite being adequately protected (>=", 100 * target,
  "%) by rarity based selection -", nrow(left_behind), "species:\n"
)
print(left_behind, n = 30)

write_csv(left_behind, here("results", "species_left_behind.csv"))

# 04 Plot ----
report_long <- report |>
  pivot_longer(-q, names_to = "selection", values_to = "n_species") |>
  mutate(selection = factor(selection, levels = c("random", "richness", "rarity")))

ggplot(report_long, aes(100 * q, n_species, colour = selection)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(
    values = c(random = "grey60", richness = "#3b4994", rarity = "#be64ac"),
    labels = c(random = "Random", richness = "Ranked by richness", rarity = "Ranked by rarity value")
  ) +
  labs(
    title = paste0("Species with at least ", 100 * target, "% of their range protected"),
    x = "Share of cells protected (%)", y = "Number of species (of 221)", colour = NULL
  ) +
  theme_light()

ggsave(here("results", "species_representation.png"), width = 7, height = 5, dpi = 200)

cat("\nwritten to results/species_representation.csv, species_left_behind.csv, species_representation.png\n")
