library(tidyverse)
library(here)
library(readxl)
library(cluster)

# two standard functional diversity indices - rao's quadratic entropy (botta-dukat 2005)
# and functional dispersion, fdis (laliberte & legendre 2010). both are weighted by each
# species' predicted presence probability at a cell, the same quantity the richness maps
# already sum, so a cell dominated by one likely species scores low on both regardless of
# how many unlikely species are also present. functional richness (convex hull volume) is
# the other common index but needs per cell geometry - intractable at 935k grid cells -
# so it is left out
norm_name <- function(x) str_replace_all(str_to_lower(x), "[^a-z0-9]+", "_")

# 01 Species list and trait table ----
fit_log <- read_csv(here("results", "fit_log.csv"), show_col_types = FALSE)
species_group <- fit_log |> filter(converged) |> select(species, group)

# the species column only holds the epithet (eg "bahianus") - taxon has the full binomial
traits_raw <- read_excel(
  here("data", "traits", "TraitCollectionFishNAtlanticNEPacificContShelf.xlsx"),
  sheet = "Trait values"
) |>
  mutate(species = norm_name(taxon))

candidate_traits <- c(
  "habitat", "feeding.mode", "tl", "body.shape", "AR", "offspring.size",
  "spawning.type", "age.maturity", "fecundity", "length.infinity",
  "growth.coefficient", "length.max", "age.max"
)

# several rows per species, one per region the trait was measured in - the first
# non missing value is kept rather than averaged, since the categorical traits
# (habitat, feeding.mode, body.shape, spawning.type) cannot be averaged at all
species_traits <- traits_raw |>
  filter(species %in% species_group$species) |>
  summarise(
    across(all_of(candidate_traits), \(x) {
      nonmissing <- x[!is.na(x)]
      if (length(nonmissing) == 0) NA else nonmissing[1]
    }),
    .by = species
  )

cat(nrow(species_traits), "of", nrow(species_group), "species matched to the trait table\n")

# 02 Gower distance and a euclidean trait space ----
# gower handles the mix of continuous and categorical traits directly, and is not
# generally euclidean, so add=TRUE applies the cailliez correction during the pcoa
# rather than assuming it already is one
trait_matrix <- species_traits |>
  column_to_rownames("species") |>
  mutate(across(c(habitat, feeding.mode, body.shape, spawning.type), factor))

gower_dist <- daisy(trait_matrix, metric = "gower")

pcoa <- cmdscale(gower_dist, k = nrow(trait_matrix) - 1, eig = TRUE, add = TRUE)
# $eig is the full n-1 spectrum, $points only has the k requested columns, so the
# keep mask has to be truncated to match before indexing
axes_kept <- pcoa$eig[seq_len(ncol(pcoa$points))] > 1e-8
coords <- pcoa$points[, axes_kept]
cat("pcoa axes kept:", ncol(coords), "of", length(pcoa$eig), "\n")

species_list <- species_traits$species
trait_dist <- as.matrix(gower_dist)[species_list, species_list]
trait_coords <- coords[species_list, ]

# 03 Species by grid cell probability matrix ----
probability <- sapply(species_list, \(sp) {
  group <- species_group$group[species_group$species == sp]
  readRDS(here("results", "predictions", group, paste0(sp, ".rds")))$mean_prob
})

# relative weight per cell - each species' share of that cell's total presence
# probability, so a cell's score depends on which species are likely there, not on
# how many species the model happens to list
weight <- probability / rowSums(probability)

# 04 Rao's quadratic entropy ----
# q = sum_a sum_b w_a w_b d_ab, the expected trait distance between two individuals
# drawn at random from the weighted community - a matrix quadratic form, computed for
# every cell at once rather than looping
rao_q <- rowSums((weight %*% trait_dist) * weight)

# 05 Functional dispersion ----
# weighted mean distance of each species to the cell's weighted centroid in the pcoa
# space. centroids for every cell come from one matrix multiply (weight %*% trait_coords);
# the species to centroid distances reuse the standard ||a-b||^2 = |a|^2 - 2a.b + |b|^2
# expansion so the cross term is also one matrix multiply, not a per cell loop
centroid <- weight %*% trait_coords
species_sq <- rowSums(trait_coords^2)
centroid_sq <- rowSums(centroid^2)
cross <- trait_coords %*% t(centroid)
dist_sq <- species_sq - 2 * cross + matrix(centroid_sq, nrow = length(species_sq), ncol = length(centroid_sq), byrow = TRUE)
dist_sq[dist_sq < 0] <- 0
species_to_centroid <- sqrt(dist_sq)
fdis <- colSums(t(weight) * species_to_centroid)

# 06 Sanity check against a toy case, not a real cell ----
# no real cell is ever dominated by exactly two species with nothing else present, so
# checking against one (as an earlier version of this script did) always "fails" even
# though the matrix form is correct - every other species still contributes a little.
# instead build a toy two-species, equal-weight community by hand, where q is known in
# closed form to be 0.5 * their trait distance, and confirm the matrix form matches it
toy_w <- c(0.5, 0.5, rep(0, ncol(weight) - 2))
toy_q <- sum((toy_w %*% trait_dist) * toy_w)
cat(
  "\ncheck: rao's q for a toy two species, equal weight community =", round(toy_q, 4),
  "| 0.5 * their trait distance =", round(0.5 * trait_dist[1, 2], 4),
  "| match:", isTRUE(all.equal(toy_q, 0.5 * trait_dist[1, 2])), "\n"
)

# 07 Assemble and save ----
grid <- read_csv(here("data", "prediction_grid.csv"), show_col_types = FALSE)

trait_diversity <- grid |>
  select(X, Y, lon, lat, year) |>
  mutate(rao_q = rao_q, fdis = fdis)

write_csv(trait_diversity, here("results", "trait_diversity_grid.csv"))

cat("\nrao's q summary:\n"); print(summary(trait_diversity$rao_q))
cat("\nfdis summary:\n"); print(summary(trait_diversity$fdis))
cat("\nwritten to results/trait_diversity_grid.csv\n")
