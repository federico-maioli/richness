library(tidyverse)
library(here)
library(readxl)
library(cluster)

# standardised effect size for rao's q and fdis, against a null model that shuffles
# which species carries which traits while holding each cell's predicted community
# (its richness, its evenness, everything about who is how likely to be there) exactly
# as observed. species richness and these indices are known to be positively
# correlated for a purely combinatorial reason - more species sampled means more trait
# combinations sampled, independent of any real ecological signal - and holding the
# community weights fixed in the null means that inflation shows up identically in the
# null draws too, and cancels out of the standardised score. what is left is whether
# the specific species actually present are more or less functionally spread than a
# random relabelling would be, at that same richness
set.seed(1)
n_null <- 500

norm_name <- function(x) str_replace_all(str_to_lower(x), "[^a-z0-9]+", "_")

# 01 Species list, traits, gower distance and pcoa ----
# identical to 11_trait_diversity.R's first half - these objects were not saved there
fit_log <- read_csv(here("results", "fit_log.csv"), show_col_types = FALSE)
species_group <- fit_log |> filter(converged) |> select(species, group)

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

species_traits <- traits_raw |>
  filter(species %in% species_group$species) |>
  summarise(
    across(all_of(candidate_traits), \(x) {
      nonmissing <- x[!is.na(x)]
      if (length(nonmissing) == 0) NA else nonmissing[1]
    }),
    .by = species
  )

trait_matrix <- species_traits |>
  column_to_rownames("species") |>
  mutate(across(c(habitat, feeding.mode, body.shape, spawning.type), factor))

gower_dist <- daisy(trait_matrix, metric = "gower")
pcoa <- cmdscale(gower_dist, k = nrow(trait_matrix) - 1, eig = TRUE, add = TRUE)
axes_kept <- pcoa$eig[seq_len(ncol(pcoa$points))] > 1e-8
coords <- pcoa$points[, axes_kept]

species_list <- species_traits$species
n_species <- length(species_list)
trait_dist <- as.matrix(gower_dist)[species_list, species_list]
trait_coords <- coords[species_list, ]

cat(n_species, "species |", ncol(trait_coords), "pcoa axes\n")

# 02 Species by node probability matrix, averaged across years ----
grid <- read_csv(here("data", "prediction_grid.csv"), show_col_types = FALSE) |>
  mutate(node = paste(X, Y))

nodes <- grid |> distinct(node, X, Y) |> arrange(node)
node_idx <- match(grid$node, nodes$node)
n_nodes <- nrow(nodes)

probability <- sapply(species_list, \(sp) {
  group <- species_group$group[species_group$species == sp]
  prob <- readRDS(here("results", "predictions", group, paste0(sp, ".rds")))$mean_prob
  node_sum <- rowsum(prob, node_idx)[, 1]
  node_count <- rowsum(rep(1, length(prob)), node_idx)[, 1]
  node_sum / node_count
})

weight <- probability / rowSums(probability)
richness <- rowSums(probability)

cat("weight matrix:", nrow(weight), "nodes x", ncol(weight), "species\n")

# 03 Observed rao's q and fdis ----
rao_q_of <- function(D) rowSums((weight %*% D) * weight)

fdis_of <- function(X) {
  centroid <- weight %*% X
  species_sq <- rowSums(X^2)
  centroid_sq <- rowSums(centroid^2)
  cross <- X %*% t(centroid)
  dist_sq <- species_sq - 2 * cross + matrix(centroid_sq, nrow = nrow(X), ncol = nrow(centroid), byrow = TRUE)
  dist_sq[dist_sq < 0] <- 0
  colSums(t(weight) * sqrt(dist_sq))
}

rao_q_obs <- rao_q_of(trait_dist)
fdis_obs <- fdis_of(trait_coords)

# 04 Null model - shuffle species identity against the trait data, keep weight fixed ----
rao_sum <- rao_sumsq <- numeric(n_nodes)
fdis_sum <- fdis_sumsq <- numeric(n_nodes)

t_loop <- system.time({
  for (iter in seq_len(n_null)) {
    perm <- sample(n_species)

    rao_null <- rao_q_of(trait_dist[perm, perm])
    fdis_null <- fdis_of(trait_coords[perm, ])

    rao_sum <- rao_sum + rao_null
    rao_sumsq <- rao_sumsq + rao_null^2
    fdis_sum <- fdis_sum + fdis_null
    fdis_sumsq <- fdis_sumsq + fdis_null^2

    if (iter %% 25 == 0) cat("null iteration", iter, "of", n_null, "\n")
  }
})
print(t_loop)

rao_null_mean <- rao_sum / n_null
rao_null_sd <- sqrt(pmax(rao_sumsq / n_null - rao_null_mean^2, 0))
fdis_null_mean <- fdis_sum / n_null
fdis_null_sd <- sqrt(pmax(fdis_sumsq / n_null - fdis_null_mean^2, 0))

rao_q_ses <- (rao_q_obs - rao_null_mean) / rao_null_sd
fdis_ses <- (fdis_obs - fdis_null_mean) / fdis_null_sd

# 05 Save and check whether the richness correlation survives standardising ----
result <- nodes |>
  mutate(
    richness = richness,
    rao_q = rao_q_obs, rao_q_ses = rao_q_ses,
    fdis = fdis_obs, fdis_ses = fdis_ses
  )

write_csv(result, here("results", "trait_diversity_ses.csv"))

cat("\nrao_q_ses summary:\n"); print(summary(rao_q_ses))
cat("\nfdis_ses summary:\n"); print(summary(fdis_ses))

cat(
  "\ncor(richness, rao_q):", round(cor(richness, rao_q_obs), 3),
  "-> cor(richness, rao_q_ses):", round(cor(richness, rao_q_ses), 3), "\n"
)
cat(
  "cor(richness, fdis):", round(cor(richness, fdis_obs), 3),
  "-> cor(richness, fdis_ses):", round(cor(richness, fdis_ses), 3), "\n"
)

cat("\nwritten to results/trait_diversity_ses.csv\n")
