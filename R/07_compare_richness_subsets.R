library(tidyverse)
library(here)

# richness of two species subsets - the dryad trait compilation and the wkabsens
# sensitive species list - against the full 230 species richness already built in
# 08_predict_grid.R. summed per-species mean probability equals the mean of the richness
# draws (expectation is linear), so this reuses the saved per-species means rather than
# re-simulating, which is both fast and exactly consistent with how richness_grid.csv
# itself was built
norm_name <- function(x) str_replace_all(str_to_lower(x), "[^a-z0-9]+", "_")

shelf_species <- list.files(here("results", "predictions", "shelf"), pattern = "\\.rds$") |> str_remove("\\.rds$")
deep_species <- list.files(here("results", "predictions", "deep"), pattern = "\\.rds$") |> str_remove("\\.rds$")

species_group <- bind_rows(
  tibble(species = shelf_species, group = "shelf"),
  tibble(species = deep_species, group = "deep")
)

# 01 Read the two species lists and match them to the fitted species ----
dryad_species <- read_delim(here("data", "traits", "dryad_traits.csv"), delim = ";", show_col_types = FALSE) |>
  transmute(species = norm_name(Species)) |>
  distinct() |>
  filter(species %in% species_group$species) |>
  pull(species)

sensitive_species <- read_csv(here("data", "sensitive_species_aphia_ids.csv"), show_col_types = FALSE) |>
  transmute(species = norm_name(species)) |>
  distinct() |>
  filter(species %in% species_group$species) |>
  pull(species)

cat("dryad trait species matched to the fitted models:", length(dryad_species), "of", nrow(species_group), "\n")
cat("wkabsens sensitive species matched to the fitted models:", length(sensitive_species), "of", nrow(species_group), "\n")

# 02 Sum mean predicted probability across a species subset ----
subset_richness <- function(species_list) {
  total <- NULL
  for (sp in species_list) {
    group <- species_group$group[species_group$species == sp]
    mean_prob <- readRDS(here("results", "predictions", group, paste0(sp, ".rds")))$mean_prob
    total <- if (is.null(total)) mean_prob else total + mean_prob
  }
  total
}

dryad_richness <- subset_richness(dryad_species)
sensitive_richness <- subset_richness(sensitive_species)

# 03 Compare against the full richness already built ----
richness_grid <- read_csv(here("results", "richness_grid.csv"), show_col_types = FALSE)

comparison <- richness_grid |>
  select(X, Y, lon, lat, year, richness_full = richness_mean) |>
  mutate(
    richness_dryad = dryad_richness,
    richness_sensitive = sensitive_richness,
    pct_dryad = 100 * richness_dryad / richness_full,
    pct_sensitive = 100 * richness_sensitive / richness_full
  )

write_csv(comparison, here("results", "richness_subset_comparison.csv"))

cat("\nfull richness:\n")
print(summary(comparison$richness_full))

cat("\ndryad subset richness (", length(dryad_species), "species ):\n")
print(summary(comparison$richness_dryad))
cat("share of full richness - median:", round(median(comparison$pct_dryad), 1), "% mean:", round(mean(comparison$pct_dryad), 1), "%\n")

cat("\nsensitive subset richness (", length(sensitive_species), "species ):\n")
print(summary(comparison$richness_sensitive))
cat("share of full richness - median:", round(median(comparison$pct_sensitive), 1), "% mean:", round(mean(comparison$pct_sensitive), 1), "%\n")

cat("\ncorrelation with full richness - dryad:", round(cor(comparison$richness_full, comparison$richness_dryad), 3),
    "| sensitive:", round(cor(comparison$richness_full, comparison$richness_sensitive), 3), "\n")
