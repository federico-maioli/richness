library(tidyverse)
library(here)
library(sf)
library(sdmTMB)
library(sdmTMBextra)
library(rnaturalearth)

# 01 Read the prepared data ----
survey_wide <- read_csv(here("data", "survey_presence_absence.csv"), show_col_types = FALSE)

haul_vars <- c(
  "survey", "gear", "haul_id", "year", "month", "quarter", "day",
  "swept_area", "lat", "lon", "X", "Y", "depth", "temperature"
)

species_all <- setdiff(names(survey_wide), haul_vars)

# 02 Drop the species too rare to model ----
# below roughly thirty presences a binomial model has nothing to estimate a slope from
min_hauls <- 30

species_kept <- species_all[colSums(survey_wide[species_all]) >= min_hauls]

model_data <- survey_wide[c(haul_vars, species_kept)]

# 03 Drop hauls left with nothing ----
model_data <- model_data |>
  filter(rowSums(across(all_of(species_kept))) > 0)

# 04 Drop hauls missing a covariate ----
model_data <- model_data |>
  drop_na(depth, temperature) |>
  mutate(year_factor = factor(year))

# 05 Split species at the shelf break ----
# deep water barely varies in temperature, so it carries no information there and is
# dropped from that formula - see the project notes for the longer version of why
depth_cutoff <- 200

med_depth <- map_dbl(species_kept, \(sp) median(model_data$depth[model_data[[sp]] == 1], na.rm = TRUE))

shelf_species <- species_kept[med_depth <= depth_cutoff]
deep_species <- species_kept[med_depth > depth_cutoff]

cat("shelf species:", length(shelf_species), "| deep species:", length(deep_species), "\n")

# 06 Mesh, with a barrier so correlation can't cross land ----
# the north sea and skagerrak sides of jutland are close in straight line distance but
# separated by land - a plain mesh would let the spatial field leak between them
mesh <- make_mesh(model_data, xy_cols = c("X", "Y"), n_knots = 400)
cat("mesh nodes:", mesh$mesh$n, "\n")

land <- ne_countries(scale = "medium", returnclass = "sf") |>
  st_transform(3035)

# the barrier polygon has to be cropped to the survey extent, in metres since that is
# what epsg 3035 uses - the mesh itself is in km, so proj_scaling undoes that below
survey_bbox <- c(
  xmin = min(model_data$X) * 1000, ymin = min(model_data$Y) * 1000,
  xmax = max(model_data$X) * 1000, ymax = max(model_data$Y) * 1000
)
land_cropped <- suppressWarnings(st_crop(land, survey_bbox))

mesh <- add_barrier_mesh(mesh, land_cropped, range_fraction = 0.2, proj_scaling = 1000, plot = FALSE)
cat("barrier triangles:", length(mesh$barrier_triangles), "of", length(mesh$barrier_triangles) + length(mesh$normal_triangles), "\n")

# 07 Fit every species and store it ----
# spatial only - one field shared across years, with year as a random intercept rather
# than its own spatiotemporal field or a fixed factor. a fixed factor gave each of the 33
# years its own free coefficient, which blew up for any species missing from some of
# them - the same separation problem that sank the earlier hmsc attempts, just moved from
# species rarity onto year coverage. the random intercept shrinks implausible years
# instead of letting them run off to infinity. each species is independent, so one that
# fails just gets flagged and skipped rather than stalling the other 277
dir.create(here("results", "fits", "shelf"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "fits", "deep"), recursive = TRUE, showWarnings = FALSE)

fit_log <- tibble(
  species = character(), group = character(),
  elapsed_sec = numeric(), converged = logical(), error = character()
)

for (sp in shelf_species) {
  formula <- as.formula(paste(
    sp, "~ poly(log(depth), 2, raw = TRUE) + poly(temperature, 2, raw = TRUE) + (1 | year_factor)"
  ))

  elapsed <- NA_real_
  converged <- FALSE
  error_msg <- NA_character_

  t <- system.time(
    result <- tryCatch(
      sdmTMB(
        formula = formula, data = model_data, mesh = mesh,
        family = binomial(link = "logit"), spatial = "on", spatiotemporal = "off"
      ),
      error = \(e) e
    )
  )
  elapsed <- as.numeric(t["elapsed"])

  if (inherits(result, "sdmTMB")) {
    converged <- sanity(result, gradient_thresh = 0.001)$all_ok
    saveRDS(result, here("results", "fits", "shelf", paste0(sp, ".rds")))
  } else {
    error_msg <- conditionMessage(result)
  }

  fit_log <- bind_rows(fit_log, tibble(
    species = sp, group = "shelf", elapsed_sec = elapsed, converged = converged, error = error_msg
  ))
  write_csv(fit_log, here("results", "fit_log.csv"))
  cat("shelf", sp, round(elapsed, 1), "sec, converged:", converged, "\n")
}

for (sp in deep_species) {
  formula <- as.formula(paste(sp, "~ poly(log(depth), 2, raw = TRUE) + (1 | year_factor)"))

  elapsed <- NA_real_
  converged <- FALSE
  error_msg <- NA_character_

  t <- system.time(
    result <- tryCatch(
      sdmTMB(
        formula = formula, data = model_data, mesh = mesh,
        family = binomial(link = "logit"), spatial = "on", spatiotemporal = "off"
      ),
      error = \(e) e
    )
  )
  elapsed <- as.numeric(t["elapsed"])

  if (inherits(result, "sdmTMB")) {
    converged <- sanity(result, gradient_thresh = 0.001)$all_ok
    saveRDS(result, here("results", "fits", "deep", paste0(sp, ".rds")))
  } else {
    error_msg <- conditionMessage(result)
  }

  fit_log <- bind_rows(fit_log, tibble(
    species = sp, group = "deep", elapsed_sec = elapsed, converged = converged, error = error_msg
  ))
  write_csv(fit_log, here("results", "fit_log.csv"))
  cat("deep", sp, round(elapsed, 1), "sec, converged:", converged, "\n")
}

cat(
  "\ndone -", sum(fit_log$converged), "of", nrow(fit_log), "converged |",
  "total time:", round(sum(fit_log$elapsed_sec, na.rm = TRUE) / 3600, 1), "hours\n"
)
