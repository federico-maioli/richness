library(tidyverse)
library(here)

# the same capture efficiency question as 14_capture_efficiency_functional.R, but
# against the standardised (null model corrected) indices rather than the raw ones -
# 15_trait_diversity_null.R found that the raw indices were substantially inflated by
# richness itself (cor ~0.3-0.45), so the earlier ~90-99% capture efficiency for
# functional diversity was likely mostly that same confound. a capture value here can
# be negative, since ses itself is mostly negative (environmental filtering, not
# richness, dominates almost everywhere) - "captured" is then a share of a negative
# total, which still answers the same question: does richness ranked selection end up
# with cells whose ses is as extreme as the most extreme cells overall
ses <- read_csv(here("results", "trait_diversity_ses.csv"), show_col_types = FALSE)

# the fastest way to reach a POSITIVE total is to take the largest values first, same
# as the rarity and raw functional diversity versions of this analysis. ses sums to a
# large negative number here (almost every cell shows environmental filtering, not
# niche complementarity), so for the "optimal" curve specifically, the fastest way to
# reach that negative total is to take the most negative values first - ranking by
# descending value there would pick the least negative cells, which is the slowest
# way to reach it, not the fastest, and silently inflates "efficiency" past 100%
capture <- function(score, value, q, decreasing = TRUE) {
  top <- order(score, decreasing = decreasing)[seq_len(ceiling(q * length(score)))]
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
  optimal_direction <- sum(target) >= 0

  curve_random <- qs
  curve_rich <- sapply(qs, \(q) capture(ses$richness, target, q))
  curve_opt <- sapply(qs, \(q) capture(target, target, q, decreasing = optimal_direction))

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
  build_curves(ses$rao_q_ses, "Rao's Q SES"),
  build_curves(ses$fdis_ses, "FDis SES")
)

write_csv(result, here("results", "capture_efficiency_ses.csv"))

result_long <- result |>
  pivot_longer(c(random, richness, optimal), names_to = "selection", values_to = "captured") |>
  mutate(selection = factor(selection, levels = c("random", "richness", "optimal")))

ggplot(result_long, aes(100 * q, 100 * captured, colour = selection)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  facet_wrap(~metric, scales = "free_y") +
  scale_colour_manual(
    values = c(random = "grey60", richness = "#21918c", optimal = "#440154"),
    labels = c(random = "Random", richness = "Ranked by richness", optimal = "Ranked by the index itself (optimal)")
  ) +
  labs(
    title = "Capture efficiency of richness based cell selection, for standardised functional diversity",
    x = "Share of cells selected (%)", y = "Share of total index value captured (%)", colour = NULL
  ) +
  theme_light() +
  theme(strip.background = element_blank())

ggsave(here("results", "capture_efficiency_ses.png"), width = 11, height = 5, dpi = 200)

cat("\nwritten to results/capture_efficiency_ses.csv and results/capture_efficiency_ses.png\n")
