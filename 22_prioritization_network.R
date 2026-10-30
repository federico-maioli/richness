library(tidyverse)
library(here)
library(tidygraph)
library(ggraph)

# the same pairwise jaccard numbers as prioritization_agreement.png, as a network
# instead of a matrix. a 7x7 grid of numbers asks the reader to scan and compare cells;
# a network with weak edges dropped just shows the shape directly - one tight cluster
# and one node sitting off on its own, which is exactly the finding
agreement <- read_csv(here("results", "prioritization_agreement.csv"), show_col_types = FALSE)

# edges below this barely register against the strong ones (max is 0.97) and only add
# visual clutter - dropping them is what lets richness visibly separate from the rest
# rather than being tied to everyone by a faint, meaningless thread
edge_threshold <- 0.1

edges <- agreement |>
  filter(jaccard >= edge_threshold) |>
  rename(from = metric_a, to = metric_b, weight = jaccard)

nodes <- tibble(name = unique(c(agreement$metric_a, agreement$metric_b))) |>
  mutate(label = str_replace_all(name, "_", "\n"))

graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

set.seed(1)
ggraph(graph, layout = "fr") +
  geom_edge_link(aes(width = weight, alpha = weight), colour = "#3b4994") +
  geom_node_point(size = 16, colour = "#be64ac") +
  geom_node_text(aes(label = label), colour = "white", size = 3, lineheight = 0.85) +
  scale_edge_width(range = c(0.5, 4), guide = "none") +
  scale_edge_alpha(range = c(0.3, 0.9), name = "Jaccard\noverlap") +
  labs(title = "Which prioritization metrics agree on what to protect") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 10)))

ggsave(here("results", "prioritization_network.png"), width = 8, height = 7, dpi = 200)

cat("written to results/prioritization_network.png\n")
