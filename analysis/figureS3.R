#' Figure S3
#'
#' Upset plot for library preparation methods
#'
#' Author: Markus Riedl
#' Created: May 2024
#'

library(here)
library(RColorBrewer)
library(ComplexUpset)
library(tidyverse)

load(
  here("results/R/04_expressed-genes/environment.RData")
)

library_strat_upset_plot <- upset(
  library_strat_upset_data,
  rev(colnames(library_strat_upset_data)),
  labeller = (\(name) str_split(test, "_") %>% map(paste, collapse = " ")),
  base_annotations = list(
    "Intersection size" = intersection_size() +
      theme(axis.title.y = element_text(vjust = -40))
  ),
  set_sizes = upset_set_size() + theme(axis.text.x = element_text(
    angle = 90, vjust = 0.5, hjust = 1
  )),
  sort_sets = FALSE
)
library_strat_upset_plot

ggsave(here("plots/paper/figureS3.pdf"),
  width = 7, height = 4,
  plot = library_strat_upset_plot
)
