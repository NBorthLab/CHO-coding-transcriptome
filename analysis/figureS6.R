#' Figure S6
#'
#' Author: Markus Riedl
#' Created: May 2024
#'


library(tidyverse)
library(here)
library(ComplexHeatmap)

load(here("results/R/06_chromatin-states/environment.RData"))

rupp_emissions <- read_tsv(
  here("resources/chromatin-states-PICR/emissions_11_all.txt")
) %>%
  rename("state" = "state (Emission order)") %>%
  column_to_rownames("state") %>%
  as.matrix() %>%
  # Make matrix columns same order as our emissions
  `[`(, colnames(emissions_matrix))

emission_heatmap_defaults <- list(
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  column_names_gp = gpar(fontsize = base_fontsize),
  row_names_gp = gpar(fontsize = base_fontsize - 1),
  row_labels = state_list,
  col = colour_mapping,
  rect_gp = gpar(col = "white", lwd = 1.5),
  border_gp = gpar(col = "black", lwd = 1),
  show_row_names = FALSE,
  show_heatmap_legend = FALSE,
  row_names_side = "left",
  # row_names_gp = gpar(fontsize = base_fontsize),
  # column_names_gp = gpar(fontsize = base_fontsize),
  left_annotation = rowAnnotation(
    r_gap = anno_empty(border = FALSE, width = gap_size)
  ),
  bottom_annotation = columnAnnotation(
    b_gap = anno_empty(border = FALSE, height = gap_size)
  ),
  width = unit(cell_size * 6, "cm"),
  height = unit(cell_size * 11, "cm")
)

emission_compare_heatmap <- local({
  heatmap_args <- utils::modifyList(
    emission_heatmap_defaults,
    list(
      matrix = rupp_emissions,
      name = "rupp",
      column_title = "PICR",
      show_row_names = TRUE
    )
  )
  picr_emissions <- do.call(Heatmap, heatmap_args)

  heatmap_args <- utils::modifyList(
    emission_heatmap_defaults,
    list(
      matrix = emissions_matrix,
      name = "emissions",
      column_title = "PICRH"
    )
  )
  picrh_emissions <- do.call(Heatmap, heatmap_args)

  grid.grabExpr(draw(picr_emissions + picrh_emissions))
})

ggsave(
  here("plots/paper/figureS6.pdf"),
  plot = emission_compare_heatmap,
  width = 4.5,
  height = 3.2
)
