#' Figure 4
#'
#' Chromatin States
#'
#' Author: Markus Riedl
#' Created: April 2024
#'


library(here)
library(ComplexHeatmap)
library(tidyverse)

# x11(width = 7)

# neighbor_hmaps <- ggarrange(
#   NULL, tss_hmaps, NULL, tes_hmaps,
#   ncol = 4,
#   widths = c(-0.05, 0.5, -0.25, 0.5)
# )

# ggarrange(overlap_hmaps, neighbor_hmaps, nrow = 2)


overlap_hmaps <- ggarrange(overlap_hmap,
  nrow = 1
)
# ggsave(here("plots/paper/figure4.pdf"),
#   plot = overlap_hmaps,
#   height = 3.3
# )


# expr_hmaps <- grid.grabExpr({
#   draw(
#     (neighbor_heatmaps$expr_tss + neighbor_heatmaps$expr_tes),
#     column_title = "Expressed"
#   )
# })
# indet_hmaps <- grid.grabExpr({
#   draw(
#     (neighbor_heatmaps$igno_tss + neighbor_heatmaps$igno_tes),
#     column_title = "Indeterminate"
#   )
# })
# nonexpr_hmaps <- grid.grabExpr({
#   draw(
#     (neighbor_heatmaps$noexpr_tss + neighbor_heatmaps$noexpr_tes),
#     column_title = "Non-expressed"
#   )
# })
# neighbor_hmaps2 <- ggarrange(expr_hmaps, indet_hmaps, nonexpr_hmaps,
#   ncol = 3, nrow = 1
# )

neighbor_hmaps2 <- grid.grabExpr({
  draw(
    (neighbor_heatmaps$expr_tss + neighbor_heatmaps$expr_tes +
      neighbor_heatmaps$noexpr_tss + neighbor_heatmaps$noexpr_tes +
      neighbor_heatmaps$react_tss + neighbor_heatmaps$react_tes
    ),
    gap = unit(c(1, 2, 1, 2, 1), "mm"),
    padding = unit(c(3, 0, 0, 0), "mm")
  )
  decorate_heatmap_body("expr_tss", {
    grid.text("Expressed", 1, 1.1,
      gp = gpar(fontsize = base_fontsize, fontface = "bold")
    )
  })
  decorate_heatmap_body("react_tss", {
    grid.text("Reactive", 1, 1.1,
      gp = gpar(fontsize = base_fontsize, fontface = "bold")
    )
  })
  decorate_heatmap_body("noexpr_tss", {
    grid.text("Non-expressed", 1, 1.1,
      gp = gpar(fontsize = base_fontsize, fontface = "bold")
    )
  })
})

figure4 <- ggarrange(NULL, overlap_hmap, neighbor_hmaps2,
  nrow = 3,
  heights = c(-0.02, 0.5, 0.2),
  labels = c("", "a", "b"),
  label.y = c(0, 0.95, 1.3),
  label.x = 0.03
)

ggsave(here("plots/paper/figure4.pdf"), width = 7, height = 5.1, plot = figure4)
