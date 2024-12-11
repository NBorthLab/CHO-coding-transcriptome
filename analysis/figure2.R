#' Figure 2
#'
#' Plot of intersections between expressed, indeterminate and non-expressed
#' genes.
#'
#' Author: Markus Riedl
#' Created: April 2024
#'

library(here)
library(RColorBrewer)
library(ComplexHeatmap)
library(grid)
library(ggpubr)
library(tidyverse)

load(here("results/R/04_expressed-genes/environment.RData"))
reactive_upset_matrix <- readRDS(
  here("results/R/indeterminate-genes/reactive_upset_matrix.rds")
)


# fraction_colors <- c(
#   rev(brewer.pal(6, "Reds")[c(3, 4, 5)]),
#   rev(brewer.pal(6, "Blues")[c(3, 4, 5)]),
#   rev(brewer.pal(6, "Greens")[c(3, 4, 5)])
# )

fraction_colors <- c(
  c("#66c2a5", "#71d8b8", "#80f3cf"), # K1
  c("#8da0cb", "#a2b9e9", "#b0c9fe"), # DXB11
  c("#c76e4d", "#e47e58", "#fc8d62") # S
)

base_fontsize <- 9

# ================================================
#   Fractions plot
# ================================================

p <- expression_fractions %>%
  ggplot(aes(x = cline, y = size, fill = interaction(gene_fraction, cline))) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  labs(x = "", y = "Number of genes") +
  scale_fill_manual(values = fraction_colors) +
  scale_x_discrete(labels = c("K1", "DXB11", "S")) +
  scale_y_continuous(
    limits = c(0, 23000),
    expand = c(0, 0.1),
    sec.axis = sec_axis(~ . / 1e4,
      breaks   = c(0.45, 1.05, 1.65),
      labels   = c("Non-\nexpressed", "Indeterminate", "Expressed")
    )
  ) +
  annotation_custom(segmentsGrob(
    x0 = unit(1.02, "npc"), x1 = unit(1.02, "npc"),
    y0 = unit(0.94, "npc"), y1 = unit(0.51, "npc")
  )) +
  annotation_custom(segmentsGrob(
    x0 = unit(1.02, "npc"), x1 = unit(1.02, "npc"),
    y0 = unit(0.43, "npc"), y1 = unit(0.48, "npc")
  )) +
  annotation_custom(segmentsGrob(
    x0 = unit(1.02, "npc"), x1 = unit(1.02, "npc"),
    y0 = unit(0.40, "npc"), y1 = unit(0.013, "npc")
  )) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = base_fontsize - 2, margin = margin(3, 0, 0, 0)),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = base_fontsize - 1),
    axis.text.y = element_text(size = base_fontsize - 2),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_text(
      size = base_fontsize - 1,
      margin = margin(0, 0, 0, 8)
    ),
    # axis.line.y.left = element_line(linewidth = 0.3),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(2, 4, 0, 10)
  )
g <- ggplotGrob(p)
g$layout$clip[g$layout$name == "panel"] <- "off"
cell_line_fractions_plot <- grid.grabExpr(grid.draw(g))

# ================================================
#   Intersection plot
# ================================================



# ------------------------------------------------
# Abandoned
#
# if (FALSE) {
#   test <- list(
#     "K1expr" = get_geneset(expressions_cell_lines, "e.."),
#     "Snonexpr1" = get_geneset(expressions_cell_lines, ".n."),
#     "DXBnonexpr1" = get_geneset(expressions_cell_lines, "..n"),
#     "Sexpr" = get_geneset(expressions_cell_lines, ".e."),
#     "K1nonexpr2" = get_geneset(expressions_cell_lines, "n.."),
#     "DXBnonexpr2" = get_geneset(expressions_cell_lines, "..n"),
#     "DXBexpr" = get_geneset(expressions_cell_lines, "..e"),
#     "K1nonexpr3" = get_geneset(expressions_cell_lines, "n.."),
#     "Snonexpr3" = get_geneset(expressions_cell_lines, ".n.")
#   )
#
#   test_comb <- make_comb_mat(test, mode = "intersect") %>%
#     `[`(,
#       c("110000000", "101000000", "111000000",
#         "000110000", "000101000", "000111000",
#         "000000110", "000000101", "000000111"
#       )
#     )
#
#   UpSet(
#     test_comb,
#     set_order = names(test),
#     comb_order = 1:9,
#     row_split = c(rep("K1", 3), rep("S", 3), rep("DXB", 3)),
#     column_split = c(rep("k1", 3), rep("s", 3), rep("d", 3)),
#     # bg_col = c(fraction_colors[4], fraction_colors[1], fraction_colors[9]),
#     left_annotation = NULL,
#     right_annotation = NULL,
#     top_annotation = upset_top_annotation(
#       test_comb,
#       add_numbers = TRUE
#     )
#   )
#
#   list_components()
#
# }
# ------------------------------------------------

# p <- tibble(
#   K1_vs_S = get_geneset_size(expressions_cell_lines, "en."),
#   K1_vs_DXB = get_geneset_size(expressions_cell_lines, "e.n"),
#   S_vs_K1 = get_geneset_size(expressions_cell_lines, "ne."),
#   S_vs_DXB = get_geneset_size(expressions_cell_lines, ".en"),
#   DXB_vs_K1 = get_geneset_size(expressions_cell_lines, "n.e"),
#   DXB_vs_S = get_geneset_size(expressions_cell_lines, ".ne"),
# ) %>%
#   pivot_longer(
#     col = everything(),
#     names_to = "comparisons",
#     values_to = "count"
#   ) %>%
#   mutate(
#     comparisons = factor(comparisons,
#       levels = c(
#         "K1_vs_S", "K1_vs_DXB", "S_vs_K1", "S_vs_DXB", "DXB_vs_K1",
#         "DXB_vs_S"
#       )
#     ),
#     expr = str_split_i(comparisons, "_", 1),
#     nonexpr = str_split_i(comparisons, "_", 3),
#     expr = factor(expr, levels = c("K1", "S", "DXB")),
#     nonexpr = factor(nonexpr, levels = c("K1", "S", "DXB"))
#   ) %>%
#   ggplot(aes(x = nonexpr, y = count)) +
#   geom_bar(stat = "identity") +
#   geom_text(aes(y = count + 6, label = count), size = 5) +
#   facet_grid(cols = vars(expr), scales = "free", switch = "x") +
#   scale_y_continuous(
#     limits = c(0, 194),
#     expand = c(0, 0.1), breaks = c(50, 100, 150)
#   ) +
#   coord_cartesian(clip = "off") +
#   theme_bw() +
#   theme(
#     # plot.margin = margin(10, 10, 30, 10),
#     axis.ticks.x = element_blank(),
#     axis.ticks = element_line(linewidth = 0.5),
#     axis.ticks.length = unit(2, "mm"),
#     strip.background = element_rect(
#       colour = "black",
#       fill = NA,
#       linewidth = 1
#     ),
#     strip.text = element_text(size = 14),
#     axis.text = element_text(size = 14),
#     panel.border = element_rect(linewidth = 0.5),
#     axis.title = element_blank(),
#     plot.margin = margin(10, 10, 10, 70)
#   ) +
#   annotation_custom(rectGrob(),
#     xmin = 0.4, xmax = 1.5, ymin = -45, ymax = -20
#   ) +
#   annotation_custom(rectGrob(),
#     xmin = 1.5, xmax = 2.6, ymin = -45, ymax = -20
#   )
# print(p)
# grid.text("Expressed",
#   x = unit(0.25, "npc"), y = unit(0.16, "npc"), just = "right"
# )
# grid.text("Non-expressed",
#   x = unit(0.25, "npc"), y = unit(0.08, "npc"), just = "right"
# )
# intersect_plot <- grid.grab()


# ComplexHeatmap Appraoch ------------------------
p <- Heatmap(
  matrix(c("K1_D", "K1_S", "D_K", "D_S", "S_K1", "S_D"), byrow = T, nrow = 1),
  column_split = factor(
    c("k", "k", "d", "d", "s", "s"),
    levels = c("k", "d", "s")
  ),
  column_title = NULL,
  top_annotation = HeatmapAnnotation(
    `Intersection size` = anno_barplot(
      c(
        get_geneset_size(expressions_cell_lines, "e.n"),
        get_geneset_size(expressions_cell_lines, "en."),
        get_geneset_size(expressions_cell_lines, "n.e"),
        get_geneset_size(expressions_cell_lines, ".ne"),
        get_geneset_size(expressions_cell_lines, "ne."),
        get_geneset_size(expressions_cell_lines, ".en")
      ),
      gp = gpar(fill = "grey40"),
      add_numbers = TRUE,
      axis_param = list(
        gp = gpar(fontsize = base_fontsize)
      ),
      numbers_gp = gpar(fontsize = base_fontsize - 2),
      numbers_rot = 0,
      height = unit(2.5, "cm")
    ),
    # blank = anno_empty(height = unit(-1, "mm"), border = FALSE),
    annotation_name_rot = 90,
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = base_fontsize - 1),
    annotation_name_offset = unit(9, "mm")
  ),
  bottom_annotation = HeatmapAnnotation(
    `Expressed` = anno_block(
      labels = c("K1", "DXB11", "S"),
      gp = gpar(fill = fraction_colors[c(1, 4, 7)]),
      height = unit(0.5, "cm"),
      show_name = TRUE,
      labels_gp = gpar(fontsize = base_fontsize - 1)
    ),
    `Non-\nexpressed` = anno_text(
      c("DXB11", "S", "K1", "S", "K1", "DXB11"),
      gp = gpar(
        border = "black",
        fill = fraction_colors[c(6, 9, 3, 9, 3, 6)],
        fontsize = base_fontsize - 4
      ),
      rot = 0,
      just = "center",
      location = 0.5,
      show_name = TRUE,
      height = unit(0.5, "cm")
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = base_fontsize - 3)
  ),
  show_heatmap_legend = FALSE,
  width = unit(40, "mm"),
  height = unit(0, "cm")
)
intersect_plot <- grid.grabExpr(draw(p))


# ================================================
#   Expressed/Indeterminate/Non-expressed
# ================================================

# Global title padding options
ht_opt$TITLE_PADDING <- unit(c(1.5, 8.5), "points")

get_upset_plot <- function(state = "e") {
  patterns <- c(
    paste0(state, ".", "."), paste0(".", state, "."), paste0(".", ".", state)
  )

  if (state == "e") {
    col <- fraction_colors[c(1, 4, 7)] %>%
      `names<-`(c("K1", "S", "DXB"))
    title <- "Expressed"
    breaks <- c(5000, 10000)
  } else if (state == "i") {
    col <- fraction_colors[c(2, 5, 8)] %>%
      `names<-`(c("K1", "S", "DXB"))
    title <- "Indeterminate"
    breaks <- c(500, 1000)
  } else if (state == "r") {
    col <- fraction_colors[c(2, 5, 8)] %>%
      `names<-`(c("K1", "S", "DXB"))
    title <- "Reactive"
    breaks <- c(500)
  } else {
    col <- fraction_colors[c(3, 6, 9)] %>%
      `names<-`(c("K1", "S", "DXB"))
    title <- "Non-expressed"
    breaks <- c(5000, 10000)
  }

  col_annot <- rowAnnotation(
    labs = anno_text(c("K1", "S", "DXB11"),
      just = "right",
      location = 0.9,
      gp = gpar(fontsize = base_fontsize - 3)
    ),
    cline = c("K1", "DXB", "S"),
    col = list(cline = col),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    simple_anno_size = unit(2, "mm")
  )

  if (state == "r") {
    upset_data <- make_comb_mat(reactive_upset_matrix)
  } else {
    upset_data <- get_geneset(expressions_cell_lines, patterns) %>%
      `names<-`(c("CHO-K1", "CHO-S", "CHO-DXB11")) %>%
      make_comb_mat()
  }

  upset_plot <- UpSet(upset_data,
    set_order = c("CHO-K1", "CHO-DXB11", "CHO-S"),
    row_names_gp = gpar(fontsize = base_fontsize - 2),
    show_row_names = FALSE,
    column_title = title,
    column_title_gp = gpar(fontsize = base_fontsize - 1),
    comb_col = "grey20",
    pt_size = unit(2, "mm"),
    lwd = 1.5,
    top_annotation = upset_top_annotation(
      upset_data,
      gp = gpar(fill = "grey40"),
      add_numbers = TRUE,
      numbers_gp = gpar(fontsize = base_fontsize - 3),
      axis_param = list(
        at = breaks
      ),
      height = unit(2, "cm"),
      annotation_name_rot = 90,
      annotation_name_gp = gpar(fontsize = base_fontsize - 1),
      annotation_name_offset = unit(1.1, "cm")
    ),
    left_annotation = col_annot,
    right_annotation = NULL,
    height = unit(10, "mm"),
    width = unit(38, "mm")
  )
  grid.grabExpr(draw(upset_plot))
}


# ================================================
#   Panels
# ================================================


upset_plots_panels <- ggarrange(
  get_upset_plot("e"),
  get_upset_plot("n"),
  get_upset_plot("i"),
  get_upset_plot("r"),
  ncol = 2,
  nrow = 2,
  labels = c("c", "d", "e", "f")
)
upset_plots_panels

sizes_panels <- ggarrange(
  ggarrange(cell_line_fractions_plot, NULL, widths = c(0.9, 0.05)),
  NULL,
  ggarrange(intersect_plot, ggplot() +
    theme_void(), widths = c(0.6, 0.05)),
  nrow = 3,
  heights = c(0.5, -0.0, 0.5),
  widths = c(0.5, 0.5),
  labels = c("a", "", "b")
  # label.y = c(1, 0, 0.95)
)
sizes_panels

if (FALSE) {
  x11(width = 7.25)
}

figure <- ggarrange(
  plotlist = list(sizes_panels, NULL, upset_plots_panels),
  # widths = c(0.5, 0.48),
  widths = c(0.35, -0.02, 0.65),
  ncol = 3
)
figure

# x11()
ggsave(here("plots/paper/figure2.pdf"),
  plot = figure,
  width = 180,
  height = 90,
  units = "mm"
)
