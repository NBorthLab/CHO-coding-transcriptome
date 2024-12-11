#' Figure 1
#'
#' Plot for data exploration
#'
#' Author: Markus Riedl
#' Created: April 2024
#'

library(here)
library(ggforce)
library(DESeq2)
library(ggh4x)
library(tidyverse)
library(ggpubr)

load(here("results/R/03_data-exploration/environment.RData"))


# TODO: Continue here
base_fontsize <- 12


# ================================================
#   Library sizes plot
# ================================================

(library_sizes_plot <- colData(normalized_counts_SE) %>%
  as.data.frame() %>%
  mutate(
    LibraryStrategy = as.character(LibraryStrategy),
    CellLine = factor(CellLine, levels = c("CHO-S", "CHO-K1", "CHO-DXB11"))
  ) %>%
  ggplot(aes(
    x = interaction(DatasetAbbrev, LibraryStrategy, CellLine),
    y = LibrarySize,
    fill = DatasetName,
    shape = interaction(CellLine, LibraryStrategy)
  )) +
  geom_boxplot(alpha = 0.3) +
  geom_jitter(position = position_jitter(width = .12), size = 1.9) +
  scale_fill_manual(values = dataset_colors) +
  scale_shape_manual(values = libstrat_shapes) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, 1.5e8),
    breaks = scales::breaks_extended(7),
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  guides(x = "axis_nested") +
  labs(y = "Library size") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 11),
    axis.text.x = element_text(size = 12, hjust = 0.5, vjust = 1),
    axis.title.x = element_blank(),
    ggh4x.axis.nestline.x = element_line(linewidth = 1),
    ggh4x.axis.nesttext.x = element_text(
      color = "black", angle = 0, size = 9, hjust = 0.5, vjust = 0.5
    ),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.title = element_blank()
  )
)

# ================================================
#   MDS Plot function
# ================================================

make_mds_plot <- function(mds_data,
                          x_lim = c(-50, 35),
                          y_lim = c(-30, 30),
                          axis_labels = "PC",
                          coord_fixed = TRUE) {
  legend_shapes <- mds_data %>%
    mutate(shape = case_when(
      CellLine == "CHO-K1" & LibraryStrategy == "rRNA" ~ libstrat_shapes["CHO-K1.rRNA"],
      CellLine == "CHO-K1" & LibraryStrategy == "PolyA" ~ libstrat_shapes["CHO-K1.PolyA"],
      CellLine == "CHO-S" & LibraryStrategy == "PolyA" ~ libstrat_shapes["CHO-S.PolyA"],
      CellLine == "CHO-DXB11" & LibraryStrategy == "rRNA" ~ libstrat_shapes["CHO-DXB11.rRNA"]
    )) %>%
    select(DatasetName, shape) %>%
    unique() %>%
    mutate(lowercase = tolower(DatasetName)) %>%
    arrange(lowercase) %>%
    pull(shape, name = DatasetName)

  percent_variance <- round(attr(mds_data, "percentVar") * 100, 2)

  # label_pos <- mds_data %>%
  #   group_by(DatasetName) %>%
  #   summarize(
  #     x_mean = mean(component1),
  #     y_mean = mean(component2)
  #   ) %>%
  #   right_join(mds_data)

  mds_plot <- mds_data %>%
    ggplot(
      aes(
        x = component1, y = component2,
        shape = interaction(CellLine, LibraryStrategy),
        fill = DatasetName
      )
    ) +
    geom_mark_hull(
      aes(
        label = DatasetAbbrev
        # color = DatasetName
        # x0 = rep(10, 293),
        # y0 = rep(10, 293)
      ),
      alpha = .2,
      linewidth = 0.3,
      # color = NA,
      concavity = 30,
      expand = unit(2.5, "mm"),
      label.fontsize = 9,
      label.margin = margin(0.1, 0, 0, 0, "mm"),
      label.buffer = unit(0.1, "mm"),
      label.minwidth = unit(0, "mm"),
      label.width = NULL,
      label.lineheight = 1,
      con.type = "straight",
      con.size = 0.3,
      # con.colour = dataset_colors,
      con.cap = 0
    ) +
    geom_point(size = 2.5, alpha = 0.7) +
    scale_shape_manual(values = libstrat_shapes) +
    scale_fill_manual(values = dataset_colors) +
    scale_color_manual(values = dataset_colors) +
    lims(x = x_lim, y = y_lim) +
    labs(
      x = paste0(axis_labels, "1: ", percent_variance[1], "% variance"),
      y = paste0(axis_labels, "2: ", percent_variance[2], "% variance"),
      fill = "Dataset",
      shape = "Cell line/Library"
    ) +
    guides(
      shape = guide_legend(
        override.aes = list(size = 3),
        nrow = 2,
        ncol = 2
      ),
      fill = guide_legend(
        override.aes = list(shape = legend_shapes, size = 3),
        nrow = 5,
        ncol = 3
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.key.size = unit(7, "mm"),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12)
    )

  if (coord_fixed) {
    mds_plot <- mds_plot + coord_fixed()
  }

  mds_plot
}


# pca_data_vst %>%
#   group_by(DatasetName) %>%
#   summarize(
#     x_mean = mean(component1),
#     y_mean = mean(component2)
#   ) %>%
#   right_join(pca_data_vst, by = "DatasetName")

pca_plot_vst <- make_mds_plot(pca_data_vst,
  x_lim = c(-60, 45),
  y_lim = c(-38, 36),
  coord_fixed = FALSE
)
# pca_plot_vst


pcoa_plot <- make_mds_plot(pcoa_data,
  axis_labels = "leading logFC dim",
  x_lim = c(-3.1, 4),
  y_lim = c(-3.6, 3),
  coord_fixed = FALSE
)

figure_legend <- get_legend(pca_plot_vst)

label_font <- list(size = 20)

figure_mds_plots <- ggarrange(
  pca_plot_vst + rremove("legend"),
  pcoa_plot + rremove("legend"),
  align = "hv",
  labels = c("b", "c"),
  font.label = label_font
)

figure_librarysize <- ggarrange(
  library_sizes_plot,
  NULL,
  labels = c("a", ""),
  font.label = label_font
)

figure1 <- ggarrange(
  figure_librarysize,
  figure_mds_plots,
  nrow = 2,
  heights = c(.42, .5),
  align = "hv",
  common.legend = TRUE,
  legend = "none"
) +
  annotation_custom(ggplotGrob(as_ggplot(figure_legend)),
    xmin = 0.5, xmax = 1,
    ymin = 0.5, ymax = 1
  )

# x11(width = 7.24409, height = 5.70866)
# figure1
# dev.off()

ggsave(here("plots/paper/figure1.pdf"),
  plot = figure1,
  width = 184,
  height = 145,
  scale = 1.5,
  units = "mm"
)
