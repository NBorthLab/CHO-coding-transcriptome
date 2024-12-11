#' Figure S2
#'
#' Normalized expressions histogram
#'
#' Author: Markus Riedl
#' Created: May 2024
#'

library(here)
library(RColorBrewer)
library(ggpubr)
library(tidyverse)

load(here("results/R/02_normalization/environment.RData"))

sizes_labels <- paste0(format(
  c(low_depth_size, middle_depth_size, high_depth_size) / 1e6,
  digits = 1
), "M")

figureS2a <- different_depths_histogram +
  geom_vline(xintercept = 0) +
  scale_fill_manual(
    labels = sizes_labels,
    values = brewer.pal(6, "YlOrRd")[4:6]
  ) +
  scale_color_manual(
    labels = sizes_labels,
    values = brewer.pal(6, "YlOrRd")[4:6]
  ) +
  labs(
    x = "GeTMM logCPM",
    y = "Count",
    color = "Library size",
    fill = "Library size"
  ) +
  theme(
    legend.position = c(0.8, 0.8),
    text = element_text(size = 9),
    legend.key.width = unit(3, "mm"),
    legend.key.height = unit(4, "mm")
    # axis.title.y = element_text(vjust = -3)
  )

figureS2b <- ggecdf(different_depths_samples,
  x = "logcpm",
  color = "sample",
  xlab = "GeTMM logCPM"
) +
  geom_vline(xintercept = 0) +
  scale_color_manual(
    labels = sizes_labels,
    values = brewer.pal(6, "YlOrRd")[4:6]
  ) +
  guides(color = "none") +
  theme(
    legend.position = c(0.8, 0.8),
    text = element_text(size = 9),
    legend.key.width = unit(3, "mm"),
    legend.key.height = unit(4, "mm")
    # axis.title.y = element_text(vjust = -2)
  )

figuresS2 <- ggarrange(figureS2a, figureS2b,
  nrow = 1, align = "h", labels = c("a", "b")
  # label.x = c(0.1, 0.08)
)
figuresS2

ggsave(
  here("plots/paper/figureS2.pdf"),
  width = 6,
  height = 3,
  plot = figuresS2
)
