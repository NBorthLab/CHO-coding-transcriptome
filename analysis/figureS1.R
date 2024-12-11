#' Figure S1
#'
#' Additional plots of the data exploration
#'
#' Author: Markus Riedl
#' Created: May 2024
#'

library(here)
library(ggforce)
library(DESeq2)
library(ggh4x)
library(tidyverse)
library(ggpubr)

load(here("results/R/03_data-exploration/environment.RData"))

figure_themes <- function() {
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    plot.title = element_text(size = 8, vjust = -10, hjust = 0.08),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.box.margin = margin(1, 1, 1, 1),
    legend.box = "horizontal",
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(2, "mm"),
    legend.box.spacing = unit(4, "mm"),
    legend.margin = margin(4, 1, 4, 1),
    plot.margin = margin(4, 4, 4, 4)
  )
}

# ------------------------------------------------
# PCA

exploration_counts_vst$Strandedness <- recode(
  exploration_counts_vst$Strandedness,
  Forward = "Stranded",
  Reverse = "Stranded",
  Unstranded = "Unstranded"
)

endedness_pca_data <- plotPCA(exploration_counts_vst,
  intgroup = c("PairedEnd", "Strandedness"),
  returnData = TRUE
)

pca_pct_variance <- attr(endedness_pca_data, "percentVar") %>% round(2) * 100

endedness_pca <- ggplot(endedness_pca_data) +
  geom_point(
    aes(x = PC1, y = PC2, shape = Strandedness, fill = PairedEnd),
    alpha = 0.7,
    size = 1.8,
    stroke = 0.3
  ) +
  scale_shape_manual(values = c(21, 25)) +
  scale_fill_brewer(
    # start = 0,
    # end = 0.4,
    palette = "Accent",
    labels = c("Single end", "Paired end"),
    guide = guide_legend(
      override.aes = list(shape = 22, stroke = 0.3, size = 2, alpha = 1)
    )
  ) +
  labs(
    x = str_glue("PC1: {pca_pct_variance[1]}% variance"),
    y = str_glue("PC2: {pca_pct_variance[2]}% variance")
  ) +
  theme_bw() +
  figure_themes() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.position = c(0.25, 0.15),
    legend.key.height = unit(4, "mm")
  )
endedness_pca

# ------------------------------------------------
# Mean SD Plots

meansd_raw_plot <- vsn::meanSdPlot(assays(exploration_counts_SE)$counts,
  bins = 100,
  rank = TRUE,
  plot = FALSE
)$gg +
  ylim(c(0, 3e3)) +
  xlim(c(0, 1.2e4)) +
  scale_x_continuous(breaks = c(0, 5000, 10000)) +
  labs(title = "Raw counts") +
  theme_bw() +
  figure_themes() +
  theme(legend.position = "none")

meansd_vst_plot <- vsn::meanSdPlot(assay(exploration_counts_vst),
  bins = 100,
  rank = TRUE,
  plot = FALSE
)$gg +
  labs(title = "VST") +
  theme_bw() +
  figure_themes() +
  theme(legend.position = "none")

meansd_rlog_plot <- vsn::meanSdPlot(assay(exploration_counts_rlog),
  bins = 100,
  rank = TRUE,
  plot = FALSE
)$gg +
  labs(title = "rlog") +
  theme_bw() +
  figure_themes() +
  theme(legend.position = "none")

meansd_getmm_plot <- vsn::meanSdPlot(exploration_logcpm,
  bins = 100,
  rank = TRUE,
  plot = FALSE
)$gg +
  labs(title = "GeTMM + log") +
  theme_bw() +
  figure_themes() +
  theme(legend.position = "none")

# ------------------------------------------------
# PCA plots


pca_plot_vst <- make_mds_plot(pca_data_vst,
  coord_fixed = FALSE,
  point_size = 1.5,
  stroke_size = 0.2
) +
  figure_themes() +
  labs(title = "VST")
pca_plot_vst$layers[[1]] <- NULL

pca_plot_rlog <- make_mds_plot(pca_data_rlog,
  x_lim = c(-50, 50), y_lim = c(-60, 30),
  coord_fixed = FALSE,
  point_size = 1.5,
  stroke_size = 0.2
) + figure_themes() +
  labs(title = "rlog")
pca_plot_rlog$layers[[1]] <- NULL

pca_plot_logcpm <- make_mds_plot(pca_data_logcpm,
  x_lim = c(-90, 65), y_lim = c(-60, 60),
  coord_fixed = FALSE,
  point_size = 1.5,
  stroke_size = 0.2
) + figure_themes() +
  labs(title = "GeTMM + log")
pca_plot_logcpm$layers[[1]] <- NULL

percent_variance <- round(attr(pca_data_vst, "percentVar") * 100, 2)

# Combined

pca_plot_arranged <- ggarrange(pca_plot_logcpm, pca_plot_vst, pca_plot_rlog,
  nrow = 1,
  common.legend = TRUE,
  legend = "bottom"
)
pca_plot_arranged

figureS1 <- ggarrange(
  ggarrange(
    endedness_pca,
    meansd_raw_plot,
    nrow = 1,
    widths = c(0.55, 0.45),
    labels = c("a", "b")
  ),
  ggarrange(
    meansd_getmm_plot,
    meansd_vst_plot,
    meansd_rlog_plot,
    nrow = 1,
    # common.legend = TRUE,
    legend = "none"
  ),
  pca_plot_arranged,
  nrow = 3,
  heights = c(0.3, 0.23, 0.33),
  labels = c("", "c", "d")
)

ggsave(here("plots/paper/figureS1.pdf"), height = 8, width = 7)
