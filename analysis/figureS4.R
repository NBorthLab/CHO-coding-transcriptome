#' Figure S4
#'
#' GO Terms
#'
#' Author: Markus Riedl
#' Created: May 2024
#'

library(here)
library(clusterProfiler)
library(ggpubr)
library(tidyverse)

load(
  here("results/R/05_functional-enrichment/environment.RData")
)


common_theme <- function() {
  theme(
    text = element_text(size = 8),
    axis.text.y = element_text(size = 7), # , margin = margin(3, 3, 3, 6)),
    axis.title.x = element_text(size = 7)
  )
}

n_genes <- go_term_enrichment %>%
  `[`(c("K1_vs_S", "DXB11_vs_K1", "K1_vs_DXB11", "only_DXB11")) %>%
  map(function(x) {
    str_split(x@result$GeneRatio[1], "/")[[1]][2]
  })


ggarrange(
  go_dotplots$K1_vs_S +
    common_theme() +
    labs(title = "K1 vs S"),
  go_dotplots$DXB11_vs_K1 +
    common_theme() +
    labs(title = "DXB11 vs K1"),
  go_dotplots$K1_vs_DXB11 +
    common_theme() +
    labs(title = "K1 vs DXB11"),
  go_dotplots$only_DXB11 +
    common_theme() +
    labs(title = "only DXB11"),
  nrow = 2,
  ncol = 2,
  align = "hv",
  common.legend = TRUE
)

ggsave(here("plots/paper/figureS4.pdf"), width = 6, height = 6)
