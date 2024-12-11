#' ---
#' title: Indeterminate genes analysis
#' author: Markus Riedl
#' date: now
#' date-modified: last-modified
#' ---

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(clusterProfiler)
  library(AnnotationHub)
  library(AnnotationDbi)
  library(diptest)
  library(pathview)
  library(ComplexUpset)
  library(gghalves)
  library(ggpubr)
  library(scales)
  library(here)
  library(tidyverse)
})

prefix <- "indeterminate-genes"


attach(load_utils())


#' # Description
#'
#' We would like to investigate indeterminate a bit closer. We don't know which
#' of these genes are actually on-off switchers, depending on context, and which
#' ones are just expressed at a low level at or around the expression threshold
#' of logCPM = 0.
#'
#' # Data loading and preprocessing

log_info("Loading data")

# Load normalized expression values
normalized_counts_SE <- readRDS(
  here("results/R/02_normalization/normalized_counts_SE.rds")
)

# Load expression classes
expr_class <- readRDS(
  here("results/R/04_expressed-genes/expressions_cell_lines.rds")
)

(expr_class_df <- expr_class %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "CellLine", values_to = "class")
)

samples_metadata <- colData(normalized_counts_SE)

samples_to_cellline <- colData(normalized_counts_SE) %>%
  as_tibble() %>%
  select(SampleAccession, CellLine)

(expr_values_df <- assay(normalized_counts_SE, "normalized_logcpm") %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene, names_to = "SampleAccession", values_to = "expression"
  ) %>%
  left_join(samples_to_cellline, by = "SampleAccession")
)


(indet_expr_df <- left_join(expr_values_df, expr_class_df,
  by = c("gene", "CellLine")
) %>%
  filter(class == "i")
)

all_indet_genes <- indet_expr_df$gene %>% unique()


#' # Utility functions
#'
#' ## `plot_expression` of one gene in all cell lines

plot_expression <- function(genes) {
  indet_expr_df %>%
    filter(gene %in% {{ genes }}) %>%
    ggboxplot(
      x = "CellLine",
      y = "expression",
      add = "jitter",
      facet.by = "gene"
    )
}


#' ## `plot_top12` for top indeterminate genes

plot_top12 <- function(data) {
  top12 <- data %>%
    mutate(gene_cline = interaction(gene, CellLine)) %>%
    head(12) %>%
    pull(gene_cline)

  indet_expr_df %>%
    mutate(gene_cline = interaction(gene, CellLine)) %>%
    filter(gene_cline %in% top12) %>%
    mutate(gene_cline = factor(gene_cline, levels = top12)) %>%
    ggviolin(
      x = "gene_cline",
      y = "expression",
      add = "jitter"
    ) +
    rotate_x_text(50)
}

#' ## `plot_bottom12` for genes around the threshold

plot_bottom12 <- function(data) {
  bottom12 <- data %>%
    mutate(gene_cline = interaction(gene, CellLine)) %>%
    tail(12) %>%
    pull(gene_cline)

  indet_expr_df %>%
    mutate(gene_cline = interaction(gene, CellLine)) %>%
    filter(gene_cline %in% bottom12) %>%
    mutate(gene_cline = factor(gene_cline, levels = bottom12)) %>%
    ggboxplot(
      x = "gene_cline",
      y = "expression",
      add = "jitter"
    ) +
    rotate_x_text(50)
}

#'
#' # Ranking indeterminate genes
#'

log_info("Computing measures for indeterminate genes")

indet_expr_summary <- indet_expr_df %>%
  group_by(gene, CellLine) %>%
  summarize(
    range = max(expression) - min(expression),
    iqr = IQR(expression),
    mad = mad(expression),
    sd = sd(expression)
  )


#' ## Absolute Range
#'
#' This might be very susceptible to outliers but is able to capture also small
#' subsets of samples in contrast to other measures.

indet_by_range <- arrange(indet_expr_summary, desc(range))

plot_top12(indet_by_range)
plot_bottom12(indet_by_range)

#' Distribution of range values:

gghistogram(indet_by_range, x = "range", bins = 100)

#' How does the expression look at the 'peak' of range values?

indet_by_range %>%
  filter(range > 4.99 & range < 5.01) %>%
  {
    # subset <- .$gene
    subset <- interaction(.$gene, .$CellLine)
    indet_expr_df %>%
      filter(interaction(gene, CellLine) %in% subset)
  } %>%
  ggviolin(
    x = "gene",
    y = "expression"
  ) +
  geom_jitter(alpha = 0.2) +
  facet_grid(~CellLine, space = "free", scales = "free_x") +
  rotate_x_text(70)

#' It looks like these are just regular variable genes and not indeterminate
#' genes that turn on and off given certain conditions.
#'
#'
#' ## Inter-quartile range
#'
#' This might be favourable since this considers the distribution of the
#' expression of the genes across all samples and therefore might be less
#' problematic for outliers. However, if a small subset of all samples has a
#' different condition, to which the gene reacts, then this could not be
#' detected.
#'
#' Here, the negative control (bottom 12) shows that this measure is not useful.

indet_by_iqr <- arrange(indet_expr_summary, desc(iqr))

plot_top12(indet_by_iqr)
plot_bottom12(indet_by_iqr)


#' ## Range between Whiskers
#'
#' This is similar in that it considers the global gene expression in the
#' samples to the range and only excludes very extreme outliers.
#' A similar caveat also applies here than with the IQR.
#'
#'
#' ## Mean absolute deviation (MAD)
#'
#' Here, the negative control (bottom 12) shows that this measure is not useful.

indet_by_mad <- arrange(indet_expr_summary, desc(mad))

plot_top12(indet_by_mad)
plot_bottom12(indet_by_mad)


#' ## Standard deviation

indet_by_sd <- arrange(indet_expr_summary, desc(sd))

plot_top12(indet_by_sd)
plot_bottom12(indet_by_sd)

gghistogram(indet_by_sd, x = "sd", bins = 50)

local({
  g <- indet_by_sd %>%
    filter(sd < 2.01 & sd > 1.99)
  plot_expression(g$gene)
})


#' ## Bi/Multimodality esimation
#'
#' In theory, the interesting genes will have a bimodal distribution. Meaning
#' that there will be samples with close to zero expression and samples with a
#' fairly high expression of the gene. This yields a bimodal distribution which
#' is different from the unimodal distribution of genes that are only around the
#' expression threshold. I could use this principle to get the genes that are
#' bimodal or multimodal in their expression distribution across all samples of
#' a cell line.

log_info("Running dip test")

get_dip_pvalue <- function(gene, cell_line) {
  gene_exprs <- indet_expr_df %>%
    filter(gene == {{ gene }}, CellLine == {{ cell_line }}) %>%
    pull(expression)
  dip_result <- dip.test(gene_exprs)
  return(dip_result$p.value)
}
# get_dip_pvalue("LOC100764653", "CHO-K1")


dip_pvalues <- xfun::cache_rds(
  {
    indet_expr_df %>%
      dplyr::select(gene, CellLine) %>%
      unique() %>%
      mutate(pval = map2_vec(.$gene, .$CellLine, get_dip_pvalue))
  },
  file = "dip_pvalues",
  dir = str_glue("results/R/cache/{prefix}/"),
  hash = list(indet_expr_df)
)

dip_signif_level <- 0.05

hist(dip_pvalues$pval, breaks = 50)

reactive_genes <- dip_pvalues %>%
  mutate(
    pval_adjust = p.adjust(pval, method = "BH")
  ) %>%
  filter(pval_adjust < dip_signif_level) %>%
  mutate(
    gene_cline = paste(gene, CellLine, sep = ".")
  )


# Export reactive genes of CHO-K1 for chromatin states analysis
local({
  k1_reactive <- reactive_genes %>%
    dplyr::select(gene, CellLine) %>%
    filter(CellLine == "CHO-K1") %>%
    pull(gene)

  reactive_file <- file(
    here("results/chromatin-states/annotation/reactive-genes.txt"),
    open = "w"
  )
  writeLines(k1_reactive, con = reactive_file, sep = "\n")
  close(reactive_file)

  k1_non_reactive <- expr_class_df %>%
    filter(CellLine == "CHO-K1" & class == "i") %>%
    pull(gene) %>%
    setdiff(., k1_reactive)

  non_reactive_file <- file(
    here("results/chromatin-states/annotation/non-reactive-genes.txt"),
    open = "w"
  )
  writeLines(k1_non_reactive, con = non_reactive_file, sep = "\n")
  close(non_reactive_file)
})

#' ## Non-reactive genes

non_reactive_genes <- setdiff(
  all_indet_genes,
  reactive_genes$gene %>% unique()
)




# Visualize a random subset of the "reactive" genes.
sample_n(reactive_genes, size = 10) %>%
  {
    x <- mutate(., gene_cline = interaction(gene, CellLine))
    indet_expr_df %>%
      mutate(gene_cline = interaction(gene, CellLine)) %>%
      filter(gene_cline %in% x$gene_cline) %>%
      ggviolin(
        x = "gene",
        y = "expression"
      ) +
      geom_jitter(width = 0.2, alpha = 0.2) +
      geom_hline(yintercept = 0, alpha = 0.5) +
      facet_grid(~CellLine, scales = "free_x", space = "free") +
      rotate_x_text(70)
  }

# Visualize Genes of interest
indet_expr_df %>%
  filter(gene %in% c("Epo", "Acadm", "Isg15", "Fabp4", "Fgl2")) %>%
  ggplot(aes(x = gene, y = expression, fill = CellLine)) +
  geom_half_violin() +
  geom_half_point() +
  theme_bw()

reactive_upset_matrix <- reactive_genes %>%
  dplyr::select(gene, CellLine) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = CellLine, values_from = value, values_fill = 0) %>%
  column_to_rownames("gene")

saveobj(reactive_upset_matrix)


(reactive_genes_upset_plot <- upset(
  data = reactive_upset_matrix,
  intersect = c("CHO-K1", "CHO-S", "CHO-DXB11")
)
)

ggsave(
  here("plots/paper/plot_reactive_genes_upset_plot.pdf"),
  plot = reactive_genes_upset_plot
)

(n_reactive_genes <- reactive_genes %>%
  pull(gene) %>%
  unique() %>%
  length()
)

n_reactive_genes_cline <- reactive_genes %>%
  group_by(CellLine) %>%
  summarise(n_cline = n()) %>%
  pull(n_cline, name = CellLine)



#' # Biological context of 'interesting' indeterminate genes
#'
#' Taking only the interesting indeterminate genes into account, i.e. the genes
#' that turn on or off based on certain culture conditions or phenotypes, we can
#' ask the question if there are any biological processes enriched in those.

log_info("Starting functional enrichment")

universe <- readRDS(
  here("results/R/05_functional-enrichment/universe_Crigri.rds")
)

annotation_snapshot <- "2023-10-21"
annotation_id <- "AH114610"

annotation_hub <- AnnotationHub()
snapshotDate(annotation_hub) <- annotation_snapshot

AnnotationHub::query(
  annotation_hub,
  "Cricetulus griseus"
)[annotation_id]

annot_db_Cg <- AnnotationHub::query(
  annotation_hub,
  "Cricetulus griseus"
)[[annotation_id]]


#' ## GO Biological Process

log_info("Enriching GO terms")

reactive_genes_gobp <- enrich_GO_terms(
  map_symbol_to_geneid(reactive_genes$gene) %>% unique(),
  annot_db_Cg,
  universe
)

reactive_genes_gobp %>%
  dotplot()


#' ## KEGG Pathways

log_info("Enriching KEGG pathways")

reactive_genes_kegg <- enrich_KEGG_pathways(
  map_symbol_to_geneid(reactive_genes$gene) %>% unique(),
  "cge",
  annot_db_Cg,
  universe
)

reactive_genes_kegg_df <- reactive_genes_kegg %>% as.data.frame()

ribosome_pw_id <- reactive_genes_kegg_df %>%
  filter(Description == "Ribosome") %>%
  pull(ID)

covid_pw_id <- reactive_genes_kegg_df %>%
  filter(Description == "Coronavirus disease - COVID-19") %>%
  pull(ID)

#' ## Results

log_info("Plotting results")

dotplot_defaults <- function(p) {
  p + scale_x_continuous(
    limits = c(0, 0.15),
    expand = expansion(mult = c(0, 0.1)),
    breaks = c(0, 0.05, 0.1, 0.15)
  ) +
    labs(x = "Gene ratio of reactive genes") +
    # scale_fill_distiller(palette = "OrRd")
    scale_y_discrete(labels = label_wrap(width = 24)) +
    scale_fill_fermenter(
      palette = "Oranges",
      breaks = c(0.05, 0.01, 0.005, 0.001),
      limits = rev(c(0.0001, 0.05)),
      direction = 1,
      trans = "reverse"
    ) +
    guides(
      fill = guide_colorsteps(
        title = expression(adjusted ~ italic(p)),
        theme = theme(
          legend.text = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.ticks = element_line(),
          legend.frame = element_rect(colour = "black")
        )
      )
    ) +
    scale_size_continuous(
      limits = c(1, 90),
      breaks = c(10, 20, 40, 60, 80)
    ) +
    theme(
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      legend.direction = "horizontal",
      legend.title.position = "top",
      legend.text = element_text(size = 8, margin = margin(l = 2, t = 3)),
      legend.title = element_text(size = 8),
      legend.spacing.x = unit(0, "mm"),
      legend.key.size = unit(4, "mm"),
      legend.key.spacing.x = unit(2, "mm"),
      plot.margin = margin(l = 8, r = 8)
    )
}

gobp_dotplot <- reactive_genes_gobp %>%
  dotplot(showCategory = 8) %>%
  dotplot_defaults() +
  labs(y = "GO Biological Process")

kegg_dotplot <- reactive_genes_kegg %>%
  dotplot(showCategory = 8) %>%
  dotplot_defaults() +
  labs(y = "KEGG Pathway")

(plot_biological_context <- ggarrange(gobp_dotplot, kegg_dotplot,
  common.legend = TRUE,
  nrow = 2,
  align = "hv",
  legend = "bottom",
  labels = c("a", "b")
  # label.x = c(0.25, 0.25)
)
)

ggsave(
  here("plots/paper/figure3.pdf"),
  plot = plot_biological_context,
  width = 85,
  height = 170,
  unit = "mm"
)


(gobp_cnetplot <- reactive_genes_gobp %>%
  cnetplot(node_label = "category")
)

(kegg_cnetplot <- reactive_genes_kegg %>%
  # cnetplot(node_label = "category")
  cnetplot(
    node_label = "category",
    layout = "kk",
    showCategory = 8,
    cex.params = list(category_label = 0.8, gene_label = 0.4),
    shadowtext = "category",
    colorEdge = TRUE
  )
)

ggsave(
  here("plots/paper/figureS4.pdf"),
  plot = kegg_cnetplot,
  scale = 1.2
)

#' Visualize the pathway

ribosome_enriched_genes <- reactive_genes_kegg %>%
  as.data.frame() %>%
  filter(ID == ribosome_pw_id) %>%
  pull(geneID) %>%
  str_split("/") %>%
  `[[`(1)

ribosome_enriched_genes %>%
  bitr("SYMBOL", "ENTREZID", annot_db_Cg) %>%
  `$`(ENTREZID) %>%
  run_pathview(
    gene.data = .,
    pathway.id = ribosome_pw_id,
    species = "cge",
    kegg.dir = here(str_glue("results/R/cache/{prefix}"))
  )

# Copy the pathview figure to the paper directory
file.copy(
  here("results/R/indeterminate-genes/cge03010.pathview.png"),
  here("plots/paper/figureS5.png"),
  overwrite = TRUE
)


#' Plot the normalized expression of the ribosome enriched genes in the
#' different samples and datasets to see if there is any correlation with
#' metadata.


(ribosome_enriched_genes %>%
  bitr("SYMBOL", "GENENAME", annot_db_Cg)
)


expr_values_df %>%
  filter(gene %in% sample(ribosome_enriched_genes, 1)) %>%
  left_join(as.data.frame(samples_metadata), by = "SampleAccession") %>%
  ggplot(
    aes(
      x = DatasetReadable,
      y = expression,
      fill = Producer
    ),
    color = "black",
  ) +
  geom_jitter(shape = 21) +
  geom_violin(alpha = 0.3) +
  scale_fill_brewer(palette = "Set3") +
  facet_grid(~CellLine.y, scales = "free_x", space = "free_x") +
  rotate_x_text(90)

#' Nevermind... that's pretty useless. There does not seem to be a proper
#' correlation with anything really. And proving otherwise would be a huge pain
#' in the ass with most likely no return on investment.
#'
#'
#'
#' # Session info

rio::export(
  list(
    GO = as.data.frame(reactive_genes_gobp),
    KEGG = as.data.frame(reactive_genes_kegg)
  ),
  here("results/R/indeterminate-genes/reactive_genes_biological_context.xlsx")
)

save(
  dip_signif_level, n_reactive_genes, n_reactive_genes_cline,
  reactive_genes_kegg_df, ribosome_pw_id, covid_pw_id,
  file = here(str_glue("results/paper/{prefix}.RData"))
)

log_info("All done")

session_info()
