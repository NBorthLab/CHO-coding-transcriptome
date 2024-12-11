#' Get gene sets of (non-)expressed and ignored genes
#'
#' Cell Lines:
#'   Order: K1, S and DXB11.
#'   Characters: <e>xpressed, <n>ot expressed and <i>gnored.
#'
#' Examples:
#'   "eee" expressed everywhere
#'   "ee\[^e\]" expressed in K1 and S, not/ignored in DXB11
#'
#'
#' @param expressions Character matrix (genes x cell lines) with characters "e",
#'      "i" and "n" indicating the expression state of the genes in cell lines.
#' @param regex_pattern Regular expression to choose the gene set.
#'
#' @return Vector of gene symbols
#'
get_geneset <- function(expressions, regex_pattern = NULL) {
  if (is.null(regex_pattern)) {
    cat("Order of phenotypes:\n")
    cat("  Cell Lines: CHO-K1, CHO-S, CHO-DXB11\n")
    return()
  }
  expr_strings <- apply(expressions, 1, paste0, collapse = "")
  if (length(regex_pattern) > 1) {
    results <- list()
    for (i in seq_along(regex_pattern)) {
      is_match <- stringr::str_detect(expr_strings, regex_pattern[i])
      results[[i]] <- rownames(expressions)[is_match]
    }
  } else {
    is_match <- stringr::str_detect(expr_strings, regex_pattern)
    results <- rownames(expressions)[is_match]
  }
  return(results)
}


#' Get size of gene sets
#'
#' Cell Lines:
#'   Order: K1, S and DXB11.
#'   Characters: <e>xpressed, <n>ot expressed and <i>gnored.
#'
#' Examples:
#'   "eee" expressed everywhere
#'   "ee\[^e\]" expressed in K1 and S, not/ignored in DXB11
#'
#'
#' @param expressions Character matrix (genes x cell lines) with characters "e",
#'      "i" and "n" indicating the expression state of the genes in cell lines.
#' @param regex_pattern Regular expression to choose the gene set.
#'
#' @return Numeric vector of gene set size(s)
#'
get_geneset_size <- function(expressions, regex_pattern = NULL) {
  if (is.null(regex_pattern)) {
    get_geneset()
  } else {
    geneset <- get_geneset(expressions, regex_pattern)
    if (is.list(geneset)) {
      map_int(geneset, length)
    } else {
      length(geneset)
    }
  }
}


#' Perform Gene Ontology Over-representation Analysis on biological processes
#'
enrich_GO_terms <- function(genes, db, universe) {
  terms_raw <- enrichGO(
    gene          = genes,
    OrgDb         = db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    minGSSize     = 4,
    universe      = universe
  )

  # If possible/applicable, make gene identifiers readable
  terms_readable <- tryCatch(
    {
      setReadable(terms_raw, OrgDb = db, keyType = "ENTREZID")
    },
    error = function(e) {
      log_warn("WARNING: Can't convert to readable!\n")
      return(terms_raw)
    }
  )

  # Simplify the GO terms to reduce redundancy of resulting terms
  terms_returned <- tryCatch(
    {
      simplify(terms_readable)
    },
    error = function(e) {
      log_warn("WARNING: Terms not simplified!")
      return(terms_readable)
    }
  )

  terms_returned@result <- terms_returned@result %>%
    mutate(log_p_adjust = -log(p.adjust, 10))

  return(terms_returned)
}


#' Perform KEGG pathway enrichment on biological processes
#'
enrich_KEGG_pathways <- function(genes, organism, organism_db, universe) {
  pathways_raw <- enrichKEGG(
    gene = genes,
    organism = organism,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe = universe,
    minGSSize = 10,
    maxGSSize = 500
  )

  pathways_readable <- tryCatch(
    {
      setReadable(pathways_raw, OrgDb = organism_db, keyType = "ENTREZID")
    },
    error = function(e) {
      log_warn("WARNING: Can't convert to readable!\n")
      return(pathways_raw)
    }
  )

  pathways_readable@result <- pathways_readable@result %>%
    mutate(log_p_adjust = -log(p.adjust, 10))

  pathways_readable@result <- pathways_readable@result %>%
    mutate(
      Description = str_remove(
        Description, " - Cricetulus griseus \\(Chinese hamster\\)"
      )
    )

  return(pathways_readable)
}


#' Mapping of gene symbols to gene ids, required for GO enrichment
map_symbol_to_geneid <- function(alias, annot = annot_db_Cg) {
  bitr(alias, fromType = "SYMBOL", toType = "ENTREZID", annot)$ENTREZID
}


#' Clean up gene alias, i.e. remove the underscore something from the gene
#' alias.
clean_up_alias <- function(dirty_alias) {
  return(stringr::str_replace(dirty_alias, "_.", ""))
}


#' Visualize KEGG pathways
#' Temporarily switch working directory so that this stupid function doesn't
#' just throw the plot in the project root.
run_pathview <- function(...) {
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(here(str_glue("results/R/{prefix}/")))
  pathview(...)
}
