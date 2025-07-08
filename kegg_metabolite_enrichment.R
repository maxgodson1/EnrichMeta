#' Perform KEGG metabolite enrichment
#'
#' @param KEGGid A character vector of KEGG compound IDs (e.g., "C00160")
#' @param species A character string, KEGG species code (e.g., "hsa", "eco")
#' @param p.adjust.method Adjustment method for multiple testing, default is "BH"
#'
#' @return A data.frame with enrichment results
#' @export
kegg_metabolite_enrichment <- function(KEGGid, species, p.adjust.method = "BH") {
  pathways <- KEGGREST::keggList("pathway", species)
  path_ids <- sub("path:", "", names(pathways))

  pathway2compound <- list()
  for (pid in path_ids) {
    entry <- tryCatch(KEGGREST::keggGet(pid)[[1]], error = function(e) NULL)
    if (!is.null(entry) && !is.null(entry$COMPOUND)) {
      compounds <- names(entry$COMPOUND)
      pathway2compound[[pid]] <- compounds
    }
  }

  all_cmpds <- unique(unlist(pathway2compound))

  results <- data.frame(
    ID = character(), Description = character(),
    GeneRatio = character(), BgRatio = character(),
    pvalue = numeric(), EnrichmentRatio = numeric(),
    geneID = character(), Count = integer(),
    stringsAsFactors = FALSE
  )

  N <- length(KEGGid)
  M <- length(all_cmpds)

  for (pid in names(pathway2compound)) {
    pw_cmpds <- pathway2compound[[pid]]
    overlap <- intersect(KEGGid, pw_cmpds)
    k <- length(overlap)
    n <- length(pw_cmpds)
    if (k > 0) {
      p <- phyper(k - 1, n, M - n, N, lower.tail = FALSE)
      enrich_ratio <- (k * M) / (n * N)
      matched <- paste(overlap, collapse = "; ")
      results <- rbind(results, data.frame(
        ID = pid,
        Description = pathways[pid],
        GeneRatio = paste0(k, "/", N),
        BgRatio = paste0(n, "/", M),
        pvalue = p,
        EnrichmentRatio = enrich_ratio,
        geneID = matched,
        Count = k
      ))
    }
  }

  results$p.adjust <- p.adjust(results$pvalue, method = p.adjust.method)
  results <- results[order(results$p.adjust), ]
  return(results)
}
