#' Title: Draw KEGG dotplot
#'
#' @description
#' Draw a dotplot of KEGG enrichment results returned by `kegg_metabolite_enrichment()`.
#'
#' @param results A data.frame returned from kegg_metabolite_enrichment().
#' @return A ggplot2 object.
#' @export
keggdot <- function(results) {
  ggplot(results, aes(x = reorder(Description, p.adjust),
                      y = GeneRatio,
                      size = Count,
                      color = p.adjust)) +
    geom_point(alpha = 0.8) +
    coord_flip() +
    scale_color_gradient(low = "blue", high = "red") +
    labs(x = "Pathway", y = "GeneRatio",
         color = "p.adjust", size = "Count") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
}
