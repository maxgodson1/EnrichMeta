#' Title: Draw KEGG dotplot
#'
#' @description
#' Draw a dotplot of KEGG enrichment results returned by `kegg_metabolite_enrichment()`.
#'
#' @param results A data.frame returned from kegg_metabolite_enrichment().
#' @return A ggplot2 object.
#' @export
enrichdot <- function(results, top = 25) {
  # 确保所需包已加载
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请先安装ggplot2包：install.packages('ggplot2')")
  }
  
  # 按校正p值升序排序（最小p值在最前面）
  results <- results[order(results$p.adjust), ]
  
  # 取前top个通路
  if (nrow(results) > top) {
    results <- results[1:top, ]
  }
  
  # 关键修正：确保通路按p.adjust升序排列（最显著的在顶部）
  results$Description <- factor(
    results$Description,
    levels = rev(unique(results$Description))  # 反转顺序
  )
  
  # 创建ggplot对象
  ggplot2::ggplot(results, 
                  ggplot2::aes(x = -log10(pvalue),
                               y = Description,  # 直接使用因子化的Description
                               size = EnrichmentRatio,
                               color = pvalue)) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_color_gradient(
      low = "red", 
      high = "blue",
      breaks = c(0.00, 0.05, 0.10, 0.15),
      limits = c(0, 0.15)
    ) +
    ggplot2::scale_size_continuous(
      range = c(3, 8),
      breaks = c(10, 15, 20, 25),
      name = "Enrichment Ratio"
    ) +
    ggplot2::labs(
      x = "-log10(p-value)", 
      y = "Metabolic Pathway",
      color = "P-value"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8, face = "bold"),
      axis.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.position = "right",
      panel.grid.major = ggplot2::element_line(color = "grey90"),
      plot.margin = ggplot2::unit(c(1,1,1,1), "cm")
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0.1))
}
