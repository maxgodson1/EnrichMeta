#' Title: Draw KEGG dotplot
#'
#' @description
#' Draw a barplot of KEGG enrichment results returned by `kegg_metabolite_enrichment()`.
#'
#' @param results A data.frame returned from kegg_metabolite_enrichment().
#' @return A ggplot2 object.
#' @export
keggbar <- function(results, top = 25) {
  # 确保必要的包已加载
  library(ggplot2)
  
  # 筛选最显著的 top 条通路（按校正p值排序）
  top_results <- head(results[order(results$p.adjust), ], top)
  
  # 按校正p值升序排列（最显著的在上方）- 使用 p.adjust 排序
  top_results <- top_results[order(top_results$p.adjust, decreasing = FALSE), ]
  top_results$Description <- factor(top_results$Description, 
                                    levels = rev(top_results$Description))
  
  # 创建绘图
  p <- ggplot(top_results, aes(
    x = EnrichmentRatio,
    y = Description,
    fill = pvalue  # 使用原始p值填充颜色
  )) +
    geom_bar(stat = "identity", width = 0.8) +
    geom_text(
      aes(label = sprintf("%.1f", EnrichmentRatio)),
      hjust = -0.2,
      size = 3.5,
      color = "black",
      fontface = "bold"
    ) +
    scale_fill_gradient(
      name = "P-value",
      low = "#FF0000",  # 深红色表示显著（p值小）
      high = "#FFFF00",  # 浅黄色表示不显著（p值大）
      limits = c(0, 0.2),  # 设置颜色范围为0-0.2
      breaks = seq(0, 0.2, by = 0.05),  # 设置图例断点
      labels = function(x) sprintf("%.2f", x)
    ) +
    labs(
      x = "Enrichment Ratio", 
      y = "Metabolite Sets",
      title = "Overview of Enriched Metabolite Sets (Top 25)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 11, color = "black", face = "bold"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 2, 1, 1, "cm")
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.2)),
      breaks = scales::pretty_breaks(n = 5)  # 自动生成合适的断点
    ) +
    coord_cartesian(clip = "off")
  
  return(p)
}
