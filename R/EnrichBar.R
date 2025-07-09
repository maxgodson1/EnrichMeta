#' Title: Draw KEGG dotplot
#'
#' @description
#' Draw a barplot of KEGG enrichment results returned by `KeggEnrich()`.
#'
#' @param results A data.frame returned from KeggEnrich().
#' @param top Number of top pathways to display (default=25)
#' @param title Plot title (default="Overview of Enriched Metabolite Sets")
#' @param low_color Color for significant p-values (default="#FF0000")
#' @param high_color Color for non-significant p-values (default="#FFFF00")
#' @param bar_width Bar width (default=0.8)
#' @param base_size Base font size (default=12)
#' @param wrap_width Description wrap width (default=40)
#' @param show_ratio Show enrichment ratio values? (default=TRUE)
#' @param color_limits Color scale limits (default=c(0, 0.2))
#' @param legend_position Legend position (default="right")
#' @return A ggplot2 object.
#' @export
EnrichBar <- function(results, 
                      top = 25,
                      title = "Overview of Enriched Metabolite Sets",
                      low_color = "#FF0000",
                      high_color = "#FFFF00",
                      bar_width = 0.8,
                      base_size = 12,
                      wrap_width = 40,
                      show_ratio = TRUE,
                      color_limits = c(0, 0.2),
                      legend_position = "right") {
  
  library(ggplot2)
  library(stringr)
  
  # 筛选并排序数据
  top_results <- head(results[order(results$p.adjust), ], top)
  top_results <- top_results[order(top_results$p.adjust, decreasing = FALSE), ]
  
  # 处理长标签
  top_results$Description <- str_wrap(top_results$Description, width = wrap_width)
  top_results$Description <- factor(top_results$Description, 
                                    levels = rev(top_results$Description))
  
  # 创建基础绘图
  p <- ggplot(top_results, aes(
    x = EnrichmentRatio,
    y = Description,
    fill = pvalue
  )) +
    geom_bar(stat = "identity", width = bar_width) +
    labs(
      x = "Enrichment Ratio", 
      y = NULL,
      title = title
    ) +
    theme_minimal(base_size = base_size) +
    scale_fill_gradient(
      name = "P-value",
      low = low_color,
      high = high_color,
      limits = color_limits,
      breaks = seq(min(color_limits), max(color_limits), length.out = 5),
      labels = function(x) sprintf("%.2f", x)
    ) +
    theme(
      axis.text.y = element_text(
        size = base_size - 2, 
        color = "black", 
        face = "bold",
        margin = margin(r = 10)
      ),
      axis.text.x = element_text(size = base_size - 2, color = "black"),
      axis.title.x = element_text(
        size = base_size, 
        face = "bold", 
        margin = margin(t = 10)
      ),
      plot.title = element_text(
        size = base_size + 4, 
        face = "bold", 
        hjust = 0.5, 
        margin = margin(b = 15)
      ),
      legend.position = legend_position,
      legend.title = element_text(size = base_size - 2, face = "bold"),
      legend.text = element_text(size = base_size - 3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 2.5, 1, 1, "cm"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.3)),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    coord_cartesian(clip = "off")
  
  # 条件添加富集比标签
  if(show_ratio) {
    p <- p + geom_text(
      aes(label = sprintf("%.1f", EnrichmentRatio)),
      hjust = -0.2,
      size = base_size - 8.5,
      color = "black",
      fontface = "bold"
    )
  }
  
  return(p)
}
