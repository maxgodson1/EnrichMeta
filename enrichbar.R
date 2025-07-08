#' Title: Draw KEGG dotplot
#'
#' @description
#' Draw a barplot of KEGG enrichment results returned by `kegg_metabolite_enrichment()`.
#'
#' @param results A data.frame returned from kegg_metabolite_enrichment().
#' @return A ggplot2 object.
#' @export
enrichbar <- function(results, top = 25) {
  # 确保必要的包已加载
  library(ggplot2)
  library(stringr)  # 用于换行长标签
  
  # 筛选最显著的 top 条通路（按校正p值排序）
  top_results <- head(results[order(results$p.adjust), ], top)
  
  # 按校正p值升序排列（最显著的在上方）- 使用 p.adjust 排序
  top_results <- top_results[order(top_results$p.adjust, decreasing = FALSE), ]
  
  # 处理长标签 - 添加换行符
  top_results$Description <- str_wrap(top_results$Description, width = 40)
  
  # 创建因子水平
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
      y = NULL,  # 移除y轴标题
      title = "Overview of Enriched Metabolite Sets"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      # 调整y轴标签
      axis.text.y = element_text(
        size = 10, 
        color = "black", 
        face = "bold",
        margin = margin(r = 10)  # 右侧增加边距
      ),
      
      # 调整x轴标签
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
      
      # 标题设置
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
      
      # 图例设置
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      
      # 网格线设置
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      
      # 整体边距设置
      plot.margin = margin(1, 2.5, 1, 1, "cm"),  # 右侧增加更多空间
      
      # 绘图区域设置
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.3)),  # 右侧增加更多空间
      breaks = scales::pretty_breaks(n = 5)  # 自动生成合适的断点
    ) +
    coord_cartesian(clip = "off")  # 防止标签被截断
  
  return(p)
}
