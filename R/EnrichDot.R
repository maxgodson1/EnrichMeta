#' Title: Draw KEGG dotplot
#'
#' @description
#' Draw a dotplot of KEGG enrichment results returned by `KeggEnrich()`.
#'
#' @param results A data.frame returned from KeggEnrich().
#' @return A ggplot2 object.
#' @export
EnrichDot <- function(results, top = 25) {
  # 确保所需包已加载
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请先安装ggplot2包：install.packages('ggplot2')")
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("请先安装stringr包：install.packages('stringr')")
  }
  
  # 按校正p值升序排序（最小p值在最前面）
  results <- results[order(results$p.adjust), ]
  
  # 取前top个通路
  if (nrow(results) > top) {
    results <- results[1:top, ]
  }
  
  # 处理长标签 - 添加换行符
  results$Description <- stringr::str_wrap(results$Description, width = 40)
  
  # 确保通路按p.adjust升序排列（最显著的在顶部）
  results$Description <- factor(
    results$Description,
    levels = rev(unique(results$Description)))
    
    # 创建ggplot对象
    ggplot2::ggplot(results, 
                    ggplot2::aes(x = -log10(pvalue),
                                 y = Description,
                                 size = EnrichmentRatio,
                                 color = pvalue)) +
      ggplot2::geom_point(alpha = 0.8) +
      ggplot2::scale_color_gradient(
        low = "red", 
        high = "yellow",
        breaks = c(0.00, 0.05, 0.10, 0.15),
        limits = c(0, 0.15),
        name = "P-value"  # 添加图例名称
      ) +
      ggplot2::scale_size_continuous(
        range = c(3, 8),
        breaks = scales::pretty_breaks(n = 4),  # 自动生成断点
        name = "Enrichment Ratio"
      ) +
      ggplot2::labs(
        x = expression(-log[10](P-value)),  # 使用数学表达式
        y = NULL,  # 移除y轴标题
        title = "Enriched Metabolic Pathways"  # 优化标题
      ) +
      ggplot2::theme_minimal(base_size = 12) +  # 设置基础字体大小
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(
          size = 10,  # 增大字体大小
          face = "bold",
          margin = ggplot2::margin(r = 10)  # 右侧增加边距
        ),
        axis.text.x = ggplot2::element_text(size = 10),
        axis.title.x = ggplot2::element_text(
          size = 12, 
          face = "bold",
          margin = ggplot2::margin(t = 10)
        ),
        plot.title = ggplot2::element_text(
          size = 16, 
          face = "bold", 
          hjust = 0.5,  # 标题居中
          margin = ggplot2::margin(b = 15)  # 底部增加边距
        ),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 10, face = "bold"),
        legend.text = ggplot2::element_text(size = 9),
        panel.grid.major.y = ggplot2::element_line(color = "grey90"),  # 仅保留水平网格线
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(1, 1.5, 1, 1), "cm"),  # 增加右侧边距
        plot.background = ggplot2::element_rect(fill = "white", color = NA),  # 设置白色背景
        panel.background = ggplot2::element_rect(fill = "white", color = NA)
      ) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.1)))  # 右侧增加空间
}
