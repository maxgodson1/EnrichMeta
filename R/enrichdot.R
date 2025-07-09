#' Title: Draw KEGG dotplot
#'
#' @description
#' Draw a dotplot of KEGG enrichment results returned by `keggenrich()`.
#'
#' @param results A data.frame returned from keggenrich().
#' @param top Number of top pathways to display (default=25)
#' @param title Plot title (default="Enriched Metabolic Pathways")
#' @param low_color Color for significant p-values (default="#FF0000")
#' @param high_color Color for non-significant p-values (default="#FFFF00")
#' @param point_size_range Point size range (default=c(3,8))
#' @param point_alpha Point transparency (default=0.8)
#' @param base_size Base font size (default=12)
#' @param wrap_width Description wrap width (default=40)
#' @param color_limits Color scale limits (default=c(0,0.15))
#' @param legend_position Legend position (default="right")
#' @param show_grid Show horizontal grid lines? (default=TRUE)
#' @param x_expansion Multiplicative expansion for x-axis (default=c(0,0.2))
#' @param padding_right Additional right padding for points (default=0.3)
#' @return A ggplot2 object.
#' @export
enrichdot <- function(results, 
                      top = 25,
                      title = "Enriched Metabolic Pathways",
                      low_color = "#FF0000",
                      high_color = "#FFFF00",
                      point_size_range = c(3, 8),
                      point_alpha = 0.8,
                      base_size = 12,
                      wrap_width = 40,
                      color_limits = c(0, 0.15),
                      legend_position = "right",
                      show_grid = TRUE,
                      x_expansion = c(0, 0.2),
                      padding_right = 0.3) {
  
  # 检查必要包
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请先安装ggplot2包：install.packages('ggplot2')")
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("请先安装stringr包：install.packages('stringr')")
  }
  
  # 按校正p值升序排序
  results <- results[order(results$p.adjust), ]
  
  # 取前top个通路
  if (nrow(results) > top) {
    results <- results[1:top, ]
  }
  
  # 处理长标签
  results$Description <- stringr::str_wrap(results$Description, width = wrap_width)
  
  # 确保通路按p.adjust升序排列
  results$Description <- factor(
    results$Description,
    levels = rev(unique(results$Description)))
    
  # 计算x轴的最大值（考虑点的大小）
  max_x_value <- max(-log10(results$pvalue))
  max_point_size <- max(point_size_range)
    
  # 计算需要扩展的空间（基于点的大小）
  x_max_limit <- max_x_value + (max_point_size * 0.005) + padding_right
    
  # 创建基础绘图
  p <- ggplot2::ggplot(results, 
                        ggplot2::aes(x = -log10(pvalue),
                                    y = Description,
                                    size = EnrichmentRatio,
                                    color = pvalue)) +
       ggplot2::geom_point(alpha = point_alpha) +
       ggplot2::scale_color_gradient(
         low = low_color, 
         high = high_color,
         limits = c(0, 0.2),  # 固定颜色范围为0-0.2
         breaks = seq(0, 0.2, by = 0.05),  # 设置固定刻度
         labels = function(x) sprintf("%.2f", x),  # 格式化标签显示两位小数
         name = "P-value"
       ) +
       ggplot2::scale_size_continuous(
         range = point_size_range,
         name = "Enrichment Ratio"
       ) +
       ggplot2::labs(
         x = expression(-log[10](P-value)),
         y = NULL,
         title = title
       ) +
       ggplot2::theme_minimal(base_size = base_size) +
       ggplot2::theme(
         axis.text.y = ggplot2::element_text(
           size = base_size - 2,
           face = "bold",
           margin = ggplot2::margin(r = 10)
         ),
         axis.text.x = ggplot2::element_text(size = base_size - 2),
         axis.title.x = ggplot2::element_text(
           size = base_size,
           face = "bold",
           margin = ggplot2::margin(t = 10)
         ),
         plot.title = ggplot2::element_text(
           size = base_size + 4,
           face = "bold", 
           hjust = 0.5,
           margin = ggplot2::margin(b = 15)
         ),
         legend.position = legend_position,
         legend.title = ggplot2::element_text(size = base_size - 2, face = "bold"),
         legend.text = ggplot2::element_text(size = base_size - 3),
         panel.grid.minor = ggplot2::element_blank(),
         plot.margin = ggplot2::unit(c(1, 1.5, 1, 1), "cm"),
         plot.background = ggplot2::element_rect(fill = "white", color = NA),
         panel.background = ggplot2::element_rect(fill = "white", color = NA)
       ) +
       ggplot2::scale_x_continuous(
         expand = ggplot2::expansion(mult = x_expansion),
         limits = c(0, x_max_limit))  # 设置x轴上限
        
        # 条件添加网格线
      if(show_grid) {
          p <- p + ggplot2::theme(
            panel.grid.major.y = ggplot2::element_line(color = "grey90")
          )
        } else {
          p <- p + ggplot2::theme(
            panel.grid.major.y = ggplot2::element_blank()
          )
        }
        
        return(p)
}
