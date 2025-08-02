#' Draw a Dot Plot for KEGG Enrichment Results
#'
#' Creates a dot plot to visualize KEGG enrichment results, with point size representing enrichment ratio and color representing p-value.
#'
#' @param results A data frame containing KEGG enrichment results as returned by `enrichkegg()`
#' @param top Number of top pathways to display (default = 25)
#' @param title Plot title (default = "Enriched Metabolic Pathways")
#' @param low_color Color for significant p-values (default = "#FF0000" - red)
#' @param high_color Color for non-significant p-values (default = "#FFFF00" - yellow)
#' @param point_size_range Range of point sizes (default = c(3, 8))
#' @param point_alpha Point transparency (default = 0.8)
#' @param base_size Base font size (default = 12)
#' @param wrap_width Number of characters per line for wrapping descriptions (default = 40)
#' @param color_limits Color scale limits (default = c(0, 0.15))
#' @param legend_position Legend position ("none", "left", "right", "bottom", "top"; default = "right")
#' @param show_grid Show horizontal grid lines? (default = TRUE)
#' @param x_expansion Multiplicative expansion for x-axis (default = c(0, 0.15))
#' @param padding_right Additional right padding for points (default = 0.3)
#'
#' @return A ggplot2 object
#'
#' @details
#' This function creates a dot plot where:
#' \itemize{
#'   \item Point position (x-axis) represents -log10(p-value)
#'   \item Point size represents fold enrichment
#'   \item Point color indicates p-value significance
#'   \item Pathways are ordered by adjusted p-value
#' }
#'
#' @examples
#' \dontrun{
#' # Basic dot plot
#' plotenrichdot(enrich_results)
#'
#' # Customized dot plot
#' plotenrichdot(
#'   results = enrich_results,
#'   top = 15,
#'   title = "Top Metabolic Pathways",
#'   low_color = "darkred",
#'   point_size_range = c(4, 10),
#'   base_size = 14
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal scale_color_gradient
#' @importFrom ggplot2 scale_size_continuous theme element_text margin element_rect
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom stringr str_wrap
#' @export

plotenrichdot <- function(results,
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
                          x_expansion = c(0, 0.15),
                          padding_right = 0.3) {

  # 检查必要包
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install the ggplot2 package first: install.packages('ggplot2')")
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Please install the 'stringr' package first: install.packages('stringr')")
  }

  # 按校正p值升序排序
  results <- results[order(results$p.adjust), ]

  # 取前top个通路
  if (nrow(results) > top) {
    results <- results[1:top, ]
  }

  # 处理pvalue=0的情况（防止无限大值）
  results$pvalue_plot <- results$pvalue
  results$pvalue_plot[results$pvalue_plot == 0] <- .Machine$double.xmin

  # 处理长标签
  results$Description <- stringr::str_wrap(results$Description, width = wrap_width)

  # 确保通路按p.adjust升序排列
  results$Description <- factor(
    results$Description,
    levels = rev(unique(results$Description)))

  # 计算x轴的最大值
  max_x_value <- max(-log10(results$pvalue_plot))
  max_point_size <- max(point_size_range)
  x_max_limit <- max_x_value + (max_point_size * 0.005) + padding_right

  # 创建基础绘图
  p <- ggplot2::ggplot(results,
                       ggplot2::aes(x = -log10(pvalue_plot),
                                    y = Description,
                                    size = EnrichmentRatio,
                                    color = pvalue)) +
    ggplot2::geom_point(alpha = point_alpha) +
    ggplot2::scale_color_gradient(
      low = low_color,
      high = high_color,
      limits = color_limits,
      breaks = seq(min(color_limits), max(color_limits), length.out = 5),
      labels = function(x) sprintf("%.3f", x),
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
      limits = c(0, x_max_limit))

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
