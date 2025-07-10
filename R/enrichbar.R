#' Draw a Bar Plot for KEGG Enrichment Results
#'
#' Creates a horizontal bar plot to visualize KEGG enrichment results, with bars colored by p-value significance.
#'
#' @param results A data frame containing KEGG enrichment results as returned by `keggenrich()`.
#'                Must contain columns: `Description`, `EnrichmentRatio`, `pvalue`, `p.adjust`.
#' @param top Number of top enriched pathways to display (default = 25)
#' @param title Plot title (default = "Overview of Enriched Metabolite Sets")
#' @param low_color Color for the most significant p-values (default = "#FF0000" - red)
#' @param high_color Color for the least significant p-values (default = "#FFFF00" - yellow)
#' @param bar_width Width of the bars (default = 0.8)
#' @param base_size Base font size for the plot (default = 12)
#' @param wrap_width Number of characters per line for wrapping pathway descriptions (default = 40)
#' @param show_ratio Logical indicating whether to display enrichment ratio values on the bars (default = TRUE)
#' @param color_limits Range of p-values for the color scale (default = c(0, 0.2))
#' @param legend_position Position of the legend ("none", "left", "right", "bottom", "top", or two-element numeric vector; default = "right")
#'
#' @return A ggplot2 object representing the enrichment bar plot.
#'
#' @details
#' This function creates a horizontal bar plot where:
#' \itemize{
#'   \item Each bar represents a KEGG pathway
#'   \item Bar length corresponds to the enrichment ratio
#'   \item Bar color indicates the p-value (red = significant, yellow = non-significant)
#'   \item Pathways are ordered by adjusted p-value (most significant at the top)
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming 'enrich_results' is output from keggenrich()
#'
#' # Basic plot
#' enrichbar(enrich_results)
#'
#' # Customized plot
#' enrichbar(
#'   results = enrich_results,
#'   top = 15,
#'   title = "Top Enriched Pathways",
#'   low_color = "darkred",
#'   high_color = "lightyellow",
#'   bar_width = 0.7,
#'   base_size = 14,
#'   wrap_width = 30,
#'   show_ratio = FALSE,
#'   color_limits = c(0, 0.1),
#'   legend_position = "bottom"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_minimal scale_fill_gradient
#' @importFrom ggplot2 theme element_text margin element_rect scale_x_continuous
#' @importFrom ggplot2 coord_cartesian geom_text
#' @importFrom stringr str_wrap
#' @importFrom scales pretty_breaks
#' @export
enrichbar <- function(results,
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

  # 筛选并排序数据
  top_results <- head(results[order(results$p.adjust), ], top)
  top_results <- top_results[order(top_results$p.adjust, decreasing = FALSE), ]

  # 处理长标签
  top_results$Description <- stringr::str_wrap(top_results$Description, width = wrap_width)
  top_results$Description <- factor(top_results$Description,
                                    levels = rev(top_results$Description))

  # 创建基础绘图
  p <- ggplot2::ggplot(top_results, ggplot2::aes(
    x = EnrichmentRatio,
    y = Description,
    fill = pvalue
  )) +
    ggplot2::geom_bar(stat = "identity", width = bar_width) +
    ggplot2::labs(
      x = "Enrichment Ratio",
      y = NULL,
      title = title
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::scale_fill_gradient(
      name = "P-value",
      low = low_color,
      high = high_color,
      limits = color_limits,
      breaks = seq(min(color_limits), max(color_limits), length.out = 5),
      labels = function(x) sprintf("%.2f", x)
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        size = base_size - 2,
        color = "black",
        face = "bold",
        margin = ggplot2::margin(r = 10)
      ),
      axis.text.x = ggplot2::element_text(size = base_size - 2, color = "black"),
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
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(1, 2.5, 1, 1, "cm"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.3)),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  # 条件添加富集比标签
  if(show_ratio) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f", EnrichmentRatio)),
      hjust = -0.2,
      size = base_size - 8.5,
      color = "black",
      fontface = "bold"
    )
  }

  return(p)
}


