#' Visualize Pathway Networks Based on Shared Compounds
#'
#' Creates an interactive network visualization of biological pathways where edges represent shared metabolites.
#' 
#' @param pathway_data A list containing pathway information with two components:
#'   - `pathways`: Named vector of pathway IDs to pathway names
#'   - `pathway2compound`: List mapping pathway IDs to compound vectors
#' @param min_shared Minimum number of shared compounds required to draw an edge (default: 1)
#' @param layout Network layout algorithm. Options: "fr" (Fruchterman-Reingold), 
#'               "kk" (Kamada-Kawai), "circle", "dh" (Davidson-Harel) (default: "fr")
#' @param point_size Base scaling factor for node size (default: 3)
#' @param label_size Text size for node labels (default: 1)
#' @param edge_width Base scaling factor for edge width (default: 0.5)
#' @param vertex_color Fill color for nodes (default: "#8DA0CB" - lavender blue)
#' @param vertex_alpha Transparency for nodes (0-1, default: 0.85)
#' @param vertex_border Border color for nodes (default: "#4A4A8C" - dark blue)
#' @param edge_color Color for edges (default: "#666666" - medium gray)
#' @param edge_alpha Transparency for edges (0-1, default: 0.7)
#' @param main_title Main title for the plot (default: "KEGG Pathway Network")
#' @param legend_position Position for legend. Options: "bottomright", "bottom", 
#'                        "bottomleft", "left", "topleft", "top", "topright", 
#'                        "right", "center" (default: "bottomright")
#' @param show_legend Whether to display legend (default: TRUE)
#' @param margin Plot margins in order: bottom, left, top, right (default: c(0,0,0,0))
#' @param legend_inset Legend position adjustment (x,y). Positive values move inside, 
#'                     negative values move outside (default: c(-0.5, -0.03))
#' 
#' @return Invisibly returns the igraph network object
#' @export
#' 
#' @examples
#' # Example pathway data structure
#' pathway_data <- list(
#'   pathways = c("hsa00010" = "Glycolysis", 
#'                "hsa00020" = "Citrate cycle"),
#'   pathway2compound = list(
#'     "hsa00010" = c("C00031", "C00033"),
#'     "hsa00020" = c("C00033", "C00036")
#'   )
#' )
#' 
#' # Create network visualization
#' pathnetwork(
#'   pathway_data,
#'   min_shared = 1,
#'   vertex_color = "#5D9CEC",
#'   legend_inset = c(-0.4, -0.02)
#' )
pathnetwork <- function(pathway_data, min_shared = 1, 
                           layout = "fr", point_size = 3, 
                           label_size = 1, edge_width = 0.5,
                           vertex_color = "#8DA0CB", vertex_alpha = 0.85,
                           vertex_border = "#4A4A8C", 
                           edge_color = "#666666", edge_alpha = 0.7,
                           main_title = "KEGG Pathway Network", 
                           legend_position = "bottomright",
                           show_legend = TRUE,
                           margin = c(0, 0, 0, 0),
                           legend_inset = c(-0.5, -0.03)) {
  
  # 提取数据
  pathways <- pathway_data$pathways
  pathway2compound <- pathway_data$pathway2compound
  path_names <- names(pathway2compound)
  
  # 检查通路数量
  n_paths <- length(path_names)
  if (n_paths < 2) {
    stop("At least two pathways are required to create a network")
  }
  
  # 创建邻接矩阵（存储共享代谢物数量）
  adj_matrix <- matrix(0, nrow = n_paths, ncol = n_paths,
                       dimnames = list(path_names, path_names))
  
  # 填充共享代谢物计数
  message("Calculating shared compounds between pathways...")
  for (i in 1:(n_paths - 1)) {
    cmpds_i <- pathway2compound[[path_names[i]]]
    for (j in (i + 1):n_paths) {
      shared <- intersect(cmpds_i, pathway2compound[[path_names[j]]])
      adj_matrix[i, j] <- length(shared)
      adj_matrix[j, i] <- length(shared)  # 对称矩阵
    }
  }
  
  # 创建igraph网络
  net <- igraph::graph_from_adjacency_matrix(
    adj_matrix,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  # 设置顶点属性（添加边框美化）
  igraph::V(net)$name <- path_names
  igraph::V(net)$label <- sapply(path_names, function(pid) {
    name <- pathways[pid]
    if (is.null(name) || is.na(name)) pid else name
  })
  igraph::V(net)$label.color <- "#333333"  # 深灰色标签更清晰
  igraph::V(net)$label.family <- "sans"    # 使用无衬线字体
  
  # 计算代谢物数量（节点大小）
  compound_counts <- sapply(pathway2compound[path_names], length)
  igraph::V(net)$size <- sqrt(compound_counts) * point_size + 3
  
  # 设置边属性（美化边效果）
  igraph::E(net)$width <- igraph::E(net)$weight * edge_width
  igraph::E(net)$curved <- 0.2  # 固定曲线度
  net <- igraph::delete_edges(net, igraph::E(net)[igraph::E(net)$weight < min_shared])
  
  # 检查是否有边
  if (igraph::ecount(net) == 0) {
    stop("No shared compounds between pathways or below threshold")
  }
  
  # 设置布局函数
  layout_fun <- switch(layout,
                       "circle" = igraph::layout_in_circle,
                       "fr" = igraph::layout_with_fr,
                       "kk" = igraph::layout_with_kk,
                       "dh" = igraph::layout_with_dh,
                       igraph::layout_with_fr
  )
  
  # 计算布局
  coords <- layout_fun(net)
  
  # 准备颜色（带透明度）
  v_color <- grDevices::adjustcolor(vertex_color, alpha.f = vertex_alpha)
  v_border <- grDevices::adjustcolor(vertex_border, alpha.f = vertex_alpha)
  e_color <- grDevices::adjustcolor(edge_color, alpha.f = edge_alpha)
  
  # 设置绘图参数
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  graphics::par(mar = margin)  # 使用用户指定的边距
  
  # 绘制网络图（美化参数）
  igraph::plot.igraph(
    net,
    layout = coords,
    margin = rep(0, 4),  # 使用igraph默认边距
    vertex.label.cex = label_size,
    vertex.color = v_color,
    vertex.frame.color = v_border,  # 添加节点边框
    vertex.frame.width = 0.8,       # 边框宽度
    vertex.shape = "circle",        # 统一圆形节点
    edge.color = e_color,
    edge.curved = 0.2,
    main = main_title,
    main.color = "#2c3e50",         # 深蓝色标题
    main.font = 2,                  # 粗体标题
    sub.color = "#7f8c8d",          # 灰色副标题
    sub.cex = 0.9
  )
  
  # 添加美化后的图例
  if (show_legend) {
    # 重组图例文本（包含min_shared）
    legend_text <- c(
      "Node size: Pathway compounds",
      "Edge width: Shared compounds",
      paste0("Min shared: ", min_shared)
    )
    
    # 使用美化图例参数
    graphics::legend(
      legend_position,
      legend = legend_text,
      cex = 0.8,
      bty = "n",                # 无边框
      text.col = "#34495e",     # 深灰色文本
      inset = legend_inset,     # 使用调整后的偏移量
      xpd = TRUE
    )
  }
  
  # 返回网络对象
  return(invisible(net))
}
