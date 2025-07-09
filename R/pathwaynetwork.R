#' Plot KEGG Pathway Network
#'
#' Visualize metabolic pathways as a network based on shared compounds.
#'
#' @param pathway_data List containing pathway data from \code{getkpd} or RDS file.
#' @param min_shared Minimum number of shared compounds for an edge (default=1).
#' @param layout Network layout algorithm: "circle", "fr", "kk", or "dh" (default="fr").
#' @param point_size Base point size multiplier (default=3).
#' @param label_size Vertex label size (default=1).
#' @param edge_width Edge width multiplier (default=0.5).
#' @param vertex_color Color for vertices (default="#8DA0CB" - light purple).
#' @param vertex_alpha Vertex transparency (0-1, default=0.7).
#' @param edge_color Color for edges (default="#666666" - dark gray).
#' @param edge_alpha Edge transparency (0-1, default=0.5).
#' @param main_title Main plot title (default="KEGG Pathway Network").
#' @param sub_title Subtitle (default=NULL, auto-generated).
#' @param legend_position Legend position ("bottomright", "topleft", etc., default="bottomright").
#' @param show_legend Show legend? (default=TRUE).
#' @param margin Plot margin (default=c(0,0,0,0)).
#' @param legend_inset Legend inset from plot edges (default=c(0.02, 0.02)).
#'
#' @return Invisible igraph network object.
#'
#' @details
#' This function creates a network visualization where:
#' - Nodes represent KEGG pathways
#' - Node size is proportional to number of compounds in the pathway
#' - Edges represent shared compounds between pathways
#' - Edge width is proportional to number of shared compounds
#'
#' @examples
#' \dontrun{
#' # Load pathway data
#' pathway_data <- readRDS("kegg_pathways_hsa.rds")
#'
#' # Basic plot
#' plot_pathway_network(pathway_data)
#'
#' # Customized plot with more margin
#' plot_pathway_network(
#'   pathway_data,
#'   min_shared = 3,
#'   margin = c(0.05, 0.05, 0.1, 0.05), # bottom, left, top, right
#'   legend_inset = c(0.05, 0.05)
#' )
#' }
#'
#' @importFrom igraph graph_from_adjacency_matrix V E delete_edges ecount
#'             layout_in_circle layout_with_fr layout_with_kk layout_with_dh
#' @importFrom grDevices adjustcolor
#' @export
pathwaynetwork <- function(pathway_data, min_shared = 1, 
                                 layout = "fr", point_size = 3, 
                                 label_size = 1, edge_width = 0.5,
                                 vertex_color = "#8DA0CB", vertex_alpha = 0.7,
                                 edge_color = "#666666", edge_alpha = 0.5,
                                 main_title = "KEGG Pathway Network", 
                                 sub_title = NULL,
                                 legend_position = "bottomright",
                                 show_legend = TRUE,
                                 margin = c(0, 0, 0, 0),  # bottom, left, top, right
                                 legend_inset = c(0.02, 0.02)) {
  
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
  
  # 设置顶点属性
  igraph::V(net)$name <- path_names
  igraph::V(net)$label <- sapply(path_names, function(pid) {
    name <- pathways[pid]
    if (is.null(name) || is.na(name)) pid else name
  })
  
  # 计算代谢物数量（节点大小）
  compound_counts <- sapply(pathway2compound[path_names], length)
  igraph::V(net)$size <- sqrt(compound_counts) * point_size + 3
  
  # 设置边属性（过滤弱连接）
  igraph::E(net)$width <- igraph::E(net)$weight * edge_width
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
  
  # 自动生成副标题
  if (is.null(sub_title)) {
    sub_title <- paste("Min shared:", min_shared)  # 更短的文本
  }
  
  # 计算布局
  coords <- layout_fun(net)
  
  # 准备颜色（带透明度）
  v_color <- grDevices::adjustcolor(vertex_color, alpha.f = vertex_alpha)
  e_color <- grDevices::adjustcolor(edge_color, alpha.f = edge_alpha)
  
  # 设置绘图参数
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  graphics::par(mar = c(5, 4, 4, 2) + 0.1)  # 默认边距
  
  # 绘制网络图（添加边距）
  igraph::plot.igraph(
    net,
    layout = coords,
    margin = margin,  # 控制图形边距
    vertex.label.cex = label_size,
    vertex.label.color = "black",
    vertex.color = v_color,
    vertex.frame.color = NA,
    edge.color = e_color,
    edge.curved = 0.2,
    main = main_title,
    sub = sub_title
  )
  
  # 添加图例（避免遮挡）
  if (show_legend) {
    # 简化的图例文本
    legend_text <- c(
      "Node: Pathway compounds",
      "Edge: Shared compounds"
    )
    
    # 自动调整图例位置
    graphics::legend(
      legend_position,
      legend = legend_text,
      cex = 0.8,
      bty = "n",
      inset = legend_inset,  # 控制图例内边距
      xpd = TRUE  # 允许在绘图区域外绘制
    )
  }
  
  # 返回网络对象
  return(invisible(net))
}
