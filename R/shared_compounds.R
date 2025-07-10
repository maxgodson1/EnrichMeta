#' Calculate shared metabolites between pathways
#' 
#' @param pathway_data Data structure containing pathway information
#' @param min_shared Minimum threshold for shared metabolite count (default: 1)
#' 
#' @return A list containing the following elements:
#'   - graph: igraph network object
#'   - shared_df: data.frame containing detailed shared metabolite information
#'   - pathway_names: pathway name mapping vector
#'   - compound_counts: vector of metabolite counts for each pathway
#' @export
shared_compounds <- function(pathway_data, min_shared = 1) {
  # 提取数据
  pathways <- pathway_data$pathways
  pathway2compound <- pathway_data$pathway2compound
  path_names <- names(pathway2compound)
  
  # 检查通路数量
  n_paths <- length(path_names)
  if (n_paths < 2) {
    stop("至少需要两个通路才能进行分析")
  }
  
  # 创建邻接矩阵
  adj_matrix <- matrix(0, nrow = n_paths, ncol = n_paths,
                       dimnames = list(path_names, path_names))
  
  # 存储共享代谢物详细信息
  shared_list <- list()
  
  # 填充共享代谢物计数和详细信息
  for (i in 1:(n_paths - 1)) {
    cmpds_i <- pathway2compound[[path_names[i]]]
    for (j in (i + 1):n_paths) {
      shared <- intersect(cmpds_i, pathway2compound[[path_names[j]]])
      count <- length(shared)
      
      if (count >= min_shared) {
        adj_matrix[i, j] <- count
        adj_matrix[j, i] <- count
        
        # 存储共享代谢物详细信息
        edge_id <- paste(sort(c(path_names[i], path_names[j])), collapse = "|")
        shared_list[[edge_id]] <- data.frame(
          from = path_names[i],
          to = path_names[j],
          shared_count = count,
          compounds = paste(shared, collapse = ";")
        )
      }
    }
  }
  
  # 创建igraph网络
  net <- igraph::graph_from_adjacency_matrix(
    adj_matrix,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  # 检查是否有边
  if (igraph::ecount(net) == 0) {
    stop("通路间没有共享代谢物或低于阈值")
  }
  
  # 设置顶点属性
  igraph::V(net)$name <- path_names
  igraph::V(net)$label <- pathways[path_names]
  
  # 计算每个通路的代谢物数量
  compound_counts <- sapply(pathway2compound[path_names], length)
  names(compound_counts) <- path_names
  
  # 返回结果
  shared_df <- if (length(shared_list) > 0) do.call(rbind, shared_list) else NULL
  rownames(shared_df) <- NULL
  
  list(
    graph = net,
    shared_df = shared_df,
    pathway_names = pathways,
    compound_counts = compound_counts
  )
}

