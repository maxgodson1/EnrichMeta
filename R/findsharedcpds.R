#' Calculate Shared Metabolites Between Pathways
#'
#' Analyzes metabolite sharing between pathways without network construction.
#'
#' @param pathway_data Pathway data from `getkeggdata()` containing:
#'   \itemize{
#'     \item `pathways`: Named list of pathway names (ID:name)
#'     \item `pathscpds`: List of compounds per pathway
#'   }
#' @param min_shared Minimum shared metabolites for inclusion (default = 1)
#'
#' @return A list containing:
#'   \item{shared_df}{Data frame detailing shared metabolites between pathway pairs}
#'   \item{pathway_names}{Named list mapping pathway IDs to names}
#'   \item{compound_counts}{Named vector of metabolite counts per pathway}
#'
#' @details
#' Computes pairwise metabolite sharing between pathways. Returns a data frame with:
#' \itemize{
#'   \item from: Source pathway ID
#'   \item to: Target pathway ID
#'   \item shared_count: Number of shared metabolites
#'   \item keggID: Semicolon-separated list of shared metabolite IDs
#' }
#'
#' @examples
#' \dontrun{
#' # After enrichment analysis
#' filtered_data <- list(
#'   pathways = kegg_data$pathways[selected_pathways],
#'   pathscpds = kegg_data$pathscpds[selected_pathways]
#' )
#' shared_result <- findsharedcpds(filtered_data, min_shared = 2)
#'
#' # View shared metabolites
#' head(shared_result$shared_df)
#' }
#'
#' @export

findsharedcpds <- function(pathway_data, min_shared = 1) {
  # 提取数据
  pathways <- pathway_data$pathways
  pathscpds <- pathway_data$pathscpds
  path_names <- names(pathscpds)

  # 检查通路数量
  n_paths <- length(path_names)
  if (n_paths < 2) {
    stop("At least two pathways are required for analysis")
  }

  # 创建邻接矩阵
  adj_matrix <- matrix(0, nrow = n_paths, ncol = n_paths,
                       dimnames = list(path_names, path_names))

  # 存储共享代谢物详细信息
  shared_list <- list()

  # 填充共享代谢物计数和详细信息
  for (i in 1:(n_paths - 1)) {
    cmpds_i <- pathscpds[[path_names[i]]]
    for (j in (i + 1):n_paths) {
      shared <- intersect(cmpds_i, pathscpds[[path_names[j]]])
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
          keggID = paste(shared, collapse = ";")
        )
      }
    }
  }


  # 计算每个通路的代谢物数量
  compound_counts <- sapply(pathscpds[path_names], length)
  names(compound_counts) <- path_names

  # 返回结果
  shared_df <- if (length(shared_list) > 0) do.call(rbind, shared_list) else NULL
  rownames(shared_df) <- NULL

  list(
    shared_df = shared_df,
    pathway_names = pathways,
    compound_counts = compound_counts
  )
}
