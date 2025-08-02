#' Calculate Shared Metabolites Between Selected Pathways
#'
#' Analyzes metabolite sharing between user-specified pathways and returns a data frame of shared metabolites.
#'
#' @param pathwayid Character vector of pathway IDs to analyze (e.g., "hsa00010")
#' @param pathway_data Pathway data object from \code{\link{getkeggdata}} containing:
#'   \itemize{
#'     \item `path_cpd_map`: List of compounds per pathway
#'   }
#' @param min_shared Minimum shared metabolites for inclusion (default = 1)
#'
#' @return A data frame detailing shared metabolites between pathway pairs with columns:
#'   \itemize{
#'     \item from: Source pathway ID
#'     \item to: Target pathway ID
#'     \item shared_count: Number of shared metabolites
#'     \item keggID: Semicolon-separated list of shared metabolite IDs
#'   }
#'
#' @details
#' Computes pairwise metabolite sharing between user-selected pathways. Only pathway pairs with at least
#' `min_shared` metabolites are included. Pathways are compared based on their metabolite composition.
#'
#' @examples
#' \dontrun{
#' # After enrichment analysis
#' kegg_data <- getkeggdata("hsa", cache_dir = "./cache")
#' selected_paths <- c("hsa00010", "hsa00020", "hsa00030")
#' shared_df <- findsharedcpds(
#'   pathwayid = selected_paths,
#'   pathway_data = kegg_data,
#'   min_shared = 2
#' )
#'
#' # View shared metabolites
#' head(shared_df)
#' }
#'
#' @export

findsharedcpds <- function(pathwayid, pathway_data, min_shared = 1) {
  # 验证输入数据
  if (!"path_cpd_map" %in% names(pathway_data)) {
    stop("pathway_data must contain 'path_cpd_map' element")
  }
  if (!is.character(pathwayid)) {
    stop("pathwayid must be a character vector")
  }

  # 获取指定通路的代谢物映射
  path_cpd_map <- pathway_data$path_cpd_map

  # 查找有效通路 - 直接匹配ID
  valid_paths <- pathwayid[pathwayid %in% names(path_cpd_map)]

  # 检查无效通路ID
  if (length(valid_paths) != length(pathwayid)) {
    invalid <- setdiff(pathwayid, names(path_cpd_map))
    warning("Ignoring invalid pathway IDs: ", paste(invalid, collapse = ", "))
  }

  # 检查有效通路数量
  n_paths <- length(valid_paths)
  if (n_paths < 2) {
    message("At least two valid pathways are required. Returning empty data frame.")
    return(data.frame(
      from = character(),
      to = character(),
      shared_count = integer(),
      keggID = character(),
      stringsAsFactors = FALSE
    ))
  }

  # 预分配结果列表
  shared_list <- list()

  # 计算共享代谢物
  for (i in 1:(n_paths - 1)) {
    path1 <- valid_paths[i]
    cmpds1 <- path_cpd_map[[path1]]

    for (j in (i + 1):n_paths) {
      path2 <- valid_paths[j]
      cmpds2 <- path_cpd_map[[path2]]

      shared <- intersect(cmpds1, cmpds2)
      count <- length(shared)

      if (count >= min_shared) {
        shared_list[[length(shared_list) + 1]] <- data.frame(
          from = path1,
          to = path2,
          shared_count = count,
          keggID = paste(shared, collapse = ";"),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # 返回结果
  if (length(shared_list) > 0) {
    return(do.call(rbind, shared_list))
  } else {
    message("No pathway pairs found with >= ", min_shared, " shared metabolites.")
    return(data.frame(
      from = character(),
      to = character(),
      shared_count = integer(),
      keggID = character(),
      stringsAsFactors = FALSE
    ))
  }
}
