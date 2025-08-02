#' Retrieve Pathway Information for Metabolites
#'
#' Identifies which KEGG pathways contain the specified metabolites.
#'
#' @param keggid Character vector of metabolite KEGG IDs (e.g., c("C00160", "C00041"))
#' @param pathway_data Pathway data object from \code{\link{getkeggdata}}
#'
#' @return A data frame with two columns:
#'   \itemize{
#'     \item \code{keggID}: KEGG ID of the metabolite
#'     \item \code{pathwayID}: KEGG pathway ID containing the metabolite
#'   }
#'   Rows are ordered by keggID and pathwayID. Metabolites not found in any pathways are excluded.
#'
#' @details
#' This function checks the provided pathway data to determine which pathways contain each input metabolite.
#' The result provides a comprehensive mapping between metabolites and their associated pathways.
#'
#' For large datasets (more than 50 metabolites), a progress bar is displayed to track processing status.
#'
#' @note
#' The pathway_data object must contain:
#' \itemize{
#'   \item \code{path_cpd_map}: List mapping pathway IDs to vectors of compound IDs
#' }
#'
#' @examples
#' \dontrun{
#' # Get pathway data
#' kegg_data <- getkeggdata("hsa", cache_dir = "./cache")
#'
#' # Get pathway information for metabolites
#' metabolite_pathways <- findcpdspaths(c("C00031", "C00221"), pathway_data = kegg_data)
#'
#' # View results
#' head(metabolite_pathways)
#' }
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export

findcpdspaths <- function(keggid, pathway_data) {
  # 验证必需参数
  if (missing(keggid)) stop("Missing required argument: keggid")
  if (missing(pathway_data)) stop("Missing required argument: pathway_data")

  # 检查pathway_data中必需的组成部分
  if (!"path_cpd_map" %in% names(pathway_data)) {
    stop("pathway_data must contain 'path_cpd_map' element")
  }

  # 提取通路-代谢物映射数据
  path_cpd_map <- pathway_data$path_cpd_map

  # 初始化空的结果数据框
  results <- data.frame(
    keggID = character(),
    pathwayID = character(),
    stringsAsFactors = FALSE
  )

  # 确定是否需要显示进度条（当代谢物数量大于50时）
  show_progress <- length(keggid) > 50
  if (show_progress) {
    message("Querying pathway information for metabolites...")
    pb <- utils::txtProgressBar(min = 0, max = length(keggid), style = 3)
  }

  # 处理每个代谢物ID
  for (i in seq_along(keggid)) {
    cpd <- keggid[i]

    # 检查每个通路是否包含当前代谢物
    for (path_id in names(path_cpd_map)) {
      if (cpd %in% path_cpd_map[[path_id]]) {
        # 将通路-代谢物对添加到结果中
        results <- rbind(results, data.frame(
          keggID = cpd,
          pathwayID = path_id,
          stringsAsFactors = FALSE
        ))
      }
    }

    # 如果需要显示进度条，则更新进度
    if (show_progress) utils::setTxtProgressBar(pb, i)
  }

  # 关闭进度条（如果已开启）
  if (show_progress) close(pb)

  # 处理并返回结果
  if (nrow(results) > 0) {
    # 按代谢物ID和通路ID排序结果
    results <- results[order(results$keggID, results$pathwayID), ]

    # 重置行名以获得更整洁的输出
    rownames(results) <- NULL

    return(results)
  } else {
    message("None of the provided metabolites were found in any pathways.")
    return(results)  # 返回具有正确结构的空数据框
  }
}
