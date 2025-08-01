#' Retrieve Pathway Information for Metabolites
#'
#' Identifies which KEGG pathways contain the specified metabolites.
#'
#' @param keggid Vector of metabolite KEGG IDs (e.g., c("C00160", "C00041"))
#' @param pathway_data Pathway data obtained from `getkeggdata()`
#'
#' @return A data frame with three columns:
#'   \itemize{
#'     \item \code{keggID}: KEGG ID of the metabolite
#'     \item \code{pathwayID}: KEGG pathway ID containing the metabolite
#'     \item \code{Description}: Full name of the pathway
#'   }
#'   Rows are ordered by keggid and pathwayID. Metabolites not found in any pathways are excluded.
#'
#' @details
#' This function checks the provided pathway data to determine which pathways contain each input metabolite.
#' The result provides a comprehensive mapping between metabolites and their associated pathways.
#'
#' @note
#' The pathway_data object must contain:
#' \itemize{
#'   \item \code{pathways}: Named list mapping pathway IDs to pathway names
#'   \item \code{pathscpds}: List mapping pathway IDs to vectors of compound IDs
#' }
#'
#' @examples
#' \dontrun{
#' # Get pathway data
#' kegg_data <- getkeggdata("hsa", cache_dir = "./cache")
#'
#' # Get pathway information for metabolites
#' findcpdspaths(c("C00031", "C00221"), pathway_data = kegg_data)
#' }
#'
#' @export


findcpdspaths <- function(keggid, pathway_data) {
  # 验证输入
  if (missing(keggid)) stop("Missing required argument: keggid")
  if (missing(pathway_data)) stop("Missing required argument: pathway_data")
  if (!all(c("pathways", "pathscpds") %in% names(pathway_data))) {
    stop("Invalid pathway_data structure. Must contain 'pathways' and 'pathscpds' elements.")
  }

  # 准备通路名称映射
  pathway_names <- pathway_data$pathways
  pathscpds <- pathway_data$pathscpds

  # 创建空结果列表
  results <- list()

  # 进度条（当处理大量代谢物时）
  if (length(keggid) > 50) {
    message("Querying pathway information for metabolites...")
    pb <- utils::txtProgressBar(min = 0, max = length(keggid), style = 3)
  }

  # 查找每个代谢物的通路
  for (i in seq_along(keggid)) {
    compound_id <- keggid[i]
    found_pathways <- character(0)

    # 检查哪些通路包含当前代谢物
    for (path_id in names(pathscpds)) {
      if (compound_id %in% pathscpds[[path_id]]) {
        found_pathways <- c(found_pathways, path_id)
      }
    }

    # 如果有找到通路则添加到结果
    if (length(found_pathways) > 0) {
      for (path_id in found_pathways) {
        results[[length(results) + 1]] <- data.frame(
          keggID = compound_id,
          pathwayID = path_id,
          Description = pathway_names[[path_id]],  # 直接访问列表元素
          stringsAsFactors = FALSE
        )
      }
    }

    # 更新进度条
    if (length(keggid) > 50) utils::setTxtProgressBar(pb, i)
  }

  # 关闭进度条
  if (length(keggid) > 50) close(pb)

  # 合并结果并排序
  if (length(results) > 0) {
    df <- do.call(rbind, results)
    df <- df[order(df$keggID, df$pathwayID), ]
    rownames(df) <- NULL
    return(df)
  } else {
    message("None of the provided metabolites were found in any pathways.")
    return(data.frame(keggID = character(0),
                      pathwayID = character(0),
                      Description = character(0)))
  }
}

