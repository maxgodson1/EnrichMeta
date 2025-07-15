#' Retrieve Pathway Information for Metabolites
#'
#' Identifies which KEGG pathways contain the specified metabolites.
#'
#' @param KEGGid Vector of metabolite KEGG IDs (e.g., c("C00160", "C00041"))
#' @param pathway_data Pathway data obtained from `getkd()`
#'
#' @return A data frame with three columns:
#'   \itemize{
#'     \item \code{keggID}: KEGG ID of the metabolite
#'     \item \code{pathwayID}: KEGG pathway ID containing the metabolite
#'     \item \code{pathwayName}: Full name of the pathway
#'   }
#'   Rows are ordered by keggID and pathwayID. Metabolites not found in any pathways are excluded.
#'
#' @details
#' This function checks the provided pathway data to determine which pathways contain each input metabolite.
#' The result provides a comprehensive mapping between metabolites and their associated pathways.
#'
#' @examples
#' \dontrun{
#' # Get pathway data
#' kegg_data <- getkd("hsa", cache_dir = "./cache")
#'
#' # Get pathway information for metabolites
#' metapathinfor(c("C00031", "C00221"), pathway_data = kegg_data)
#' }
#'
#' @export

metapathinfor <- function(KEGGid, pathway_data) {
  # 验证输入
  if (missing(KEGGid)) stop("Missing required argument: KEGGid")
  if (missing(pathway_data)) stop("Missing required argument: pathway_data")
  if (!all(c("pathways", "pathway2compound") %in% names(pathway_data))) {
    stop("Invalid pathway_data structure. Must contain 'pathways' and 'pathway2compound' elements.")
  }

  # 准备通路名称映射（移除"path:"前缀作为ID）
  pathway_names <- setNames(
    unname(pathway_data$pathways),
    sub("^path:", "", names(pathway_data$pathways))
  )

  # 准备通路到代谢物的映射
  pathway2compound <- pathway_data$pathway2compound

  # 创建空结果列表
  results <- list()

  # 进度条（当处理大量代谢物时）
  if (length(KEGGid) > 50) {
    message("Querying pathway information for metabolites...")
    pb <- utils::txtProgressBar(min = 0, max = length(KEGGid), style = 3)
  }

  # 查找每个代谢物的通路
  for (i in seq_along(KEGGid)) {
    compound_id <- KEGGid[i]
    found_pathways <- character(0)

    # 检查哪些通路包含当前代谢物
    for (path_id in names(pathway2compound)) {
      if (compound_id %in% pathway2compound[[path_id]]) {
        found_pathways <- c(found_pathways, path_id)
      }
    }

    # 如果有找到通路则添加到结果
    if (length(found_pathways) > 0) {
      for (path_id in found_pathways) {
        results[[length(results) + 1]] <- data.frame(
          keggID = compound_id,  # 列名改为keggID
          pathwayID = path_id,
          pathwayName = unname(pathway_names[path_id]),  # 安全访问名称
          stringsAsFactors = FALSE
        )
      }
    }

    # 更新进度条
    if (length(KEGGid) > 50) utils::setTxtProgressBar(pb, i)
  }

  # 关闭进度条
  if (length(KEGGid) > 50) close(pb)

  # 合并结果并排序
  if (length(results) > 0) {
    df <- do.call(rbind, results)
    df <- df[order(df$keggID, df$pathwayID), ]  # 按keggID排序
    rownames(df) <- NULL
    return(df)
  } else {
    message("None of the provided metabolites were found in any pathways.")
    return(data.frame(keggID = character(0),  # 列名改为keggID
                      pathwayID = character(0),
                      pathwayName = character(0)))
  }
}

