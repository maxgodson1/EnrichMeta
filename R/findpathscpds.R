#' Retrieve Metabolites in KEGG Pathways
#'
#' Retrieves all metabolites associated with specified KEGG pathway(s).
#'
#' @param pathwayid Character vector of KEGG pathway IDs (e.g., c("hsa00010", "hsa00020"))
#' @param pathway_data Pathway data object from \code{\link{getkeggdata}}
#'
#' @return A data frame with two columns:
#'   \itemize{
#'     \item \code{pathwayID}: KEGG pathway ID
#'     \item \code{keggID}: Metabolite KEGG IDs in the pathway
#'   }
#'   Each row represents a metabolite in a pathway.
#'
#' @details
#' This function retrieves all metabolites associated with the specified KEGG pathways.
#' If a pathway ID is not found in the pathway data, a warning will be issued for that ID.
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
#' # Get metabolites in pathways
#' metabolites <- findpathscpds(c("hsa00010", "hsa00020"), pathway_data = kegg_data)
#'
#' # View metabolites in glycolysis pathway
#' subset(metabolites, pathwayID == "hsa00010")
#' }
#' @export

findpathscpds <- function(pathwayid, pathway_data) {
  # 验证必需参数
  if (missing(pathwayid)) stop("Missing required argument: pathwayid")
  if (missing(pathway_data)) stop("Missing required argument: pathway_data")

  # 检查pathway_data中必需的组成部分
  if (!"path_cpd_map" %in% names(pathway_data)) {
    stop("pathway_data must contain 'path_cpd_map' element")
  }

  # 初始化空的结果数据框
  results <- data.frame(
    pathwayID = character(),
    keggID = character(),
    stringsAsFactors = FALSE
  )

  # 提取通路-代谢物映射数据
  path_cpd_map <- pathway_data$path_cpd_map

  # 处理每个通路ID
  for (id in pathwayid) {
    # 检查通路是否存在于数据中
    if (id %in% names(path_cpd_map)) {
      metabolites <- path_cpd_map[[id]]

      # 如果当前通路中有代谢物，则添加它们
      if (length(metabolites) > 0) {
        # 为当前通路的代谢物创建数据框
        pathway_df <- data.frame(
          pathwayID = id,
          keggID = metabolites,
          stringsAsFactors = FALSE
        )

        # 添加到结果中
        results <- rbind(results, pathway_df)
      }
    } else {
      # 为缺失的通路ID发出警告
      warning(sprintf("Pathway '%s' not found in pathway_data", id))
    }
  }

  return(results)
}
