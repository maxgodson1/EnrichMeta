#' Map Pathway IDs to Pathway Names
#'
#' Retrieves pathway names for given pathway IDs using cached KEGG data.
#'
#' @param pathwayid Character vector of KEGG pathway IDs (e.g., c("hsa00010", "hsa00020"))
#' @param pathways_data Pathway data object from `getkeggdata()`
#'
#' @return A data frame with two columns:
#'   \itemize{
#'     \item \code{pathwayID}: Input pathway ID
#'     \item \code{Description}: Pathway name (or "Not found" if ID not in data)
#'   }
#'
#' @examples
#' \dontrun{
#' kegg_data <- getkeggdata("hsa", cache_dir = "./cache")
#' findpathsnames(c("hsa00010", "hsa00020"), pathways_data = kegg_data)
#' }
#' @export

findpathsnames <- function(pathwayid, pathways_data) {
  # 确保输入是字符向量
  if (!is.character(pathwayid)) {
    stop("pathwayid must be a character vector")
  }

  # 创建结果数据框
  result <- data.frame(
    pathwayID = character(length(pathwayid)),
    Description = character(length(pathwayid)),
    stringsAsFactors = FALSE
  )

  # 填充结果
  for (i in seq_along(pathwayid)) {
    id <- pathwayid[i]

    # 检查ID是否在通路数据中
    if (id %in% names(pathways_data$pathways)) {
      result$pathwayID[i] <- id
      result$Description[i] <- pathways_data$pathways[[id]]
    } else {
      result$pathwayID[i] <- id
      result$Description[i] <- "Not found"
    }
  }

  return(result)
}
