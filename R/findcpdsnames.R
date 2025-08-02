#' Map Compound IDs to Compound Names
#'
#' Retrieves compound names for given KEGG compound IDs using cached data.
#'
#' @param keggid Character vector of KEGG compound IDs (e.g., c("C00031", "C00022"))
#' @param pathway_data Pathway data object from \code{\link{getkeggdata}}
#'
#' @return A data frame with two columns:
#'   \itemize{
#'     \item \code{keggID}: Input compound ID
#'     \item \code{Name}: Compound name (or "Not found" if ID not in data)
#'   }
#'
#' @examples
#' \dontrun{
#' kegg_data <- getkeggdata("hsa", cache_dir = "./cache")
#' findcpdsnames(c("C00031", "C00022"), pathway_data = kegg_data)
#' }
#' @export

findcpdsnames <- function(keggid, pathway_data) {
  # 确保输入是字符向量
  if (!is.character(keggid)) {
    stop("keggid must be a character vector")
  }

  # 创建结果数据框
  result <- data.frame(
    keggID = character(length(keggid)),
    Name = character(length(keggid)),
    stringsAsFactors = FALSE
  )

  # 填充结果
  for (i in seq_along(keggid)) {
    id <- keggid[i]

    # 检查ID是否在化合物数据中
    if (id %in% names(pathway_data$compounds)) {
      result$keggID[i] <- id
      result$Name[i] <- pathway_data$compounds[[id]]
    } else {
      result$keggID[i] <- id
      result$Name[i] <- "Not found"
    }
  }

  return(result)
}
