#' Get and Cache KEGG Pathway Data
#'
#' Retrieves pathway data from KEGG and caches it for enrichment analysis.
#'
#' @param species KEGG species code (e.g., "hsa" for Homo sapiens, "dre" for zebrafish)
#' @param cache_dir Directory for caching the downloaded data
#'
#' @return A list containing two elements:
#'   \item{pathways}{Named list mapping pathway IDs to pathway names (accessible via `$`)}
#'   \item{pathcompounds}{List mapping pathway IDs to vectors of compound IDs}
#'
#' @details
#' This function first checks if cached data exists for the species. If found,
#' loads cached data; otherwise downloads from KEGG and saves to cache.
#'
#' The returned object has two components:
#' \itemize{
#'   \item \code{pathways}: Named list where names are pathway IDs (e.g., "dre04744")
#'         and values are pathway names
#'   \item \code{pathcompounds}: List where each element is a vector of compound IDs
#'         (e.g., "cpd:C00187") for a pathway
#' }
#'
#' @examples
#' \dontrun{
#' # Get pathway data for zebrafish
#' kegg_data <- getkd("dre", cache_dir = "./cache")
#'
#' # Access pathway names
#' pathway_name <- kegg_data$pathways$dre04744
#'
#' # Access compounds in a pathway
#' compounds <- kegg_data$pathcompounds$dre04744
#' }
#'
#' @importFrom KEGGREST keggList keggGet
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export

getkd <- function(species, cache_dir) {

  # 确保缓存目录存在
  if (!dir.exists(cache_dir)) {
    message(sprintf("创建缓存目录: %s", cache_dir))
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # 创建缓存文件名
  cache_file <- file.path(cache_dir, paste0("kegg_pathways_", species, ".rds"))

  # 尝试从缓存加载通路数据
  if (file.exists(cache_file)) {
    message("正在加载缓存中的通路数据...")
    pathway_data <- readRDS(cache_file)
    pathways_vector <- pathway_data$pathways  # 重命名变量以区分
    pathcompounds <- pathway_data$pathcompounds
    message(sprintf("从缓存加载了 %d 条通路", length(pathcompounds)))
  } else {
    # 如果没有缓存则从KEGG获取数据
    message("正在从KEGG获取通路数据...")
    pathways_vector <- KEGGREST::keggList("pathway", species)  # 重命名变量

    # 提取通路ID
    path_ids <- sub("path:", "", names(pathways_vector))

    # 进度条：获取通路中的化合物
    pb <- utils::txtProgressBar(min = 0, max = length(path_ids), style = 3)
    pathcompounds <- list()

    for (i in seq_along(path_ids)) {
      pid <- path_ids[i]
      entry <- tryCatch(KEGGREST::keggGet(pid)[[1]], error = function(e) NULL)

      # 如果存在化合物信息，添加到映射中
      if (!is.null(entry) && !is.null(entry$COMPOUND)) {
        pathcompounds[[pid]] <- names(entry$COMPOUND)
      }
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)

    # 将通路向量转换为可直接用$访问的命名列表
    path_ids <- sub("path:", "", names(pathways_vector))
    pathways_list <- as.list(pathways_vector)
    names(pathways_list) <- path_ids

    # 保存数据到缓存
    pathway_data <- list(
      pathways = pathways_list,
      pathcompounds = pathcompounds
    )
    saveRDS(pathway_data, cache_file)
    message(sprintf("已将 %d 条通路缓存到 %s", length(pathcompounds), cache_file))
  }

  return(pathway_data)

}
