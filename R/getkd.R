#' Get and Cache KEGG Pathway Data
#'
#' Retrieves pathway data from KEGG and caches it for enrichment analysis.
#'
#' @param species KEGG species code (e.g., "hsa" for Homo sapiens)
#' @param cache_dir Directory for caching the downloaded data
#'
#' @return A list containing two elements:
#'   \item{pathways}{Named vector mapping pathway IDs to pathway names}
#'   \item{pathway2compound}{List mapping pathway IDs to vectors of compound IDs}
#'
#' @details
#' This function first checks if cached data exists for the species. If found,
#' loads cached data; otherwise downloads from KEGG and saves to cache.
#'
#' @examples
#' \dontrun{
#' # Get pathway data for human
#' kegg_data <- getkd("hsa", cache_dir = "./cache")
#' }
#'
#' @importFrom KEGGREST keggList keggGet
#' @export

getkd <- function(species, cache_dir) {

  # 创建缓存文件名
  cache_file <- file.path(cache_dir, paste0("kegg_pathways_", species, ".rds"))

  # 尝试从缓存加载通路数据
  if (file.exists(cache_file)) {
    message("正在加载缓存中的通路数据...")
    pathway_data <- readRDS(cache_file)
    pathways <- pathway_data$pathways
    pathway2compound <- pathway_data$pathway2compound
    message(sprintf("从缓存加载了 %d 条通路", length(pathway2compound)))
  } else {
    # 如果没有缓存则从KEGG获取数据
    message("正在从KEGG获取通路数据...")
    pathways <- KEGGREST::keggList("pathway", species)

    # 提取通路ID
    path_ids <- sub("path:", "", names(pathways))

    # 进度条：获取通路中的化合物
    pb <- utils::txtProgressBar(min = 0, max = length(path_ids), style = 3)
    pathway2compound <- list()

    for (i in seq_along(path_ids)) {
      pid <- path_ids[i]
      entry <- tryCatch(KEGGREST::keggGet(pid)[[1]], error = function(e) NULL)

      # 如果存在化合物信息，添加到映射中
      if (!is.null(entry) && !is.null(entry$COMPOUND)) {
        pathway2compound[[pid]] <- names(entry$COMPOUND)
      }
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)

    # 保存数据到缓存
    pathway_data <- list(
      pathways = pathways,
      pathway2compound = pathway2compound
    )
    saveRDS(pathway_data, cache_file)
    message(sprintf("已将 %d 条通路缓存到 %s", length(pathway2compound), cache_file))
  }

  return(list(
    pathways = pathways,
    pathway2compound = pathway2compound
  ))
}
