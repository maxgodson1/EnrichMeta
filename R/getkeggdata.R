#' Get and Cache KEGG Pathway Data
#'
#' Retrieves pathway data from KEGG and caches it for enrichment analysis. This function
#' automatically handles caching to reduce repeated KEGG API calls.
#'
#' @param species KEGG species code (e.g., "hsa" for Homo sapiens, "dre" for zebrafish)
#' @param cache_dir Directory for caching the downloaded data
#'
#' @return A list containing three elements:
#'   \item{pathways}{Named vector mapping pathway IDs to pathway names (e.g. list("hsa00010" = "Glycolysis"))}
#'   \item{pathscpds}{List mapping pathway IDs to vectors of compound IDs (e.g. list("hsa00010" = c("C00022", "C00031")))}
#'   \item{compounds}{Named list mapping compound IDs to compound names (e.g. list("C00022" = "Pyruvate; Pyruvic acid"))}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Checks for existing cache in \code{cache_dir}
#'   \item If cached data exists, loads and returns it
#'   \item If no cache, queries KEGG API for:
#'     \itemize{
#'       \item All pathways for the species
#'       \item Compounds in each pathway
#'       \item Compound name mappings
#'     }
#'   \item Saves results to cache for future use
#' }
#'
#' @note
#' Compound names are returned as semicolon-separated strings containing all known synonyms.
#' API queries may take 2-5 minutes for species with many pathways (e.g., human has 400+ pathways).
#'
#' @examples
#' \dontrun{
#' # Get pathway data for human
#' kegg_data <- getkeggdata("hsa", cache_dir = "./kegg_cache")
#'
#' # Access glycolysis compounds
#' glycolysis_cpds <- kegg_data$pathscpds$hsa00010
#'
#' # Get pyruvate names
#' kegg_data$compounds$C00022
#' # Returns: "Pyruvate; Pyruvic acid; 2-Oxopropanoate; ..."
#' }
#'
#' @importFrom KEGGREST keggList keggGet
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export

getkeggdata <- function(species, cache_dir) {

  # 确保缓存目录存在
  if (!dir.exists(cache_dir)) {
    message(sprintf("创建缓存目录: %s", cache_dir))
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # 创建缓存文件名
  cache_file <- file.path(cache_dir, paste0("kegg_pathways_", species, ".rds"))

  # 尝试从缓存加载
  if (file.exists(cache_file)) {
    message("正在加载缓存中的通路数据...")
    pathway_data <- readRDS(cache_file)
    message(sprintf("加载了 %d 条通路和 %d 个化合物",
                    length(pathway_data$pathscpds),
                    length(pathway_data$compounds)))
  } else {
    # 从KEGG获取新数据
    message("正在从KEGG获取通路数据(首次加载可能需要较长的时间)...")
    pathways_vector <- KEGGREST::keggList("pathway", species)

    # 提取通路ID
    path_ids <- sub("path:", "", names(pathways_vector))

    # 获取通路中的化合物
    message("获取通路中的化合物...")
    pb <- utils::txtProgressBar(min = 0, max = length(path_ids), style = 3)
    pathscpds <- list()

    for (i in seq_along(path_ids)) {
      pid <- path_ids[i]
      entry <- tryCatch(KEGGREST::keggGet(pid)[[1]], error = function(e) NULL)

      if (!is.null(entry) && !is.null(entry$COMPOUND)) {
        # 移除"cpd:"前缀
        compound_ids <- sub("cpd:", "", names(entry$COMPOUND))
        pathscpds[[pid]] <- compound_ids
      }
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)

    # 转换为命名列表
    pathways_list <- as.list(pathways_vector)
    names(pathways_list) <- path_ids

    # 获取所有唯一化合物ID
    all_compound_ids <- unique(unlist(pathscpds, use.names = FALSE))

    # 创建化合物ID->名称映射（逐个获取）
    compounds <- list()

    if (length(all_compound_ids) > 0) {
      message("获取化合物名称映射...")
      pb <- utils::txtProgressBar(min = 0, max = length(all_compound_ids), style = 3)

      for (i in seq_along(all_compound_ids)) {
        cid <- all_compound_ids[i]

        # 逐个获取化合物信息
        result <- tryCatch(
          KEGGREST::keggGet(paste0("cpd:", cid))[[1]],
          error = function(e) NULL
        )

        # 提取名称
        if (!is.null(result)) {
          compounds[[cid]] <- if (!is.null(result$NAME)) {
            paste(result$NAME, collapse = "; ")
          } else {
            cid
          }
        } else {
          compounds[[cid]] <- cid
        }

        # 更新进度条
        utils::setTxtProgressBar(pb, i)
      }
      close(pb)
    } else {
      message("没有化合物需要处理")
    }

    # 构建返回对象
    pathway_data <- list(
      pathways = pathways_list,  # 通路ID->名称（列表）
      pathscpds = pathscpds,     # 通路ID->化合物ID列表（列表）
      compounds = compounds       # 化合物ID->名称（列表）
    )

    # 保存缓存
    saveRDS(pathway_data, cache_file)
    message(sprintf("已缓存 %d 条通路和 %d 个化合物",
                    length(pathscpds),
                    length(compounds)))
  }

  return(pathway_data)
}
