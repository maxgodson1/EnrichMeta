getkpd <- function(species,cache_dir) {
  
  # 创建缓存文件名
  cache_file <- file.path(cache_dir, paste0("kegg_pathways_", species, ".rds"))
  
  # 尝试从缓存加载通路数据
  if (file.exists(cache_file)) {
    message("Loading cached pathway data...")
    pathway_data <- readRDS(cache_file)
    pathways <- pathway_data$pathways
    pathway2compound <- pathway_data$pathway2compound
    message(sprintf("Loaded %d pathways from cache", length(pathway2compound)))
  } else {
    # 没有缓存则从KEGG获取
    message("Fetching pathway data from KEGG...")
    pathways <- KEGGREST::keggList("pathway", species)
    
    # 提取通路ID
    path_ids <- sub("path:", "", names(pathways))
    
    # 进度条：获取通路中的化合物
    pb <- utils::txtProgressBar(min = 0, max = length(path_ids), style = 3)
    pathway2compound <- list()
    
    for (i in seq_along(path_ids)) {
      pid <- path_ids[i]
      entry <- tryCatch(KEGGREST::keggGet(pid)[[1]], error = function(e) NULL)
      
      if (!is.null(entry) && !is.null(entry$COMPOUND)) {
        pathway2compound[[pid]] <- names(entry$COMPOUND)
      }
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # 保存到缓存
    pathway_data <- list(
      pathways = pathways,
      pathway2compound = pathway2compound
    )
    saveRDS(pathway_data, cache_file)
    message(sprintf("Cached %d pathways to %s", length(pathway2compound), cache_file))
  }
}
