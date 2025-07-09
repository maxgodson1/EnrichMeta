#' Perform KEGG metabolite enrichment with caching and progress tracking
#'
#' @param KEGGid A character vector of KEGG compound IDs (e.g., "C00160")
#' @param species A character string, KEGG species code (e.g., "hsa", "eco")
#' @param p.adjust.method Adjustment method for multiple testing, default is "BH"
#' @param use_cache Logical, whether to use cached pathway data, default TRUE
#' @param cache_dir Directory to store cache files, default is tempdir()
#'
#' @return A data.frame with enrichment results
#' @export
KeggEnrich <- function(KEGGid, species, p.adjust.method = "BH", 
                            use_cache = TRUE, cache_dir = tempdir()) {
  
  # 创建缓存文件名
  cache_file <- file.path(cache_dir, paste0("kegg_pathways_", species, ".rds"))
  
  # 尝试从缓存加载通路数据
  if (use_cache && file.exists(cache_file)) {
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
    if (use_cache) {
      pathway_data <- list(
        pathways = pathways,
        pathway2compound = pathway2compound
      )
      saveRDS(pathway_data, cache_file)
      message(sprintf("Cached %d pathways to %s", length(pathway2compound), cache_file))
    }
  }
  
  # 准备富集分析
  all_cmpds <- unique(unlist(pathway2compound))
  N <- length(KEGGid)
  M <- length(all_cmpds)
  
  # 创建ID到名称的映射（更健壮的方法）
  id_to_name <- setNames(
    unname(pathways), 
    sub("path:", "", names(pathways))
  )
  
  # 预分配结果列表
  results_list <- vector("list", length(pathway2compound))
  valid_count <- 0
  
  # 进度条：执行富集分析
  message("\nPerforming enrichment analysis...")
  pb <- utils::txtProgressBar(min = 0, max = length(pathway2compound), style = 3)
  
  for (i in seq_along(pathway2compound)) {
    pid <- names(pathway2compound)[i]
    pw_cmpds <- pathway2compound[[pid]]
    overlap <- intersect(KEGGid, pw_cmpds)
    k <- length(overlap)
    n <- length(pw_cmpds)
    
    if (k > 0) {
      valid_count <- valid_count + 1
      p <- phyper(k - 1, n, M - n, N, lower.tail = FALSE)
      enrich_ratio <- (k * M) / (n * N)
      
      # 更健壮的通路名称获取方法
      pathway_name <- if (!is.null(id_to_name[[pid]])) {
        id_to_name[[pid]]
      } else {
        # 如果无法找到名称，使用ID作为后备
        pid
      }
      
      results_list[[valid_count]] <- data.frame(
        pathwayID = pid,
        Description = pathway_name,
        MetaRatio = paste0(k, "/", N),
        BgRatio = paste0(n, "/", M),
        pvalue = p,
        EnrichmentRatio = enrich_ratio,
        keggID = paste(overlap, collapse = "; "),
        Count = k,
        stringsAsFactors = FALSE
      )
    }
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # 合并结果
  if (valid_count > 0) {
    results <- do.call(rbind, results_list[1:valid_count])
    results$p.adjust <- p.adjust(results$pvalue, method = p.adjust.method)
    results <- results[order(results$p.adjust), ]
    rownames(results) <- 1:nrow(results)
  } else {
    results <- data.frame(
      pathwayID = character(),
      Description = character(),
      MetaRatio = character(),
      BgRatio = character(),
      pvalue = numeric(),
      EnrichmentRatio = numeric(),
      keggID = character(),
      Count = integer(),
      p.adjust = numeric(),
      stringsAsFactors = FALSE
    )
    message("No significant enrichment found.")
  }
  
  return(results)
}
