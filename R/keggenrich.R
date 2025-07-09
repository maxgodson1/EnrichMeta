#' 执行KEGG代谢物富集分析
#'
#' 使用超几何检验进行KEGG代谢物富集分析
#'
#' @param KEGGid 代谢物KEGG ID向量 (e.g., "C00160")
#' @param species KEGG物种代码 (e.g., "hsa")
#' @param pathway_data 获取的通路数据
#' @param p.adjust.method p值校正方法，默认"BH"
#'
#' @return 富集结果数据框
#' @export

keggenrich <- function(KEGGid, species, pathway_data, p.adjust.method = "BH") {
  
  # 获取通路数据
  pathways <- pathway_data$pathways
  pathway2compound <- pathway_data$pathway2compound
  
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
