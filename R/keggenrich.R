#' Perform KEGG Metabolite Enrichment Analysis
#'
#' Conducts KEGG metabolite enrichment analysis using hypergeometric test
#'
#' @param KEGGid Vector of metabolite KEGG IDs (e.g., "C00160")
#' @param species KEGG species code (e.g., "hsa")
#' @param pathway_data Pathway data obtained from getkd()
#' @param p.adjust.method p-value adjustment method (default="BH")
#'
#' @return Enrichment results data frame with the following columns:
#'   \item{pathwayID}{KEGG pathway ID (e.g., "dre04744")}
#'   \item{Description}{Pathway description}
#'   \item{PathwaySize}{Number of metabolites in the pathway (background)}
#'   \item{Count}{Number of significant metabolites in pathway}
#'   \item{MetaRatio}{Ratio of significant metabolites in pathway (e.g., "2/3")}
#'   \item{BgRatio}{Background ratio (e.g., "6/203")}
#'   \item{pvalue}{Raw p-value from hypergeometric test}
#'   \item{p.adjust}{Adjusted p-value}
#'   \item{EnrichmentRatio}{Enrichment ratio (fold enrichment)}
#'   \item{keggID}{List of significant metabolite IDs in the pathway}
#'   
#' @export

keggenrich <- function(KEGGid, species, pathway_data, p.adjust.method = "BH") {
  
  # 获取通路数据
  pathways <- pathway_data$pathways
  pathcompounds <- pathway_data$pathcompounds
  
  # 准备富集分析
  all_cmpds <- unique(unlist(pathcompounds))
  N <- length(KEGGid)  # 输入的代谢物数量
  M <- length(all_cmpds)  # 背景代谢物总数
  
  # 预分配结果列表
  results_list <- vector("list", length(pathcompounds))
  valid_count <- 0
  
  # 进度条：执行富集分析
  message("\nPerforming enrichment analysis...")
  pb <- utils::txtProgressBar(min = 0, max = length(pathcompounds), style = 3)
  
  for (i in seq_along(pathcompounds)) {
    pid <- names(pathcompounds)[i]
    pw_cmpds <- pathcompounds[[pid]]
    overlap <- intersect(KEGGid, pw_cmpds)
    k <- length(overlap)  # 显著代谢物在通路中的数量
    n <- length(pw_cmpds)  # 通路中的代谢物总数
    
    if (k > 0) {
      valid_count <- valid_count + 1
      p <- phyper(k - 1, n, M - n, N, lower.tail = FALSE)
      enrich_ratio <- (k * M) / (n * N)
      
      # 获取通路名称
      pathway_name <- if (!is.null(pathways[[pid]])) {
        pathways[[pid]]
      } else {
        pid
      }
      
      results_list[[valid_count]] <- data.frame(
        pathwayID = pid,
        Description = pathway_name,
        PathwaySize = n,  # 通路中的代谢物总数
        Count = k,        # 显著代谢物数量（移动到PathwaySize后面）
        MetaRatio = paste0(k, "/", N),
        BgRatio = paste0(n, "/", M),
        pvalue = p,
        EnrichmentRatio = enrich_ratio,
        keggID = paste(overlap, collapse = "; "),
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
    
    # 按p值排序并重新编号行
    results <- results[order(results$p.adjust), ]
    
    # 重新排序列顺序（将p.adjust移到pvalue后面）
    results <- results[, c("pathwayID", "Description", "PathwaySize", "Count", 
                           "MetaRatio", "BgRatio", "pvalue", "p.adjust", 
                           "EnrichmentRatio", "keggID")]
    
    rownames(results) <- 1:nrow(results)
  } else {
    results <- data.frame(
      pathwayID = character(),
      Description = character(),
      PathwaySize = integer(),
      Count = integer(),
      MetaRatio = character(),
      BgRatio = character(),
      pvalue = numeric(),
      p.adjust = numeric(),
      EnrichmentRatio = numeric(),
      keggID = character(),
      stringsAsFactors = FALSE
    )
    message("No significant enrichment found.")
  }
  
  return(results)
}
