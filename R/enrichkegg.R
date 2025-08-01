#' Perform KEGG Metabolite Enrichment Analysis
#'
#' Conducts metabolite set enrichment analysis using KEGG pathways and hypergeometric test
#'
#' @param keggid Character vector of metabolite KEGG IDs (e.g., c("C00022", "C00031")).
#'                Must be in short format (without "cpd:" prefix).
#' @param species KEGG species code (e.g., "hsa" for human). Must match the species used
#'                in \code{pathway_data}.
#' @param pathway_data Pathway data object from \code{\link{getkeggdata}}
#' @param p.adjust.method p-value adjustment method (default="BH"). Options: "holm", "hochberg",
#'                        "hommel", "bonferroni", "BH", "BY", "fdr".
#'
#' @return A data frame with enrichment results containing these columns:
#'   \item{pathwayID}{KEGG pathway ID (e.g., "hsa00010")}
#'   \item{Description}{Pathway description}
#'   \item{PathwaySize}{Total metabolites in the pathway (background set size)}
#'   \item{Count}{Number of significant metabolites mapped to this pathway}
#'   \item{MetaRatio}{Ratio of significant metabolites (format: significant/total_input)}
#'   \item{BgRatio}{Background ratio (format: pathway_metabolites/total_background_metabolites)}
#'   \item{pvalue}{Raw p-value from hypergeometric test}
#'   \item{p.adjust}{Adjusted p-value using specified method}
#'   \item{EnrichmentRatio}{Enrichment ratio calculated as (Count/PathwaySize)/(length(keggid)/length(all_background_metabolites)}
#'   \item{keggID}{Semicolon-separated list of metabolite IDs}
#'
#' @details
#' The enrichment analysis uses the hypergeometric test with these parameters:
#' \deqn{p = 1 - \sum_{i=0}^{k-1} \frac{{\binom{n}{i} \binom{M-n}{N-i}}}{{\binom{M}{N}}}}
#' Where:
#' \itemize{
#'   \item \eqn{M}: Total metabolites in all pathways (background)
#'   \item \eqn{n}: Metabolites in target pathway
#'   \item \eqn{N}: Significant metabolites in input
#'   \item \eqn{k}: Overlap between pathway and significant metabolites
#' }
#'
#' @note
#' Only pathways with at least one significant metabolite (Count > 0) are returned.
#' Results are sorted by adjusted p-value in ascending order.
#'
#' @examples
#' \dontrun{
#' # After getting pathway data with getkeggdata()
#' sig_metabolites <- c("C00022", "C00031", "C00033")
#' enrich_results <- enrichkegg(
#'   keggid = sig_metabolites,
#'   species = "hsa",
#'   pathway_data = kegg_data
#' )
#'
#' # View top 5 enriched pathways
#' head(enrich_results, 5)
#' }
#' @export

enrichkegg <- function(keggid, species, pathway_data, p.adjust.method = "BH") {

  # 获取通路数据
  pathways <- pathway_data$pathways
  pathscpds <- pathway_data$pathscpds

  # 准备富集分析
  all_cmpds <- unique(unlist(pathscpds))
  N <- length(keggid)  # 输入的代谢物数量
  M <- length(all_cmpds)  # 背景代谢物总数

  # 预分配结果列表
  results_list <- vector("list", length(pathscpds))
  valid_count <- 0

  # 进度条：执行富集分析
  message("\nPerforming enrichment analysis...")
  pb <- utils::txtProgressBar(min = 0, max = length(pathscpds), style = 3)

  for (i in seq_along(pathscpds)) {
    pid <- names(pathscpds)[i]
    pw_cmpds <- pathscpds[[pid]]
    overlap <- intersect(keggid, pw_cmpds)
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
