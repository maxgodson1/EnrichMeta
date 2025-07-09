#' Get and cache KEGG pathway data
#'
#' This function retrieves pathway data from KEGG and caches it for enrichment analysis.
#'
#' @param species KEGG species code (e.g., "hsa")
#' @param cache_dir Cache directory
#'
#' @return A list containing two elements:
#'   - pathways: mapping from pathway IDs to pathway names
#'   - pathway2compound: mapping from pathway IDs to compound ID lists
#' @export

getkd <- function(species, cache_dir) {
  
  # Create cache filename
  cache_file <- file.path(cache_dir, paste0("kegg_pathways_", species, ".rds"))
  
  # Try to load pathway data from cache
  if (file.exists(cache_file)) {
    message("Loading cached pathway data...")
    pathway_data <- readRDS(cache_file)
    pathways <- pathway_data$pathways
    pathway2compound <- pathway_data$pathway2compound
    message(sprintf("Loaded %d pathways from cache", length(pathway2compound)))
  } else {
    # Fetch from KEGG if no cache exists
    message("Fetching pathway data from KEGG...")
    pathways <- KEGGREST::keggList("pathway", species)
    
    # Extract pathway IDs
    path_ids <- sub("path:", "", names(pathways))
    
    # Progress bar: get compounds in pathways
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
    
    # Save to cache
    pathway_data <- list(
      pathways = pathways,
      pathway2compound = pathway2compound
    )
    saveRDS(pathway_data, cache_file)
    message(sprintf("Cached %d pathways to %s", length(pathway2compound), cache_file))
  }
  
  return(list(
    pathways = pathways,
    pathway2compound = pathway2compound
  ))
}
