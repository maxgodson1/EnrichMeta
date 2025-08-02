# enrichmeta: Metabolite Enrichment Analysis Using KEGG Pathways

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

`enrichmeta` is an R package for metabolite enrichment analysis using KEGG pathways. It provides a comprehensive toolkit for:
- Retrieving and caching KEGG pathway data
- Performing metabolite set enrichment analysis
- Analyzing shared metabolites between pathways
- Visualizing enrichment results
- Mapping metabolites to pathways and vice versa

## Installation

You can install the development version from GitHub:
```r
# install.packages("devtools")
devtools::install_github("maxgodson1/enrichmeta")
```

## Quick Start

```r
library(enrichmeta)

# Step 1: Get KEGG data for a species (e.g., human 'hsa')
kegg_data <- getkeggdata("hsa", cache_dir = "./kegg_cache")

# Step 2: Perform enrichment analysis
enrich_results <- enrichkegg(
  keggid = c("C00022", "C00031", "C00041", "C00033", "C00036"),
  pathway_data = kegg_data
)

# View top 10 enriched pathways
head(enrich_results, 10)

# Step 3: Visualize results
plotenrichbar(enrich_results, top = 15, title = "Top Enriched Pathways")
plotenrichdot(enrich_results, top = 15, title = "Enriched Pathways Dot Plot")
```

## Vignette

For a detailed tutorial, see the package vignette:
```r
vignette("tutorial", package = "enrichmeta")
```

## Documentation

Full function documentation is available at:
```r
help(package = "enrichmeta")
```

## Citation

To cite `enrichmeta` in publications, please use:
> Zhang, B. (2025). enrichmeta: Metabolite Enrichment Analysis Using KEGG Pathways. R package version 1.1.1.

## Issues and Contributions

Please report any issues at: [https://github.com/maxgodson1/enrichmeta/issues](https://github.com/maxgodson1/enrichmeta/issues)

Contributions via pull requests are welcome!

