```
# enrichmeta - KEGG Metabolite Enrichment Analysis Package

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This R package provides tools for performing KEGG metabolite enrichment analysis and visualizing the results.

## Features

- Retrieve and cache KEGG pathway data
- Perform metabolite enrichment analysis using hypergeometric test
- Visualize enrichment results with bar plots and dot plots
- Analyze shared metabolites between pathways

## Installation

​```r
# Install development version from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("yourusername/enrichmeta", build_vignettes = TRUE)
```

## Usage

See the detailed tutorial:

r

```
# View the included tutorial
vignette("tutorial", package = "enrichmeta")
```

## Functions

- `getkd()`: Get and cache KEGG pathway data
- `keggenrich()`: Perform KEGG metabolite enrichment analysis
- `enrichbar()`: Draw enrichment bar plot
- `enrichdot()`: Draw enrichment dot plot
- `shared_compounds()`: Calculate shared metabolites between pathways

## Example Workflow

r

```
library(enrichmeta)

# 1. Get KEGG pathway data
kegg_data <- getkd("dre", "path/to/cache")

# 2. Perform enrichment analysis
kegg_id <- c("C00300", "C00245", "C00763", "C16308")
res <- keggenrich(kegg_id, "dre", kegg_data)

# 3. Visualize results
enrichbar(res)
enrichdot(res)

# 4. Analyze shared compounds
selected_pathways <- res$pathwayID
filtered_data <- list(
  pathways = kegg_data$pathways[selected_pathways],
  pathway2compound = kegg_data$pathway2compound[selected_pathways]
)
shared_result <- shared_compounds(filtered_data)
```

## License

GPL-3