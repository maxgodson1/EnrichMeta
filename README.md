# enrichmeta - KEGG Metabolite Enrichment Analysis Package

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

enrichmeta 是一个用于 KEGG 代谢物富集分析的 R 包，提供完整的分析工作流和可视化功能。

## 核心功能

- 获取并缓存 KEGG 通路数据
- 查询代谢物通路信息
- 执行代谢物富集分析
- 可视化富集结果
- 分析通路间共享代谢物

## 主要函数

- `getkd()`: 获取并缓存 KEGG 通路数据
- `metapathinfor()`: 查询代谢物的通路信息
- `keggenrich()`: 执行 KEGG 代谢物富集分析
- `enrichbar()`: 绘制富集条形图
- `enrichdot()`: 绘制富集点图
- `shared_compounds()`: 分析通路间共享代谢物

## 安装

```r
# 从 GitHub 安装开发版
if (!require("devtools")) install.packages("devtools")
devtools::install_github("maxgodson1/enrichmeta", build_vignettes = TRUE)
```

## 使用教程

```r
vignette("tutorial", package = "enrichmeta")
```

## 完整工作流示例

```r
library(enrichmeta)

# 1. 获取通路数据
kegg_data <- getkd("dre", "path/to/cache")

# 2. 查询代谢物通路信息
metabolites <- c("C00300", "C00245", "C00763", "C16308")
pathway_info <- metapathinfor(metabolites, kegg_data)

# 3. 执行富集分析
res <- keggenrich(metabolites, "dre", kegg_data)

# 4. 可视化结果
enrichbar(res)
enrichdot(res)

# 5. 分析共享代谢物
selected_pathways <- res$pathwayID[res$p.adjust < 0.05]
filtered_data <- list(
  pathways = kegg_data$pathways[selected_pathways],
  pathway2compound = kegg_data$pathway2compound[selected_pathways]
)
shared_result <- shared_compounds(filtered_data)
```

## 许可证
GPL-3
