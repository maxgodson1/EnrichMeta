# KEGGEnrichment
KEGG enrichment analysis and visualization tools

## 安装
```r
devtools::install_github("maxgodson1/EnrichMeta")
```

## 使用示例
```r
library(EnrichMeta)

# 富集分析
results <- KeggEnrich(c("C00031", "C00022"), species = "hsa")

# 绘制条形图
EnrichBar(results)

# 绘制点图
EnrichDot(results)
```
