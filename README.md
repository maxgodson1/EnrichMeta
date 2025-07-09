# enrichmeta
KEGG enrichment analysis and visualization tools

## 安装
```r
devtools::install_github("maxgodson1/enrichmeta")
```

## 使用示例
```r
library(enrichmeta)

# 富集分析
results <- keggenrich(c("C00031", "C00022"), species = "hsa")

# 绘制条形图
enrichbar(results)

# 绘制点图
enrichdot(results)
```
