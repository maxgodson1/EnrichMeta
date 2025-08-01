# enrichmeta: KEGG代谢物富集分析工具包

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/maxgodson1/enrichmeta/workflows/R-CMD-check/badge.svg)](https://github.com/maxgodson1/enrichmeta/actions)

## 简介

enrichmeta是一个用于KEGG代谢通路富集分析的R包，提供从数据获取、富集分析到结果可视化的完整工作流。

主要功能：
- 🧬 自动获取并缓存KEGG通路数据
- 📊 代谢物富集分析（基于超几何检验）
- 📈 富集结果可视化（点图、条形图）
- 🔗 通路间共享代谢物分析
- 🔄 代谢物/通路ID与名称的相互转换

## 安装

```r
# 从GitHub安装
devtools::install_github("yourusername/enrichmeta")
```

## 快速开始

### 1. 获取KEGG数据

```r
library(enrichmeta)

# 获取斑马鱼(dre)的KEGG通路数据
kegg_data <- getkeggdata("dre", cache_dir = "./kegg_cache")
```

### 2. 富集分析

```r
# 定义显著代谢物
sig_met <- c("C00187", "C00245", "C00641")

# 执行富集分析
enrich_results <- enrichkegg(
  keggid = sig_met,
  species = "dre",
  pathway_data = kegg_data
)

# 查看结果
head(enrich_results)
```

### 3. 结果可视化

```r
# 点图
plotenrichdot(enrich_results, title = "斑马鱼代谢通路富集分析")

# 条形图
plotenrichbar(enrich_results, top = 10, title = "Top富集通路")
```

### 4. 共享代谢物分析

```r
# 选择感兴趣的代谢通路
selected_pathways <- c("dre00010", "dre00020", "dre00030")

# 筛选数据
filtered_data <- list(
  pathways = kegg_data$pathways[selected_pathways],
  pathscpds = kegg_data$pathscpds[selected_pathways]
)

# 计算共享代谢物
shared_result <- findsharedcpds(filtered_data, min_shared = 1)
print(shared_result$shared_df)
```

## 文档

- 详细教程: [vignettes/tutorial.html](https://maxgodson1.github.io/enrichmeta/articles/tutorial.html)
- 函数参考: [参考手册](https://maxgodson1.github.io/enrichmeta/reference/)

## 引用

如果您在研究中使用了enrichmeta，请引用：

> Zhang, B. (2023). enrichmeta: KEGG Metabolite Enrichment Analysis in R. https://github.com/maxgodson1/enrichmeta

## 贡献

欢迎通过issue或pull request贡献代码：
https://github.com/maxgodson1/enrichmeta

## 许可证

本项目采用GPL-3许可证 - 详见[LICENSE](https://github.com/maxgodson1/enrichmeta/blob/master/LICENSE)文件。
