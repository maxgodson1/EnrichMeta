# enrichmeta: KEGG代谢物富集分析工具包

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## 简介

enrichmeta是一个用于KEGG代谢通路富集分析的R包，提供从数据获取、富集分析到结果可视化的完整工作流。

主要功能：
- 自动获取并缓存KEGG通路数据
- 代谢物富集分析（基于超几何检验）
- 富集结果可视化（点图、条形图）
- 通路间共享代谢物分析
- 代谢物/通路ID与名称的相互转换

## 安装

```r
# 从GitHub安装
devtools::install_github("maxgodson1/enrichmeta")

# 加载包
library(enrichmeta)
```

## 详细使用说明

### 1. 获取KEGG通路数据

使用`getkeggdata()`函数获取并缓存物种特异的KEGG通路数据：

```r
# 获取斑马鱼(dre)的KEGG通路数据
kegg_data <- getkeggdata("dre", cache_dir = "./kegg_cache")

# 查看数据结构
str(kegg_data, max.level = 1)

# 获取人类(hsa)数据
# kegg_data <- getkeggdata("hsa", cache_dir = "./kegg_cache")
```

### 2. 代谢物富集分析

使用`enrichkegg()`对显著代谢物进行富集分析：

```r
# 定义显著代谢物列表
significant_metabolites <- c("C00187", "C00245", "C00641")

# 执行富集分析
enrich_results <- enrichkegg(
  keggid = significant_metabolites,
  species = "dre",
  pathway_data = kegg_data
)

# 查看完整结果
print(enrich_results)

# 查看显著富集通路 (p.adjust < 0.05)
significant_pathways <- subset(enrich_results, p.adjust < 0.05)
print(significant_pathways)
```

### 3. 可视化富集结果

#### 3.1 点图可视化

```r
plotenrichdot(
  results = enrich_results,
  top = 15,
  title = "斑马鱼代谢通路富集分析",
  low_color = "red",
  point_size_range = c(3, 10),
  base_size = 14
)
```

#### 3.2 条形图可视化

```r
plotenrichbar(
  results = enrich_results,
  top = 10,
  title = "Top富集通路",
  low_color = "darkred",
  bar_width = 0.7,
  base_size = 14
)
```

### 4. 通路间共享代谢物分析

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

# 查看共享详情
print(shared_result$shared_df)
```

### 5. 代谢物-通路映射

#### 5.1 根据代谢物查找通路

```r
metabolite_pathways <- findcpdspaths(
  c("C00031", "C00221", "C00118"),
  pathway_data = kegg_data
)

# 查看结果
print(metabolite_pathways)
```

#### 5.2 根据通路ID查找名称

```r
path_ids <- c("dre00010", "dre00020", "dre00030")
path_names <- findpathsnames(path_ids, kegg_data)
print(path_names)
```

#### 5.3 根据代谢物ID查找名称

```r
kegg_ids <- c("C00187", "C00245", "C00641")
compound_names <- findcpdsnames(kegg_ids, kegg_data)
print(compound_names)
```

### 6. 完整工作流示例

```r
# 步骤1: 获取数据
kegg_data <- getkeggdata("hsa", cache_dir = "./kegg_cache")

# 步骤2: 定义显著代谢物
sig_met <- c("C00031", "C00033", "C00022", "C00041", "C00149")

# 步骤3: 富集分析
enrich_results <- enrichkegg(
  keggid = sig_met,
  species = "hsa",
  pathway_data = kegg_data
)

# 步骤4: 可视化
plotenrichdot(enrich_results, title = "人类代谢通路富集分析")

# 步骤5: 分析共享代谢物
selected_pathways <- head(enrich_results$pathwayID, 5)
filtered_data <- list(
  pathways = kegg_data$pathways[selected_pathways],
  pathscpds = kegg_data$pathscpds[selected_pathways]
)
shared_df <- findsharedcpds(filtered_data, min_shared = 1)$shared_df
print(shared_df)

# 步骤6: 获取代谢物名称
compound_names <- findcpdsnames(sig_met, kegg_data)
print(compound_names)
```

## 函数参考

### 核心函数

1. `getkeggdata(species, cache_dir)`
   - 获取并缓存KEGG通路数据
   - 参数:
     - `species`: KEGG物种代码 (如"hsa"表示人类)
     - `cache_dir`: 缓存目录路径

2. `enrichkegg(keggid, species, pathway_data, p.adjust.method = "BH")`
   - 执行KEGG富集分析
   - 参数:
     - `keggid`: 显著代谢物KEGG ID向量
     - `species`: KEGG物种代码
     - `pathway_data`: 从getkeggdata获取的数据
     - `p.adjust.method`: p值校正方法

3. `plotenrichdot(results, top = 25, ...)`
   - 创建富集结果点图

4. `plotenrichbar(results, top = 25, ...)`
   - 创建富集结果条形图

### 辅助函数

5. `findsharedcpds(pathway_data, min_shared = 1)`
   - 分析通路间共享代谢物

6. `findcpdspaths(keggid, pathway_data)`
   - 根据代谢物查找所属通路

7. `findpathsnames(pathwayid, pathways_data)`
   - 根据通路ID查找名称

8. `findcpdsnames(keggid, pathways_data)`
   - 根据代谢物ID查找名称

## 常见问题解答

### Q1: 如何获取物种代码？
KEGG物种代码通常是3-4个字母，例如：
- 人类: hsa
- 小鼠: mmu
- 大鼠: rno
- 斑马鱼: dre
- 大肠杆菌: eco

### Q2: 富集分析结果中没有通路怎么办？
可能原因：
1. 代谢物ID格式不正确（应使用"C00031"格式，不含"cpd:"前缀）
2. 物种代码不匹配
3. 代谢物不在任何通路中

### Q3: 如何获取代谢物的KEGG ID？
通常从代谢组学分析软件输出中获得，或使用KEGG数据库查询：
https://www.kegg.jp/kegg/compound/

## 支持与贡献

如有问题或建议，请提交至GitHub issue：
https://github.com/maxgodson1/enrichmeta/issues

欢迎贡献代码！请通过pull request提交改进：
https://github.com/maxgodson1/enrichmeta/pulls

## 引用

如果您在研究中使用了enrichmeta，请引用：

> Zhang, B. (2025). enrichmeta: KEGG Metabolite Enrichment Analysis in R. https://github.com/maxgodson1/enrichmeta

## 许可证

本项目采用GPL-3许可证 - 详见[LICENSE](https://github.com/maxgodson1/enrichmeta/blob/master/LICENSE)文件。
