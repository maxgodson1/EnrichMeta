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
enrichbar(results, 
          title = "Top Metabolic Pathways",
          low_color = "darkred",
          bar_width = 0.7)

# 绘制点图
enrichdot(results, 
          top = 15,
          point_size_range = c(4, 10),
          show_grid = FALSE)
```

## 主要函数
keggenrich()
执行KEGG代谢物富集分析，支持缓存以加速重复分析
参数:
KEGGid: KEGG化合物ID向量
species: KEGG物种代码 (如"hsa", "eco")
p.adjust.method: p值校正方法 (默认"BH")
use_cache: 是否使用缓存 (默认TRUE)
cache_dir: 缓存目录 (默认tempdir())

enrichbar()
绘制富集结果条形图
参数:
results: KeggEnrich()返回的结果
top: 显示顶部通路数量 (默认25)
title: 图表标题
low_color/high_color: 颜色梯度
bar_width: 条形宽度
base_size: 基础字体大小
wrap_width: 描述文本换行宽度
show_ratio: 是否显示富集比值
color_limits: 颜色范围
legend_position: 图例位置

enrichdot()
绘制富集结果点图
参数:
results: KeggEnrich()返回的结果
top: 显示顶部通路数量 (默认25)
title: 图表标题
low_color/high_color: 颜色梯度
point_size_range: 点大小范围
point_alpha: 点透明度
base_size: 基础字体大小
wrap_width: 描述文本换行宽度
color_limits: 颜色范围
legend_position: 图例位置
show_grid: 是否显示网格线
x_expansion: X轴扩展空间
padding_right: 右侧留白空间
