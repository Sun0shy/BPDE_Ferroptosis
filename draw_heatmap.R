## 加载热图的R包
library(pheatmap)
## 加载颜色包
library(viridisLite)
## 加载管道符
library(dplyr)

draw_heatmap <- function(allDiff,p_value,logfc,group_heat) {
  p_value = as.numeric(p_value)
  logfc = as.numeric(logfc)
  ## 传入标准筛选差异基因
  diffgene <- allDiff %>%
    filter(P.Value < p_value) %>%
    filter(abs(logFC) >logfc)
  
  heatdata <- diff_data[rownames(diffgene),]
  
  ### 将分组信息作为数据框
  annotation_col <- data.frame(group_heat)
  ### 一列是样本一列数分组，这样才能将样本和分组对应起来
  rownames(annotation_col) <- colnames(heatdata)

  ## 绘图
  #### 如果注释出界，可以通过调整格子比例和字体修正
  hot_plot <- pheatmap(heatdata, #热图的数据
                       cluster_rows = TRUE,#行聚类
                       cluster_cols = TRUE,#列聚类，可以看出样本之间的区分度
                       annotation_col =annotation_col, #标注样本分类
                       annotation_legend=TRUE, # 显示注释
                       show_rownames = T,# 显示行名
                       scale = "row", #以行来标准化，这个功能很不错
                       color = viridis(10, alpha = 1, begin = 0.5, end = 1, direction = 1),#调色
                       #filename = "heatmap_F.pdf",#是否保存
                       cellwidth = 20, # 格子宽度
                       cellheight = 5,# 格子高度
                       fontsize = 5, # 字体大小
                       treeheight_col	=20 # 列聚类树高度
  )
  
  return(hot_plot)
  
}




# # 传入用于做热图的数据
# ## 传入所有差异分析的数据
# allDiff
# ## 传入标准筛选差异基因
# diffgene <- allDiff %>%
#   filter(P.Value < 0.05) %>%
#   filter(abs(logFC) >0.7)
# heatdata <- diff_data[rownames(diffgene),]
# ## 加入注释信息
# ### 设置分组信息，按照样本的组别一个一个指定是那个组的
# group_heat
# ### 将分组信息作为数据框
# annotation_col <- data.frame(group_heat)
# ### 一列是样本一列数分组，这样才能将样本和分组对应起来
# rownames(annotation_col) <- colnames(heatdata)
# 
# ## 绘图
# #### 如果注释出界，可以通过调整格子比例和字体修正
# hot_plot <- pheatmap(heatdata, #热图的数据
#                      cluster_rows = TRUE,#行聚类
#                      cluster_cols = TRUE,#列聚类，可以看出样本之间的区分度
#                      annotation_col =annotation_col, #标注样本分类
#                      annotation_legend=TRUE, # 显示注释
#                      show_rownames = T,# 显示行名
#                      scale = "row", #以行来标准化，这个功能很不错
#                      color = viridis(10, alpha = 1, begin = 0.5, end = 1, direction = 1),#调色
#                      #filename = "heatmap_F.pdf",#是否保存
#                      cellwidth = 20, # 格子宽度
#                      cellheight = 5,# 格子高度
#                      fontsize = 5, # 字体大小
#                      treeheight_col	=20 # 列聚类树高度
# )













