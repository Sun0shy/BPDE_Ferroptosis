### 5.加载绘图包
library(factoextra)

sample_pca <- function(diff_data,group,title) {
  title_pca = as.character(title)
  # 主成分分析pca
  ## 行是样本列是基因
  ### 1.转置
  data <- t(diff_data)
  ### 2.NAN值变成0
  data[is.nan(data)] <- 0 
  ### 3.删除没有变异的基因，使用函数var
  data <-  data[,which(apply(data,2,var) != 0)]
  ### 4.分析
  res.pca <- prcomp(data, scale = TRUE)
  ### 5.绘制主成分分析图
  sample_pca <- fviz_pca_ind(
    # 数据
    res.pca,
    # 分组
    col.ind = group,
    # 是否加圈
    addEllipses = TRUE,
    # 标签名称
    legend.title = "Groups",
    # 图片标题
    title = title_pca,
    # 是否避免遮挡
    repel = TRUE,
    # 圈的大小更置信，还是更广
    ellipse.type = "confidence",
    # 点的尺寸
    pointsize = 2.5,
    # 圈的样式
    palette = 'jco'
    
  )
  return(sample_pca)
}