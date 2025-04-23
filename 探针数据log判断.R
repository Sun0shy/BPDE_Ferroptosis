is_probe_log <- function(data) {
  # 数据第一列变成行名
  rownames(data) <- data[, 1]
  # 去掉第一列
  data <- data[, -1]
  # 判断数据是否需要取log
  ## 现将数据通过四分位数查看分布情况
  qx <- as.numeric(
    quantile(
      data,
      c(0., 0.25, 0.5, 0.75, 0.99, 1.0),
      na.rm = TRUE
    )
  )
  ## 给出判断条件
  ### 如果true则需要log，否则不需要log
  is_logfc <- (
    qx[5] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2
    )
  
  ## 开始判断
  if (is_logfc) {
    ### 将count数小于0的格子赋值 NaN，因为测序不可能小于0
    data[which(data <= 0)] <- NaN
    ### 然后整体取log2
    data <- log2(data)
    print("取对数操作完成")
  } else {
    print("数据不需要取对数")
  }
  
  return(data)
}

















