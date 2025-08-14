library(limma)


gene_diff_analyse <- function(data,group) {
  ## 1.构建比较矩阵
  design <- model.matrix(~group)
  ## 2.比较矩阵命名
  colnames(design) <- levels(group)
  ## 3.线性模型拟合
  fit <- lmFit(data,design)
  ## 4.贝叶斯检验
  fit2 <- eBayes(fit)
  ## 5.输出差异分析结果，其中coef的数目不能超过design的列数
  ### coef=2代表design中第二列和第一列比较
  allDiff <- topTable(fit2,adjust='fdr',coef = 2,number = Inf)
  ## 6.加入gene列
  allDiff <- data.frame(GENE_ID=rownames(allDiff),allDiff)
  return(allDiff)
}




















