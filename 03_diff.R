rm(list = ls())


load("./output/exprSet.Rdata")


group = c(
  rep("control", 4),
  rep("BaP", 4)
)

group <- factor(group, levels = c("control", "BaP"))
# group <- factor(group, levels = c("B[a]P","control"))

res.pca <- prcomp(t(exprSet), scale = TRUE)

if(F){
  # 1.转置
  data <- t(exprSet)
  ## 2.NAN值变成0
  data[is.nan(data)] <- 0 
  ## 3.删除没有变异的基因，使用函数var
  data <-  data[,which(apply(data,2,var) != 0)]
  res.pca <- prcomp(data, scale = TRUE)
}

library(factoextra)
fviz_pca_ind(res.pca,col.ind = group, addEllipses = T)


### 1.构建比较矩阵
# model.matrix(~group) vs model.matrix(~group + 0)
# 1. 差异点
# ~group：包含截距，默认以对照组为基准。结果中的系数（logFC）表示实验组相对于对照组的差异。
# ~group + 0：不包含截距，每组有独立系数。结果中的系数表示各组的绝对表达水平，需手动计算组间差异。
# 2. 适用场景
# ~group：有明确对照组，研究组间差异（如对照 vs 实验组）。
# ~group + 0：无对照组，关注各组的独立表达水平。
# 3. 推荐
# 有对照组实验（如苯并芘染毒 vs 玉米油对照）：用 ~group，直接获得组间差异，结果解读简单。
# 无对照组或需单独评估组水平：用 ~group + 0。
design <- model.matrix(~group)
### 比较矩阵命名
colnames(design) <- levels(group)
design

library(limma)
### 2.线性模型拟合
contrast.matrix <- makeContrasts(BaP - control, levels = group)
fit <- lmFit(exprSet,design)
fit <- contrasts.fit(fit, contrast.matrix)
### 3.贝叶斯检验
fit2 <- eBayes(fit)
### 4.输出差异分析结果,其中coef的数目不能操过design的列数
### 此处的2代表的是design中第二列和第一列的比较
allDiff=topTable(fit2,adjust='fdr',coef=1,number=Inf) 
### 这个数据很重要需要保存一下
allDiff["Plin4", ]
# save(allDiff,file = "output/allDiff.Rdata")

saveRDS(allDiff,"./output/alldiff.rds")





























