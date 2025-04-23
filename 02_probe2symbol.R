
ex <- exprSet

qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)


# 判断是否需要log
if (LogC) { 
  ex[which(ex <= 0)] <- NaN
  ## 取log2
  exprSet <- log2(ex)
  print("log2 transform finished")
}else{
  print("log2 transform not needed")  
}

library(limma) 
boxplot(exprSet,outline=FALSE, notch=T, las=2)
### 该函数默认使用quntile 矫正差异 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T, las=2)
## 这步把矩阵转换为数据框很重要
exprSet <- as.data.frame(exprSet)


platformMap <- data.table::fread(
  "./data/GPL10787-9758.txt",
  data.table = F)
platformMap <- platformMap[, c(1,7)]
platformMap <- platformMap[platformMap$GENE_SYMBOL != "",]
probe2symbol_df <- platformMap
colnames(probe2symbol_df) <- c("ID", "symbol")


library(dplyr)
library(tibble)
exprSet <- exprSet %>% 
  ## 行名转列名,因为只有变成数据框的列,才可以用inner_join
  rownames_to_column("ID") %>% 
  ## 合并探针的信息
  inner_join(probe2symbol_df,by="ID") %>% 
  ## 去掉多余信息
  select(-ID) %>%  
  ## 重新排列
  select(symbol,everything()) %>%  
  ## rowMeans求出行的平均数(这边的.代表上面传入的数据)
  ## .[,-1]表示去掉出入数据的第一列，然后求行的平均值
  mutate(rowMean =rowMeans(.[,-1])) %>% 
  ## 把表达量的平均值按从大到小排序
  arrange(desc(rowMean)) %>% 
  ## 去重，symbol留下第一个
  distinct(symbol,.keep_all = T) %>% 
  ## 反向选择去除rowMean这一列
  select(-rowMean) %>% 
  ## 列名转行名
  column_to_rownames("symbol")


### 保存数据
save(exprSet,file = "./output/exprSet.Rdata")
