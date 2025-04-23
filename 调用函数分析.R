# 加载探针数据
probe_data = read.table(
  './data/probe.txt',
  comment.char = '!',
  stringsAsFactors = F,
  header = T
  )

# 判断数据是否log，如果没有，数据将被自动log
source('scripts/探针数据log判断.R')
loged_data = is_probe_log(probe_data)

# 将测序数据是否均一化，保存均一化前后的箱图，图片名字为变量名
source('scripts/图片保存到output_figure.R')

# 并将数据转换为数据框
source('scripts/测序深度均一化.R')
oned_data = normalize_sample(loged_data)

# 加载探针基因名对照数据,类型为数据框
probe2symbol <- data.table::fread('./data/ann.txt',skip = 17, data.table = F)
# 取出探针id列和symbol列
probe2symbol <- probe2symbol[,c(1,7)]

# 将oned_data的行名提取成一列一遍能merge symbol
library(tibble)
library(dplyr)
oned_data_symbol <- oned_data %>% 
  # 行名变第一列，以便能merge symobl列
  rownames_to_column(colnames(probe2symbol)[1]) %>%
  # 按照id合并
  inner_join(probe2symbol,.,by = 'ID') %>%
  # 合并后去掉id列
  select(-ID) %>%
  # 重命名
  rename(symbol = GENE_SYMBOL) %>%
  # 去掉没有对应基因名的行
  filter(symbol != '') %>%
  # 除了第一列外，每一行的其他列做均值计算，用来排除重复的基因名
  mutate(rowMeans = rowMeans(.[,-1])) %>%
  # 排序这些行，降序
  arrange(desc(rowMeans)) %>%
  # 去重symbol这一列，只保留最大的第一个，其他列也保留
  distinct(symbol,.keep_all = T) %>%
  # 去掉rowmeans列
  select(-rowMeans) %>%
  # 第一列转为行名
  column_to_rownames('symbol')


# 保存数据
diff_data <- oned_data_symbol
save(diff_data,file = 'output/outdata/diff_data.Rdata')


# 加载差异分析数据
# load('output/outdata/diff_data.Rdata')

# 调整样本顺序con con con con B(a)P B(a)P B(a)P B(a)P
# GSM1945106	Control_1
# GSM1945107	BaP_1
# GSM1945108	Control_2
# GSM1945109	BaP_2
# GSM1945110	BaP_3
# GSM1945111	Control_3
# GSM1945112	BaP_4
# GSM1945113	Control_4

diff_data <- diff_data[,c(1,3,6,8,2,4,5,7)]
# "GSM1945106" "GSM1945108" "GSM1945111" "GSM1945113"   con
# "GSM1945107" "GSM1945109" "GSM1945110" "GSM1945112"   B(a)P
# 确定分组
group <- rep(c('B(a)P','con'),each=4)
# levelsl里面把对照组放在前面
group <- factor(group,levels = c('B(a)P','con'))

# 主成分分析pca
source('scripts/样本pca分析——看测序结果好坏.R')
pca_result <- sample_pca(diff_data = diff_data,group = group,title = 'sample_pca') 
# 保存图片
figure_save(pca_result)

# 差异分析
source('scripts/差异分析.R')
allDiff <- gene_diff_analyse(diff_data,group)
# 保存r文件
save(allDiff,file = 'output/outdata/allDiff.Rdata')
# 保存为txt
write.table(
  allDiff,
  file = 'output/outdata/差异分析.txt',
  sep = '\t',
  row.names = F
  )



# 定义差异基因标准，P.Value<0.05 and logfc>0.5
library(dplyr)
diffgene <- allDiff %>% 
  filter(.$P.Value < 0.05) %>% 
  filter(abs(logFC) > 0.5)

# 保存差异基因数据
save(diffgene,file = 'output/outdata/diffgene.Rdata')
# 保存为txt
write.table(
  diffgene,
  file = 'output/outdata/所有差异表达基因.txt',
  sep = '\t',
  row.names = F
)
 
# 作图
## 行列转置
expre_set <- as.data.frame(t(diff_data))
## 添加分组信息
dd <- cbind(group=group,expre_set)
# 调整出图顺序
dd$group <- factor(dd$group,ordered = TRUE,levels=c('con','B(a)P'))
## 画图 绘制某个基因差异表达对比

# 指定比较的组，可以放多个c()
my_comparisons <- list(c('B(a)P','con'))
# 绘制基因对照箱图
source('scripts/差异基因box图对比.R')
plin4_box <- diff_boxplot(dd,'Plin4',my_comparisons)
# 保存图片
source('scripts/图片保存到output_figure.R')
figure_save(plin4_diff)



# # 绘制热图
# source('scripts/热图绘制.R')
# hot <- draw_heatmap(allDiff = allDiff,0.5,0.6,group_heat = group_heat)
# 绘制热图
## 处理数据 
heatdata <- diff_data[rownames(diffgene),]
## 制作一个分组信息用于注释
group_heat <- rep(c('con','B(a)P'),each=4)
annotation_col <- data.frame(group_heat)
rownames(annotation_col) <- colnames(heatdata)

## 删除无效基因
out_rik <- grep('Rik',diffgene$GENE_ID,ignore.case = FALSE)
out_gm <- grep('Gm',diffgene$GENE_ID,ignore.case = FALSE)
out_loc <- grep('LOC',diffgene$GENE_ID,ignore.case = FALSE)
filtered_diffgene <- diffgene[-c(out_rik,out_gm,out_loc),]

## 处理数据 
heatdata <- diff_data[rownames(filtered_diffgene),]
## 制作一个分组信息用于注释
group_heat <- rep(c('con','B(a)P'),each=4)
annotation_col <- data.frame(group_heat)
rownames(annotation_col) <- colnames(heatdata)


#如果注释出界，可以通过调整格子比例和字体修正
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

# 保存热图
source('scripts/图片保存到output_figure.R')
figure_save(hot_plot)


# 绘制火山图
library(ggplot2)
library(ggrepel)
library(dplyr)
## allDiff增加基因列，火山图数据至少三列logfc，p.value,gene
vol_data <- allDiff
vol_data$gene <- rownames(vol_data)
# 火山图竖线所在位置
logFCfilter <- 0.5
# 极值倍数指定
logFCcolor <- 0.8
# 标记上下调
index <- vol_data$P.Value < 0.05 & abs(vol_data$logFC) > logFCfilter
vol_data$group <- 0
vol_data$group[index & vol_data$logFC > 0] <- 1
vol_data$group[index & vol_data$logFC < 0] <- -1
vol_data$group <- factor(
  vol_data$group, 
  levels = c(1,0,-1),
  labels = c('Up','NS','Down')
  )

# 正式画图
vol_plot <- ggplot(data=vol_data, aes(x=logFC, y =-log10(P.Value),color=group)) +
  geom_point(alpha=0.8, size=1.2)+
  scale_color_manual(values = c("red", "grey50", "blue4"))+
  labs(x="log2 (fold change)",y="-log10 (adj.P.Val)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(-logFCfilter,logFCfilter),lty=4,lwd=0.6,alpha=0.8)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(legend.position="top")+
  geom_point(data=subset(vol_data, abs(vol_data$logFC) >= logFCcolor & vol_data$P.Value <0.05),alpha=0.8, size=3,col="green4")+
  geom_text_repel(data=subset(vol_data, abs(vol_data$logFC) >= logFCcolor & vol_data$P.Value <0.05),
                  aes(label=gene),col="black",alpha = 0.8)

# 保存火山图
source('scripts/图片保存到output_figure.R')
figure_save(vol_plot)











































