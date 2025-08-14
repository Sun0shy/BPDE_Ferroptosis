exprSet <- readRDS("./output/exprset.rds")
alldiff <- readRDS("./output/alldiff.rds")


library(dplyr)
library(tibble)

# 设置阈值
diffgene <- alldiff %>% 
  rownames_to_column() %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) > 2) %>% 
  column_to_rownames()

# 删除带有Rik和LOC的行名
diffgene <- diffgene[!grepl("Rik|LOC", rownames(diffgene)), ]
exprSet <- exprSet[!grepl("Rik|LOC", rownames(exprSet)), ]
alldiff <- alldiff[!grepl("Rik|LOC", rownames(alldiff)), ]

library(writexl)
alldiff$Gene <- rownames(alldiff)
alldiff <- alldiff[, c(7,1:6)]
GENE_alldiff <- toupper(alldiff$Gene)
alldiff$Gene <- GENE_alldiff 
write_xlsx(alldiff, "./output/GENE_alldiff.xls")
# 从表达矩阵提取diffgene绘制热图
heatdata <- exprSet[rownames(diffgene),]
# 设置分组
group <- c(rep("Control", 4), rep("B[a]P", 4))
annotation_col <- data.frame(group)
rownames(annotation_col) <- colnames(heatdata)

library(pheatmap)
library(viridisLite)

pheatmap(heatdata,
         cluster_cols = T,
         cluster_rows = T,
         annotation_col = annotation_col,
         annotation_legend = T,
         show_rownames = F,
         scale = "row",
         color = viridis(10, alpha = 1, begin = 0.5, end = 1, direction = 1),
         cellwidth = 60,
         cellheight = 0.4,
         fontsize = 10)

# 重新选择阈值
diffgene2 <- alldiff %>% 
  rownames_to_column() %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) >3) %>% 
  arrange(abs(logFC)) %>% 
  column_to_rownames()
# 从表达矩阵提取diffgene绘制热图
heatdata2 <- exprSet[rownames(diffgene2),]

## 生成图片大小
options(repr.plot.width = 20, repr.plot.height = 100)
pheatmap(heatdata2,
         cluster_cols = F,
         cluster_rows = T,
         annotation_col = annotation_col,
         annotation_legend = T,
         show_rownames = F,
         scale = "row",
         cellwidth = 30,
         cellheight = 0.1,
         fontsize = 10)


















#载入ComplexHeatmap和circlize包；
library(ComplexHeatmap)
library(circlize)


df <- heatdata2
#对数据进行归一化；
#由于scale函数默认对列进行归一化，因此这里做了两次转置；
df_scaled <- t(scale(t(df)))
#查看归一化后的数据前6行；
head(df_scaled)

#初步尝试绘制热图；
Heatmap(df_scaled,row_names_gp = gpar(fontsize = 6)
        ,column_names_gp = gpar(fontsize = 8),
        name = "Exp")

range(df_scaled)
#然后，根据数据范围建立自定义颜色映射关系；
#green-red;
col_fun = colorRamp2(c(-1.6,0.6,2.8), c("greenyellow","white", "red"))
#green-purple;
col_fun = colorRamp2(c(-1.6,0.6,2.8), c("greenyellow","white", "purple"))
#purple-orange;
col_fun = colorRamp2(c(-1.6,0.6,2.8), c("purple","white", "orange"))
#使用自定义渐变色绘制热图；
Heatmap(df_scaled,row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8),
        col = col_fun,
        name = "Exp")


#根据聚类树添加gap;
Heatmap(df_scaled,row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8),
        col = col_fun,
        column_title = NULL,
        show_row_dend = FALSE,
        name = "Exp")


#生成分组颜色条注释；
class = anno_block(gp = gpar(fill = c("white","red"),
                             col="white"),height = unit(5, "mm"),
                   labels = c("TESA", "TESB"),
                   labels_gp = gpar(col = "white", fontsize = 8,fontface="bold"))
group= HeatmapAnnotation(group=class)

#为热图添加分组颜色条；
Heatmap(df_scaled,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8),
        col = col_fun,
        column_title = NULL,
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        top_annotation =group,
        name = "Exp")

# 假设要把表格中的前7个基因标记在热图上；
mark <- rownames(df)[1:7]
mark
#[1] "ERMP1" "SEC31A" "LOC101795171" "ATOX1" "CORO2B" "ZEB1"
#[7] "PGAM1"
#anno_mark()至少需要两个参数，其中at是原始数据矩阵的索引，标签是对应的文本。
lab = rowAnnotation(ano = anno_mark(at = 1:7,
                                    labels = mark,
                                    labels_gp = gpar(fontsize = 8)))

#在热图行方向添加标记；
Heatmap(df_scaled,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8),
        col = col_fun,
        column_split = 4,
        column_title = NULL,
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        show_row_names = FALSE,
        top_annotation =group,
        right_annotation = lab,
        row_names_side = "left",
        name = "Exp")


#当然，也可以指定基因集后再获取对应的原矩阵的位置；
genelist <- c("Gfap")
index <- which(rownames(df) %in% genelist)

#得到对应的文本标签；
labs <- rownames(df)[index]

#使用labels_gp调整字体大小；
lab2 = rowAnnotation(foo = anno_mark(at = index,
                                     labels = labs,
                                     labels_gp = gpar(fontsize = 8),
                                     lines_gp = gpar()))

#添加白色格子线，再看下热图的效果；
Heatmap(df_scaled,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8),
        col = col_fun,
        rect_gp = gpar(col="white"),
        column_split = 2,
        column_title = NULL,
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        top_annotation =group,
        right_annotation = lab2,
        show_row_names = FALSE,
        width = 12,
        name = "Exp")
