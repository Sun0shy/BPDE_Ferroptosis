#### 示例代码 ##################################################################

# 导入pheatmap库用于绘制热图
library(pheatmap)

# 设置随机种子，确保结果可重复
set.seed(123)

# 创建一个20行10列的随机数据矩阵，其中包含200个标准正态分布的随机数
test = matrix(rnorm(200), 20, 10)

# 对前10行和奇数列的数据加上3，模拟部分数据具有更高的值
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3  

# 对第11到20行和偶数列的数据加上2，模拟部分数据具有中等值
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2  

# 对第15到20行和偶数列的数据再加上4，模拟更多不同的数据分布
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4  

# 设置数据矩阵的列名为"Test1"到"Test10"
colnames(test) = paste("Test", 1:10, sep = "")  

# 设置数据矩阵的行名为"Gene1"到"Gene20"
rownames(test) = paste("Gene", 1:20, sep = "")

# 创建列注释数据：包含细胞类型(CellType)和时间(Time)两列
annotation_col = data.frame(
  CellType = factor(rep(c("CT1", "CT2"), 5)),  # 交替设置前5列为CT1，后5列为CT2
  Time = 1:5  # 设置列的时间标签为1到5
)
rownames(annotation_col) = paste("Test", 1:10, sep = "")  # 设置列注释行名，匹配数据矩阵的列名

# 创建行注释数据：包含基因类别(GeneClass)和表达水平(ExpressionLevel)两列
annotation_row = data.frame(
  GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6))),  # 按照比例为基因分配类别
  ExpressionLevel = factor(rep(c("High", "Medium", "Low"), c(10, 6, 4)))  # 按照比例为基因分配表达水平
)
rownames(annotation_row) = paste("Gene", 1:20, sep = "")  # 设置行注释行名，匹配数据矩阵的行名

# 使用pheatmap函数绘制热图，包含列和行的注释
pheatmap(test, 
         annotation_col = annotation_col,  # 添加列注释数据
         annotation_row = annotation_row,  # 添加行注释数据
         col = colorRampPalette(c("white", "darkred"))(100),  # 设置热图的颜色渐变，从白色到深红色
         show_rownames = TRUE,  # 显示行名
         show_colnames = TRUE,  # 显示列名
         main = "Heatmap with Multiple Row Annotations",  # 设置热图标题
         annotation_rows_side = "left",  # 将行注释放置于热图的左侧
         annotation_names_row = T,  # 显示行注释的名称
         cluster_cols = T,  # 对列进行聚类
         cluster_rows = T,  # 对行进行聚类
         cutree_row = 3,  # 将行分成3个簇
         cutree_cols = 2   # 将列分成2个簇
)



#### 真实数据绘图代码 ##########################################################
# 导入所需的库
library(pheatmap)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(writexl)
library(openxlsx)

# 读取表达矩阵数据
Martix <- readRDS("./output/exprset.rds")

Martix$gene_name <- rownames(Martix)
# 将数据保存为 Excel 文件
write.xlsx(Martix, file = "./output/exprset.xlsx")

# 从Excel文件中读取选定基因的列表
gene <- read_excel("./data/热图选定通路.xlsx")

# 将geneID列拆分为多个基因ID，展开为多行
long_data <- gene %>%
  mutate(geneID = strsplit(as.character(geneID), "/")) %>%
  unnest(geneID)

# 读取差异分析结果数据
all_diff <- readRDS("./output/alldiff.rds")

# 删除基因名中包含'Rik'或'LOC'的基因
all_diff <- all_diff[!grepl("Rik|LOC", rownames(all_diff)), ]

# 筛选差异表达基因，条件是logFC的绝对值大于3，且调整后的p值小于0.05
filtered_data <- all_diff %>%
  filter(abs(logFC) > 3, adj.P.Val < 0.05)

# 将基因符号转换为大写并重命名为geneID
filtered_data$Gene_symbol <- toupper(rownames(filtered_data))  
colnames(filtered_data)[7] <- "geneID"

# 将差异分析结果与选定基因数据进行合并
merge_gene <- merge(filtered_data, long_data, by = "geneID")

# 使用distinct去除重复的geneID，保留每个geneID的第一行数据，按logFC的绝对值排序，取前70个基因
merge_gene_unique <- merge_gene %>%
  distinct(geneID, .keep_all = TRUE) %>%
  arrange(abs(logFC)) %>%
  head(70) %>%
  mutate(UpDown = ifelse(logFC > 0, "Up", "Down"))

# 按照Pathway_En列的值排序
merge_gene_unique <- merge_gene_unique %>% 
  arrange(Pathway_En)

# 修改Pathway_En列，针对特定基因进行路径名更新
merge_gene_modified <- merge_gene_unique %>%
  mutate(Pathway_En = ifelse(geneID %in% c("PLIN4", "PLIN2", "NFE2L2"), 
                             "Ferroptosis + Lipid metabolism", 
                             Pathway_En))

# 创建行注释数据框，包含基因ID、路径和表达水平
merge_row_annotation <- data.frame(merge_gene_modified$geneID, merge_gene_modified$Pathway_En, merge_gene_modified$UpDown)
rownames(merge_row_annotation) <- merge_row_annotation$merge_gene_modified.geneID
merge_row_annotation <- merge_row_annotation[, c(-1)]
colnames(merge_row_annotation) <- c("Pathway", "Expression")

# 生成列注释数据，设置“Group”列为因子类型并指定水平
merge_col_annotation <- data.frame(Group = factor(c(rep("Control", 4), rep("B[a]P", 4)), levels = c("Control", "B[a]P")))
rownames(merge_col_annotation) = colnames(Martix)

# 将Martix中的基因ID转换为大写，确保与其他数据一致
Martix$geneID <- toupper(rownames(Martix))
rownames(Martix) <- Martix$geneID

# 从merge_gene_modified中提取基因ID列，选择对应基因的数据
gene_ids_to_select <- merge_gene_modified$geneID
Martix_selected <- Martix[gene_ids_to_select, ]

# 移除选定数据的最后一列（可能是无关的数据）
Martix_selected <- Martix_selected[, -ncol(Martix)]

# 设置列注释颜色
groupcolor <- c("#2C7BB6", "#F9DB40") 
names(groupcolor) <- c("Control", "B[a]P")

# 设置行注释颜色：根据不同的路径设置颜色
Pathwaycolor <- c("#349CE2", "#F7AA1F", "#EDDE4A", "#AF55C1", "#00E64D")
names(Pathwaycolor) <- c("Ferroptosis", "Ferroptosis + Lipid metabolism", "Lipid metabolism", 
                         "Neurodegenerative diseases", "Oxidative stress")

# 设置表达水平的颜色：根据UpDown列设置颜色
Expressioncolor <- c("#3548FD", "#E60000")
names(Expressioncolor) <- c("Down", "Up")

# 将所有注释颜色合并为一个列表
ann_colors <- list(
  Group = groupcolor,
  Pathway = Pathwaycolor,
  Expression = Expressioncolor
)

# 定义首字母大写函数
capitalize <- function(x) {
  sapply(x, function(word) {
    paste(toupper(substring(word, 1, 1)), tolower(substring(word, 2)), sep = "")
  })
}

# 修改 Martix_selected 的行名
rownames(Martix_selected) <- capitalize(rownames(Martix_selected))

# 修改 merge_row_annotation 的行名和列
rownames(merge_row_annotation) <- capitalize(rownames(merge_row_annotation))

# 如果 merge_row_annotation 有 geneID 列
merge_row_annotation$geneID <- capitalize(merge_row_annotation$geneID)

# 修改 merge_gene_modified 的 geneID
merge_gene_modified$geneID <- capitalize(merge_gene_modified$geneID)

# 使用pheatmap绘制热图，添加多个行和列注释
heatmap <- pheatmap(Martix_selected, 
                    # 设置热图颜色渐变，从蓝色到红色，通过colorRampPalette指定渐变颜色
                    col = colorRampPalette(c("#2C7BB6", "lightblue", "white", "pink", "darkred"))(100),
                    
                    # 设置热图标题
                    main = "Heatmap with GSE75206",
                    
                    # 设置是否显示行名和列名
                    show_rownames = TRUE, 
                    show_colnames = FALSE,
                    
                    # 设置列注释数据，包含了列的分组信息
                    annotation_col = merge_col_annotation, 
                    
                    # 设置行注释数据，包含基因的通路信息和表达水平
                    annotation_row = merge_row_annotation,  
                    
                    # 设置注释颜色的方案
                    annotation_colors = ann_colors,
                    
                    # 设置行注释的位置为左侧
                    annotation_rows_side = "left",
                    
                    # 是否显示行注释的名称，设为TRUE表示显示
                    annotation_names_row = T,
                    
                    # 是否对列进行聚类，设为FALSE表示不进行列聚类
                    cluster_cols = F,
                    
                    # 是否对行进行聚类，设为TRUE表示对行进行聚类
                    cluster_rows = T,
                    
                    # 行聚类树高度
                    treeheight_row = 25,
                    
                    # 设置单元格的宽度和高度，单位为像素
                    cellwidth = 10,
                    cellheight =5,
                    
                    # 设置字体大小，列名、行名、显示数值的字体大小
                    fontsize_col = 12,
                    fontsize_row = 6,
                    fontsize_number = 5,
                    
                    # 是否在每个单元格中显示数值，设为TRUE表示显示
                    display_numbers = F,
                    
                    # 设置数值的颜色为黑色
                    number_color = "black",
                    
                    # 设置单元格的边框颜色为灰色
                    border_color = "NA",
                    
                    # 设置行聚类的分割数量为4
                    cutree_row = 4, 
                    
                    # 设置列聚类的分割数量为2
                    cutree_cols = 2,
                    
                    # 设置图片大小
                    width = 8,
                    height = 8
)

# 保存热图为PDF格式
ggsave("./output/复杂热图/heatmap调整.pdf", plot = heatmap, width = 8, height = 8, dpi = 600)

# 保存热图为TIFF格式
ggsave("./output/复杂热图/heatmap调整.tiff", plot = heatmap, width = 8, height = 8, dpi = 600, units = "in", compression = "lzw")


