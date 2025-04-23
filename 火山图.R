library(ggplot2)
library(ggrepel)
library(dplyr)
library(grid)

# 读取差异分析结果数据框
data_frame <- readRDS("./output/alldiff.rds")

# 删除基因名中包含 "Rik" 或 "LOC" 的行，这些基因通常是未知或假基因
data_frame <- data_frame[!grepl("Rik|LOC", rownames(data_frame)), ]

# 设置logFC和adj.P.Val的阈值，用于筛选显著基因
logFC_threshold <- 2
adjPVal_threshold <- 0.05

# 为数据框添加新的列“Significance”，用于标记基因的调控状态
data_frame$Significance <- "Not Significant"
data_frame$Significance[data_frame$logFC > logFC_threshold & data_frame$adj.P.Val < adjPVal_threshold] <- "Upregulated"
data_frame$Significance[data_frame$logFC < -logFC_threshold & data_frame$adj.P.Val < adjPVal_threshold] <- "Downregulated"

# 筛选出上调和下调基因，并选择前25个最显著的基因
top_upregulated <- data_frame %>%
  filter(Significance == "Upregulated") %>%
  arrange(desc(logFC)) %>%
  slice(1:25)

top_downregulated <- data_frame %>%
  filter(Significance == "Downregulated") %>%
  arrange(logFC) %>%
  slice(1:25)

# 合并两个数据集，用于标注火山图中的上调和下调基因
top_genes <- bind_rows(top_upregulated, top_downregulated)

# 筛选Plin2, Plin3, Plin4基因并加入到标注列表中
plin_genes <- data_frame[rownames(data_frame) %in% c("Plin2", "Plin3", "Plin4"), ]
top_genes <- bind_rows(top_genes, plin_genes)

# 将基因名转换为大写，并作为新列“Gene_symbol”
top_genes$Gene_symbol <- (rownames(top_genes))

# 绘制火山图
volcano_plot <- ggplot(data_frame, aes(x = logFC, y = -log10(adj.P.Val), color = -log10(adj.P.Val))) +
  
  # 绘制散点图，设置透明度和点的大小
  geom_point(alpha = 0.8, size = 1.5) +
  
  # 设置颜色渐变从蓝色到红色
  scale_color_gradient(low = "#2C7BB6", high = "#D73027") +
  
  # 添加垂直阈值线（LogFC阈值）
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), 
             linetype = "dashed", color = "black") +
  
  # 添加水平阈值线（adjusted p-value 阈值）
  geom_hline(yintercept = -log10(adjPVal_threshold), 
             linetype = "dashed", color = "black") +
  
  # 设置坐标轴标签
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  
  # 使用简洁主题
  theme_minimal() +
  
  # 定制化图表样式
  theme(
    axis.line = element_line(color = "black", size = 0.5),  # 坐标轴线样式
    plot.title = element_text(hjust = 0.56, size = 16, face = "bold"),  # 标题居中并加粗
    axis.title = element_text(size = 12, hjust = 0.56),  # 坐标轴标题居中
    axis.text = element_text(size = 10),  # 坐标轴文字大小
    legend.title = element_text(size = 12),  # 图例标题样式
    legend.text = element_text(size = 10),  # 图例文本样式
    legend.key = element_rect(color = "white", fill = "white", size = 0.5),  # 图例背景样式
    legend.background = element_rect(color = "black", size = 1),  # 图例外框样式
    legend.key.size = unit(0.8, "cm"),  # 图例键的大小
    legend.spacing.y = unit(0.3, "cm"),  # 图例项之间的垂直间距
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank()  # 去除次网格线
  ) +
  
  # 设置图表标题
  ggtitle("Volcano Plot with GSE75206") +
  
  # 添加标签（标注显著基因），避免标签重叠
  geom_label_repel(data = top_genes, aes(label = Gene_symbol, fill = logFC), 
                   color = "black", fontface = "bold",  # 设置标签的样式
                   box.padding = 0.5, point.padding = 0.2, 
                   segment.color = "gray", label.size = 0.5, 
                   label.r = 0.5, alpha = 0.8, 
                   max.overlaps = 50) +
  
  # 设置标签的填充色渐变
  scale_fill_gradient2(low = "#56B4E9", mid = "#E69F00", high = "#D73027", midpoint = 0) +
  
  # 设置颜色条（colorbar）样式
  guides(
    color = guide_colorbar(title = "-Log10(adj.P.Val)", direction = "vertical", 
                           barwidth = 0.5, barheight = 10, label.position = "right"),
    
    # 设置LogFC颜色条（fill）样式
    fill = guide_colorbar(title = "LogFC", direction = "horizontal",
                          barwidth = 3.7, barheight = 0.5,
                          label.position = "bottom")
  )


# 显示火山图
print(volcano_plot)

# 保存火山图为PDF格式，设置输出路径和图像大小
ggsave("./output/火山图/volcano_plot.pdf", plot = volcano_plot, width = 10, height = 10, dpi = 600)

# 保存火山图为TIFF格式，设置输出路径、图像大小、分辨率和压缩方式
ggsave("./output/火山图/volcano_plot.tiff", plot = volcano_plot, width = 10, height = 10, dpi = 600, units = "in", compression = "lzw")


#### 简洁火山图代码 ####
# 加载必要的库
library(ggplot2)
library(ggrepel)
library(dplyr)

# 创建火山图，基因标注为斜体
volcano_plot <- ggplot(data_frame, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  
  # 绘制散点图
  geom_point(alpha = 0.8, size = 1.5) +
  
  # 设置颜色映射
  scale_color_manual(values = c("Upregulated" = "#D73027", 
                                "Downregulated" = "#4575B4", 
                                "Not Significant" = "grey")) +
  
  # 添加基因标注，使用红色字体和斜体
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)), 
                  color = "black", fontface = "italic", size = 3, 
                  box.padding = 0.5, point.padding = 0.2, 
                  segment.color = "gray", max.overlaps = 50) +
  
  # 添加阈值线
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), 
             linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(adjPVal_threshold), 
             linetype = "dashed", color = "black") +
  
  # 设置坐标轴标签
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value", 
       title = "Volcano Plot with Italicized Labels") +
  
  # 简化主题，去除图例，并添加外框
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.grid = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "none"  # 去除图例
  )

# 打印火山图
print(volcano_plot)

# 保存火山图
ggsave("volcano_plot_italic_labels.png", plot = volcano_plot, 
       width = 5, height = 5, units = "in", dpi = 300)


# 保存火山图，设置图像大小为 5x5 英寸
ggsave("volcano_plot_simple_with_border.pdf", plot = volcano_plot, 
       width = 5, height = 5, units = "in", dpi = 600)
