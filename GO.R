#1、加载包
library(AnnotationDbi)
library(org.Hs.eg.db)#基因注释包
library(clusterProfiler)#富集包
library(dplyr)
library(ggplot2)#画图包


# 读取差异基因
## 读取logFC倍数为2的全部差异基因
all_diff <- readRDS("./output/alldiff.rds")
logover2_diff <- all_diff[abs(all_diff$logFC) > 2 & all_diff$adj.P.Val < 0.05,]
### 将所有基因大写
logover2_diff$gene_symbol <- toupper(rownames(logover2_diff))

## 读取logFC倍数为2的铁死亡差异基因数据
### 读取所有的铁死亡基因
all_ferr <- read.table("./output/ferroptosis.txt")
### 交集所有差异铁死亡基因
logover2_diff_gene_symbol <- as.data.frame(logover2_diff$gene_symbol)
index <- intersect(all_ferr$V1, logover2_diff_gene_symbol$`logover2_diff$gene_symbol`)
logover2_diff_ferr <- logover2_diff[logover2_diff$gene_symbol %in% index, ]

# 报错所有loglogFC倍数为2的差异基因和铁死亡差异基因
saveRDS(logover2_diff, "./output/rdata/logover2_diff")
saveRDS(logover2_diff_ferr, "./output/rdata/logover2_diff_ferr")


# 基因symbol --> id
diff_id <- bitr(logover2_diff$gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
diff_id_ferr <- bitr(logover2_diff_ferr$gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# 获取所有的ENTREZID
diff_ENTREZID <- diff_id$ENTREZID
diff_ferr_ENTREZID <- diff_id_ferr$ENTREZID

# 要绘制铁死亡的就把下面的diff_ENTREZID改为diff_ferr_ENTREZID

ego_ALL <- enrichGO(gene = diff_ENTREZID,#我们上面定义了
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",#富集的GO类型
                    pAdjustMethod = "BH",#这个不用管，一般都用的BH
                    minGSSize = 1,
                    pvalueCutoff = 0.01,#P值可以取0.05
                    qvalueCutoff = 0.05,
                    readable = TRUE)
ego_CC <- enrichGO(gene = diff_ENTREZID,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_BP <- enrichGO(gene = diff_ENTREZID,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_MF <- enrichGO(gene = diff_ENTREZID,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

#4、将结果保存到当前路径
ego_ALL <- as.data.frame(ego_ALL)
ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- as.data.frame(ego_CC)
ego_result_MF <- as.data.frame(ego_MF)
ego <- rbind(ego_result_BP,ego_result_CC,ego_result_MF)#或者这样也能得到ego_ALL一样的结果
write.csv(ego_ALL,file = "./output/rdata/ego_ALL.csv",row.names = T)
write.csv(ego_result_BP,file = "./output/rdata/ego_result_BP.csv",row.names = T)
write.csv(ego_result_CC,file = "./output/rdata/ego_result_CC.csv",row.names = T)
write.csv(ego_result_MF,file = "./output/rdata/ego_result_MF.csv",row.names = T)
write.csv(ego,file = "./output/rdata/ego.csv",row.names = T)

#5、但是有时候我们富集到的通路太多了，不方便绘制图片，可以选择部分绘制，这时候跳过第4步，来到第5步
display_number = c(22, 22, 22)#这三个数字分别代表选取的BP、CC、MF的通路条数，这个自己设置就行了
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]

##将以上我们摘取的部分通路重新组合成数据框
go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                         Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number[1]), 
                rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process", "cellular component","molecular function" )))

##通路的名字太长了，我选取了通路的前五个单词作为通路的名字
for(i in 1:nrow(go_enrich_df)){
  description_splite=strsplit(go_enrich_df$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
  go_enrich_df$Description[i]=description_collapse
  go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
}

##开始绘制GO柱状图
###横着的柱状图
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")#设定颜色

ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
  scale_fill_manual(values = COLS) + ###颜色
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()

###竖着的柱状图 
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  theme_bw() + 
  xlab("GO term") + 
  ylab("Num of Genes") + 
  labs(title = "The Most Enriched GO Terms")+ 
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle是坐标轴字体倾斜的角度，可以自己设置






















