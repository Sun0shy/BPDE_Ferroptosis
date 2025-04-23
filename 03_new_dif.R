library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)


load("./output/exprSet.Rdata")
rt <- exprSet
rt <- as.matrix(rt)
exp <- rt
dimnames <- list(rownames(exp), colnames(exp))
rt <- matrix(as.numeric(exp), nrow = nrow(exp), dimnames = dimnames)



rt1 <- normalizeBetweenArrays(rt)

data <- rt1

afcon <- 4
conData <- data[,as.vector(colnames(data)[1:afcon])]
aftrate <- afcon + 1
treatData <- data[,as.vector(colnames(data)[aftrate:ncol(data)])]
rt <- cbind(conData, treatData)
conNum <- ncol(conData)
treatNum <- ncol(treatData)

Type <- c(rep("con", conNum), rep("treat", treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con", "treat")
fit <- lmFit(rt, design)
con.matrix <- makeContrasts(treat-con, levels = design)
fit2 <- contrasts.fit(fit, con.matrix)
fit2 <- eBayes(fit2)

Diff <- topTable(fit2,adjust = "fdr", number = length(rownames(data)))

DIFFOUT <- rbind(id = colnames(Diff), Diff)
write.table(DIFFOUT, file = "./output/DIFF_all.xls", sep = "\t", quote = F, col.names = F)
s
Diff <- Diff[order(as.numeric(as.vector(Diff$logFC))),]
diffGene <- as.vector(rownames(Diff))
diffLength <- length(diffGene)
afGene = c()
if (diffLength > (100) ) {
  afGene <- diffGene[c(1:50, (diffLength - 50 + 1) : diffLength)]
} else {
  afGene <- diffGene
}
afExp <- rt[afGene,]

Type <- c(rep("con", conNum), rep("treat", treatNum))
names(Type) <- colnames(rt)
Type <- as.data.frame(Type)

anncolor <- list(Type = c(con = pal_npg()(1), treat = pal_npg()(2)[2]))

pheatmap(afExp,
         annotation = Type,
         color = colorRampPalette(c(pal_npg()(2)[2], "white", pal_npg()(1)))(50),
         cluster_cols = F,
         show_colnames = F,
         scale = "row",
         fontsize = 8,
         fontsize_col = 8,
         fontsize_row = 6,
         annotation_colors = anncolor
)

