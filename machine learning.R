score <- read.csv('output/outdata/得分表格.csv')
score_genelist <- as.data.frame(score[,c(1,4)])
colnames(score_genelist) <- c('score','gene')

socre_gene_log <- allDiff[allDiff[,1] %in% score_genelist[,2],c(1,2)]
score_gene_log <- as.data.frame(socre_gene_log)

colnames(score_gene_log) <- c('gene','logFC')

score_log_merged <- merge(score_gene_log,score_genelist,by='gene')

write.table(
  score_log_merged,
  file = 'output/outdata/机器学习.txt',
  sep = '\t',
  row.names = F
)
