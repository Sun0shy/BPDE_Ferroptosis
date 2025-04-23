# 1. 读取矩阵数据 ####
exprSet <- read.table(
  "./data/GSE75206_series_matrix.txt.gz",
  comment.char = "!",
  stringsAsFactors = F,
  header = T)

# 2. 第一列变为行名 ####
rownames(exprSet) <- exprSet[, 1]

# 3. 去掉第一行 ####
exprSet <- exprSet[, -1]

# 4. 调整样本顺序
# "GSM1945106" "GSM1945107" "GSM1945108" "GSM1945109"
#   control        B[a]P      control       B[a]P
# "GSM1945110" "GSM1945111" "GSM1945112" "GSM1945113"
#    B[a]P        control      B[a]P        control
exprSet <- exprSet[,c(1,3,6,8,2,4,5,7)]