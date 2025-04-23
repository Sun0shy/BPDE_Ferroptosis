library(export)


figure_save <- function(data) {
  # 将变量名称变为字符串
  name =as.character(deparse(substitute(data)))
  
  # 不同后缀
  ## 导成PPT可编辑的格式
  pptx = paste0('output/figure/',name,'.pptx')
  ## 导成pdf
  pdf = paste0('output/figure/',name,'.pdf')
  ## 导成tiff
  tif = paste0('output/figure/',name,'.tif')
  ## 导成png
  png = paste0('output/figure/',name,'.png')
  ## 导成eps
  eps = paste0('output/figure/',name,'.eps')
  
  graph2ppt(file=pptx)
  graph2pdf(file=pdf)
  graph2tif(file=tif)
  graph2png(file=png)
  graph2eps(file=eps)
  
  print('图片全部保存成功到output/figure中')
}



