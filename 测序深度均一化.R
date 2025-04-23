source("./scripts/图片保存到output_figure.R")

normalize_sample <- function(geo_data) {
  # 未均一化箱图
  no_one_box <- boxplot(geo_data, outline = FALSE, notch = TRUE, las = 2)
  # 保存未均一化的图
  figure_save(no_one_box)
  print("未均一化的图保存成功")
  
  # 矫正测序深度使其一致
  ## 默认使用quntile矫正差异
  geo_data <- normalizeBetweenArrays(geo_data)
  # 均一化箱图
  oned_box <- boxplot(geo_data, outline = FALSE, notch = TRUE, las = 2)
  # 保存均一化的图
  figure_save(oned_box)
  print("均一化的图保存成功")
  geo_data <- as.data.frame(geo_data)
  print("探针数据转为数据框成功")
  return(geo_data)
}