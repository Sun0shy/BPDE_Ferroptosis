library(ggpubr)

diff_boxplot <- function(data,gene,comparisons) {
  
  
  diff_box <- ggboxplot(
    data, 
    x = "group", 
    y = gene,
    color = "group", 
    palette = c("#00AFBB", "#E7B800"),
    add = "jitter",
    # grids = TRUE,
    # bgcolor('grey')
    ) + 
    # 加t检验
    stat_compare_means(
      comparisons = comparisons, 
      method = "t.test")
  return(diff_box)
}