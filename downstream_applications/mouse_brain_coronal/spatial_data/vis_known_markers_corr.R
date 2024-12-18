library(dplyr)

setwd("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/MouseBrain2022/01_cell_perc")
rm(list=ls())

# 评价细胞丰度与空间区域划分的能力。
data_df = read.csv("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/MouseBrain2022/01_cell_perc/analysis/known_markers_high_corr.csv",header=T,row.names=1)

corr_df = cor(t(data_df), method="pearson")

# set break range for pheatmap 
range_all <- range(c(-1, 1))    
my_palette <- colorRampPalette(c("navy", "white","firebrick3"))(n=100)
breaks = seq(range_all[1], range_all[2], length.out = 101)


p1=pheatmap::pheatmap(corr_df, color=my_palette, breaks=breaks,display_numbers=T)

pdf("figures/panel_marker_corr_spatialDC.pdf",height=5,width=5.5)
print(p1)
dev.off()

write.csv(corr_df[p1$tree_row$labels[p1$tree_row$order],p1$tree_col$labels[p1$tree_col$order]],"figures/panel_marker_corr_spatialDC.csv")
