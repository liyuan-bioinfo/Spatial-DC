library(dplyr)
library(ggplot2)
library(RColorBrewer)

library(clusterProfiler)
# library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(limma)
library(pheatmap)
library(ggpubr)

# Basic treat
 {
      rm(list=ls())
      setwd("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/MouseKPC2023/02_proteomic_profiles/01_data")
      data_df = read.csv("SpatialDC_reconstruct_norm_filter.csv", header=T, row.names=1, check.names=FALSE) %>% t() %>% as.data.frame()
      meta_df = read.csv("SpatialDC_reconstruct_norm_filter_meta.csv", header=T, row.names=1)
      var_df = read.csv("SpatialDC_reconstruct_norm_var.csv", header=T)
      row.names(var_df) = var_df$name2

      dim(meta_df) # 658 * 5
      ct_order = c("PCC","CAF",'T4','B',"NEU","DC")
      meta_df = meta_df %>% filter(cellperc > 0.01) %>% dplyr::filter(celltype %in% ct_order)#6918
      meta_df$SampleID = row.names(meta_df)
      meta_df$celltype = factor(meta_df$celltype, levels = ct_order)
      meta_df = meta_df %>% arrange(celltype)

      data_df = data_df[,meta_df$SampleID]
      dim(data_df) # 2854 * 655

      # 计算每个ct的mean value
      data_mean_df = data.frame(row.names=row.names(data_df))
      for (i in ct_order){
          temp_df = rowMeans(data_df[meta_df$celltype==i])
          data_mean_df[i] = temp_df

      }

      dep_df = read.csv("ct10_impute_v2_dep.csv",header=T)
      dep_df$pid = dep_df$feature
      

      # 保留显著并且最富集的蛋白
      # dep_df = dep_df #%>% filter(pvals_adj < 0.05)
      filter_dep_df = data.frame()
      for (i  in unique(dep_df$pid)){
          temp_dep_df = dep_df %>% filter(pid == i)
          temp_max = which(temp_dep_df$logfoldchanges == max(temp_dep_df$logfoldchanges))
          temp_dep_df = temp_dep_df[temp_max,]
          filter_dep_df = rbind(filter_dep_df, temp_dep_df)

      }

      filter_dep_df$celltype = factor(filter_dep_df$celltype, levels = ct_order)

      # save Object
      obj_list = list()
      obj_list$data_df = data_df
      obj_list$data_mean_df = data_mean_df
      obj_list$meta_df = meta_df
      obj_list$anno_df = var_df
      obj_list$ct_order = ct_order
      obj_list$dep_df = filter_dep_df

      # ct_order_color = c(
      #     "#84CEB7","#66C2A5",  "#A3DAC9", "#C2E6DB", "#E1F3ED", # B
      #     "#FC8D62", "#FCA481", 
      #     "#377EB8",#ILC
      #     "#4DAF4A",  #NK 
      #     "#8DA0CB", "#A4B3D5", # T4
      #     "#984EA3", "#AC71B5", # T8
      #     "#A6D854", "#B7DF76", "#C9E799",#DC
      #     "#FFD92F", "#E5C494", #MO and GN
      #     "#E41A1C" #Epi
      # )
      # obj_list$ct_order_color = ct_order_color        

      saveRDS(obj_list, file="MouseKPC_spatialdc_v7_ct10_v2_20240812.RDS")
  }

# 统计分析
# -------------------------------------------------------------------------
#                   II - Signficance of samples (n >= 3)                   #
# -------------------------------------------------------------------------
# Desc: ANOVA analysis of celltypes
# Input: Obj_list$impute_df; Obj_list$impute_mean_df
# Output: Obj_list$aov_ct3_df    
{
    rm(list=ls())
    Obj_list = readRDS(file="01_data/MouseKPC_spatialdc_v7_ct10_v2_20240812.RDS") 
    meta_df = Obj_list$meta_df        
    data_df = Obj_list$data_df
    data_mean_df = Obj_list$data_mean_df    
    ct_num = 6
    enrich_fc = 1
    
    # to find enriched cell-types with mean abundance
    dep_df = data.frame()
    for(i in 1:ct_num){
        temp_df = data_mean_df[,i] - data_mean_df[,1:ct_num]
        temp_pid = names(which(rowSums(temp_df>log2(enrich_fc)) == (ct_num-1)))
        temp_df2 = data.frame(pids=temp_pid)
        temp_df2$celltype = names(data_mean_df)[i]
        temp_df2$log2FC = apply(temp_df[temp_pid,-i],1,mean)        
        dep_df = rbind(dep_df,temp_df2)            
    }

    ## pvalue, with each sample
    row.names(dep_df) = dep_df$pids    
    data_df = na.omit(data_df[dep_df$pids,])
    protein_num = dim(data_df)[1] #update for enriched proteins

    Pvalue = c()
    for(i in 1:protein_num){
        pid_df = data_df[i,] %>% t() %>% as.data.frame()
        names(pid_df) = "pid"
        pid_df$SampleId = row.names(pid_df)
        pid_df$celltype = meta_df$celltype
        model = aov(data=pid_df,pid~celltype)
        # temp_pvalue = summary(model)[[1]]$`Pr(>F)`[1]
        sig_df = TukeyHSD(model)$celltype #which ct is enriched
        enriched_ct = dep_df[i,"celltype"]
        enriched_ct_sig_df = sig_df[grep(row.names(sig_df),pattern=paste0("-",enriched_ct,"$|","^",enriched_ct,"-")),] #whethor this ct is sig. with other cts
        temp_pvalue = RecordTest::fisher.method(as.numeric(enriched_ct_sig_df[,"p adj"]))$"p.value"[[1]]  
        # temp_pvalue = max(enriched_ct_sig_df[,"p adj"])                                    
        Pvalue = c(Pvalue,temp_pvalue)            
    }        

        dep_df$pvalue = Pvalue
        dep_df$fdr = p.adjust(Pvalue,method = "BH")

        Obj_list$dep_df_fc1 = dep_df
        saveRDS(Obj_list, file="01_data/MouseKPC_spatialdc_v7_ct10_v2_20240812.RDS")    
}


