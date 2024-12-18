library(dplyr)
library(ggplot2)
library(RColorBrewer)

library(clusterProfiler)
# library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(limma)
library(pheatmap)
library(ggpubr)

    # Part 1, PPI plot
  # 展示分泌蛋白与膜蛋白组成的PPI网络特征
  {
      rm(list=ls())        
      obj_list = readRDS(file = "01_data/MouseKPC_spatialdc_v7_ct10_v2_20240812.RDS")
      row.names(obj_list$anno_df) = obj_list$anno_df$pid
      anno_df = obj_list$anno_df
      head(anno_df)

      loci_order = c("Sec","PM")
                  
      # node_df 
      node_df = obj_list$identified_SP_df
      node_df = node_df[!is.na(node_df$loci),]

      node_df = node_df %>% dplyr::select(pid,loci) %>% unique()

      dim(node_df) # 763 * 2
      node_df$genename = anno_df[node_df$pid,"gene"]

      # lind_df => all sec - pm pairs
      from_df = node_df %>% filter(loci == "Sec") %>% unique()# %>% filter(pid %in% ppi_db$ligand_pid)
      to_df = node_df %>% filter(loci == "PM") %>% unique()# %>% filter(pid %in% ppi_db$receptor_pid)

      link_df = data.frame()
      for (i in 1:dim(from_df)[1]){
          temp_from_df = from_df[i,]
          names(temp_from_df) = paste0("from_",names(temp_from_df))

          temp_to_df = to_df
          names(temp_to_df) = paste0("to_",names(temp_to_df))

          temp_link_df = cbind(temp_from_df, temp_to_df)
          link_df = rbind(link_df, temp_link_df)                        
      }
      dim(link_df) #93330                
      link_df$from_to_pid = paste0(link_df$from_pid,"_",link_df$to_pid)
      link_df$report_pair = "not_reported"

      # load known ppi db
      ppi_db = read.delim("01_data/unique_lr_mouse_1928_20240417.csv",header=T,sep=",")        
      link_df[which(link_df$from_to_pid %in% ppi_db$lr_pair), "report_pair"] = "reported"
      
      table(link_df$report_pair) # 77
      write.table(link_df,file="Identified_Sec_PM.txt",sep="\t",row.names = FALSE)
      # 只保留鉴定到的Lig以及Rec

      reported_link_df = link_df %>% filter(from_pid %in% ppi_db$ligand_pid) %>% filter(to_pid %in% ppi_db$receptor_pid)
      dim(reported_link_df)
      table(reported_link_df$report_pair) # 77
      
      write.table(reported_link_df,file="Identified_PPI_1056.txt", sep="\t", row.names=FALSE,quote=FALSE)
      # write.csv(reported_link_df,file="Identified_PPI_reported_77.csv")
      head(reported_link_df)
}
