# Visualization of spatially correlated and cell-type significant proteins from both spatial and cell-type data
# Yuan
# 20241213
# Fig. 5e and Fig. S10g

library(dplyr)
# ---------------------------------------------------------------------------------------------------
# 1 - Create RDS file for storing the sig. proteins from reference data and highly corr. proteins from spatial data
{
    setwd("")
    # Preprocess of Reference data
    # Create RDS file
    {
        rm(list=ls())
        # load group_meta
        meta_df = read.csv("01_data/MousePDAC2023_impute_CellType10_meta.csv", header=T, row.names=1) # 50 * 4        
        ct_order = c('PCC', 'CAF', 'T4', 'T8', 'Treg', 'B', 'NEU', 'MAC', 'MO', 'DC')
        meta_df = meta_df %>% dplyr::filter(celltype %in% ct_order)
        meta_df$SampleID = row.names(meta_df)
        meta_df$celltype = factor(meta_df$celltype, levels = ct_order)
        meta_df = meta_df %>% arrange(celltype)

        # load intensity file
        data_df = read.csv("01_data/ct10_adata_impute_v2.csv", header=T, row.names=1, check.names=FALSE) %>% as.data.frame() # 5788 * 50
        data_df = log1p(data_df[,meta_df$SampleID])
        dim(data_df) # 5788 * 12

        # calculate the mean value of cell type
        data_mean_df = data.frame(row.names=row.names(data_df))
        for (i in ct_order){
            temp_df = rowMeans(data_df[meta_df$celltype==i])
            data_mean_df[i] = temp_df
        }
        
        # load annotation file
        anno_df = read.csv("01_data/pid_gene_17207_mmu.csv", header=T, row.names=1, check.names=FALSE)
        anno_df$pid = row.names(anno_df)
        
        # save Object
        obj_list = list()
        obj_list$data_df = data_df
        obj_list$data_mean_df = data_mean_df
        obj_list$meta_df = meta_df
        obj_list$ct_order = ct_order
        obj_list$anno_df = anno_df
     
        saveRDS(obj_list,file = "01_data/MousePDAC_reference_ct10.RDS")
    }

    # Preprocess of Reference data
    # Calculate signficance of all proteins with one-way ANOVA
    {   
        rm(list=ls())
        obj_list = readRDS(file = "01_data/MousePDAC_reference_ct10.RDS")
        meta_df = obj_list$meta_df        
        data_df = obj_list$data_df
        data_mean_df = obj_list$data_mean_df    
        protein_num = dim(data_df)[1]

        enrich_fc = 1.2
        ct_num = 10
        
        # to find enriched cell-types with mean abundance
        dep_df = data.frame()
        for(i in 1:ct_num){
            temp_df = data_mean_df[,i] - data_mean_df[,1:ct_num]
            temp_pid = names(which(rowSums(temp_df>log2(enrich_fc)) == (ct_num-1)))
            temp_df2 = data.frame(pid=temp_pid)
            temp_df2$celltype = names(data_mean_df)[i]
            temp_df2$log2FC = apply(temp_df[temp_pid,-i],1,mean)        
            dep_df = rbind(dep_df,temp_df2)            
        }
        
        ## pvalue, with each sample
        row.names(dep_df) = dep_df$pid
        data_df = na.omit(data_df[dep_df$pid,]) #1475
        protein_num = dim(data_df)[1] #update for enriched proteins

        Pvalue = c()        
        for(i in 1:protein_num){
            pid_df = data_df[i,] %>% t() %>% as.data.frame()
            names(pid_df) = "pid"
            pid_df$SampleId = row.names(pid_df)
            pid_df$celltype = meta_df$celltype

            model = aov(data=pid_df,pid~celltype)            
            sig_df = TukeyHSD(model)$celltype #which ct is enriched
            enriched_ct = dep_df[i,"celltype"]
            enriched_ct_sig_df = sig_df[grep(row.names(sig_df),pattern=paste0("-",enriched_ct,"$|","^",enriched_ct,"-")),] #whethor this ct is sig. with other cts
            temp_pvalue = RecordTest::fisher.method(as.numeric(enriched_ct_sig_df[,"p adj"]))$"p.value"[[1]]                    
            Pvalue = c(Pvalue,temp_pvalue)         
        }        

        dep_df$pvalue = Pvalue
        dep_df$fdr = p.adjust(Pvalue,method = "BH")                        
        dep_df$gene = obj_list$aov_ct10_df$anno_df[dep_df$pid,"gene"]

        obj_list$aov_ct10_df = dep_df
        saveRDS(obj_list, file = "01_data/MousePDAC_reference_ct10.RDS")
    }

    # Preprocess of Predicted Corr. data
    {
        rm(list=ls())
        data_df = read.csv("04_analysis/spatial_potential_markers/spatial_corr_proteins_SpatialDC_v20241127.csv",row.names=1)        
        data_df$pid = gsub(x=data_df$Protein,pattern=".*_",replacement="")
        dim(data_df)
        
        data_df_spread = data_df %>% tidyr::spread(key="CellType",value="Corr")
        ct_order = c('PCC', 'CAF', 'T4', 'B', 'NEU', 'DC')
        row.names(data_df_spread) = gsub(data_df_spread$Protein, pattern=".*_",replacement="") # 3607
        data_df_spread = na.omit(data_df_spread[,ct_order])
                
        enrich_corr = 0.2 # PCC cut-off 
        ct_num = 6
        
        # to find enriched cell-types with mean abundance
        dep_df = data.frame()
        for(i in 1:ct_num){
            temp_df = data_df_spread[,i] - data_df_spread[,1:ct_num]
            temp_pid = names(which(rowSums(temp_df>enrich_corr) == (ct_num-1)))
            if(length(temp_pid)==0){
                next
            }
            temp_df2 = data.frame(pids=temp_pid)
            temp_df2$celltype = names(data_df_spread)[i]
            temp_df2$corr = data_df_spread[temp_pid,i]
            dep_df = rbind(dep_df,temp_df2)            
        }
        
        obj_list = readRDS(file = "01_data/MousePDAC_reference_ct10.RDS")
        obj_list$pred_corr = dep_df
        obj_list$pred_corr_spread = data_df_spread
        saveRDS(obj_list, file = "01_data/MousePDAC_reference_ct10.RDS") # cp to 04_analysis
    }
}

# ---------------------------------------------------------------------------------------------------
# 2 - Visualization of highly corr. proteins for both spatial data and reference data.
# Fig. 5e and Fig. S10g
{
    rm(list=ls())
    setwd("")

    # vis of pred_corr
    {
        obj_list = readRDS(file = "01_data/MousePDAC_reference_ct10.RDS")        

        # keep intersected pid from ct 6
        dep_df = obj_list$aov_ct10_df %>% filter(log2FC > log2(1.2)) %>% filter(fdr < 0.05) # 1334
        pred_df = obj_list$pred_corr %>% filter(corr > 0.7) # 233
        pred_df$gene = obj_list$anno_df[pred_df$pids,"gene"]
        
        ct_order = c('PCC', 'CAF', 'DC', 'NEU', 'T4', 'B')        
        dep_df$ID = paste0(dep_df$celltype, "_", dep_df$pid)
        pred_df$ID = paste0(pred_df$celltype, "_", pred_df$pids)
        
        intersected_ID = intersect(dep_df$ID, pred_df$ID) # 51        
        pred_df = pred_df[which(pred_df$ID %in% intersected_ID),]
        pred_df$celltype = factor(pred_df$celltype, levels=c("PCC","CAF","DC"))
        pred_df = pred_df %>% arrange(celltype, desc(corr)) # 51
                                
        # plot pred corr
        pred_plot_data = obj_list$pred_corr_spread
        pred_plot_data = pred_plot_data[pred_df$pids,ct_order]
        row.names(pred_plot_data) = pred_df$gene

        # vis of pred_corr as heat map
        range_all <- range(c(-1, 1))    
        my_palette <- colorRampPalette(c("navy", "grey","white","#ffb74d","red"))(n=100)
        # my_palette <- colorRampPalette(c("navy", "white","firebrick3"))(n=100)        
        breaks = seq(range_all[1], range_all[2], length.out = 101)
        p1=pheatmap::pheatmap(pred_plot_data,scale="none", color=my_palette, 
        breaks=breaks,cluster_cols=F,cluster_rows=F,silence=T,show_rownames=T,fontsize=10)

        # Fig. 5e
        pdf(paste0("04_analysis/vis_pred_markers_ct6_corr7_cutoff2_",Sys.Date(),".pdf"),height=3,width=3)
        print(p1)
        dev.off()

        # vis of related protein as heat map
        ref_plot_data = obj_list$data_mean_df
        ref_plot_data = ref_plot_data[pred_df$pids,]

        scale_ref_plot_data = t(apply(ref_plot_data, 1, scale)) %>% as.data.frame()
        names(scale_ref_plot_data) =  names(ref_plot_data)    
        head(scale_ref_plot_data)
        row.names(scale_ref_plot_data) = pred_df$gene
        
        # my_palette = colorRampPalette(c("navy", "white","firebrick3"))(n=100)
        range_all = range(c(-2, 2))            
        breaks = seq(range_all[1], range_all[2], length.out = 101)

        p2=pheatmap::pheatmap(scale_ref_plot_data,scale="none",cluster_rows = F,cluster_cols = F,show_colnames=T,fontsize=10,
                    show_rownames = T,#border_color = "white",                    
                    # cellwidth = 5,cellheight = 0.15,
                    silent = T,color = my_palette,breaks = breaks,                
        )

        # Fig. 10g
        pdf(paste0("04_analysis/vis_ref_markers_aov_ct10_",Sys.Date(),".pdf"),height=4,width=4)
        print(p2)
        dev.off()    

    }

    # vis of pred_corr
    # save table
    {
        rm(list=ls())
        setwd("")
        obj_list = readRDS(file = "01_data/MousePDAC_reference_ct10.RDS")        

        # keep intersected pid from ct 6
        dep_df = obj_list$aov_ct10_df %>% filter(log2FC > log2(1.2)) %>% filter(fdr < 0.05) # 1334
        pred_df = obj_list$pred_corr
        pred_df$gene = obj_list$anno_df[pred_df$pids,"gene"]
                
        ct_order = c('PCC', 'CAF', 'T4', 'B', 'NEU', 'DC')        
        dep_df$ID = paste0(dep_df$celltype, "_", dep_df$pid)
        pred_df$ID = paste0(pred_df$celltype, "_", pred_df$pids)
        
        intersected_ID = intersect(dep_df$ID, pred_df$ID) # 51        
        
        pred_df = pred_df[which(pred_df$ID %in% intersected_ID),]
        dep_df = dep_df[which(dep_df$ID %in% intersected_ID),]
        
        output_df = merge(pred_df, dep_df,by="ID",all=TRUE)

        write.csv(output_df, file="04_analysis/sig_corr_ct6.csv")
                            
    }    
}


