
library(dplyr)
library(ggplot2)
library(RColorBrewer)

library(clusterProfiler)
# library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(limma)
library(pheatmap)
library(ggpubr)

# Part 2, Circos plot
{
    library(circlize)
    rm(list=ls())        
    obj_list = readRDS(file = "01_data/MouseKPC_spatialdc_v7_ct10_v2_20240812.RDS")
    row.names(obj_list$anno_df) = obj_list$anno_df$pid
    anno_df = obj_list$anno_df
    head(anno_df)

    ct_order = c("PCC","CAF","T4","B","NEU","DC")
    loci_order = c("Sec","PM")

    ct_loci_order = c()
    ct_loci_color = c()
    ct_loci_link_color = c()
    for (i in 1:length(ct_order)){
        ct_loci_order = c(ct_loci_order, paste0(ct_order[i], "_", loci_order))
        for (j in 1:length(loci_order)){            
            ct_loci_color = c(ct_loci_color, ggsci::pal_npg("nrc",alpha = 1)(10)[i])
            ct_loci_link_color = c(ct_loci_link_color, ggsci::pal_npg("nrc",alpha = 0.5)(10)[i])
        }
    }
            

    # node_df 
    node_df = obj_list$identified_SP_df
    node_df = node_df[!is.na(node_df$loci),]
    table(node_df$celltype)
    node_df$ct_loci = paste0(node_df$celltype,"_",node_df$loci)
    node_df$ct_loci = factor(node_df$ct_loci, levels=ct_loci_order)
    node_df = node_df %>% arrange(ct_loci)
    # node_df = node_df %>% filter(abundance > 5)

    order_id = c()
    for (i in ct_loci_order){
        temp_df = node_df %>% filter(ct_loci == i)
        order_id = c(order_id, 1:dim(temp_df)[1])
    }
    node_df$order_id = order_id

    # link_df 
    # add ppi information
    ppi_db = read.delim("01_data/unique_lr_mouse_1928_20240417.csv",header=T,sep=",")
    head(ppi_db)
    node_df$genename = anno_df[node_df$pid,"gene"]

    from_df = node_df %>% filter(loci == "Sec") %>% filter(pid %in% ppi_db$ligand_pid)
    to_df = node_df %>% filter(loci == "PM") %>% filter(pid %in% ppi_db$receptor_pid)

    
    link_df = data.frame()
    for (i in 1:dim(from_df)[1]){
        temp_from_df = from_df[i,]
        names(temp_from_df) = paste0("from_",names(temp_from_df))

        temp_to_df = to_df
        names(temp_to_df) = paste0("to_",names(temp_to_df))

        temp_link_df = cbind(temp_from_df, temp_to_df)
        link_df = rbind(link_df, temp_link_df)
                    
    }
    

    link_df = link_df[link_df$from_celltype != link_df$to_celltype,]
    link_df$from_to_pid = paste0(link_df$from_pid,"_",link_df$to_pid)
    link_df = link_df %>% filter(from_to_pid %in% ppi_db$lr_pair)
    dim(link_df) # 1980
    head(link_df)

    link_df$sig_order = 1
    link_df$sig_order[which(is.na(link_df$from_sig_label) | is.na(link_df$to_sig_label))] = 0
    tail(link_df)
    link_df_order = link_df %>% arrange(sig_order)

    # begin_plot
    pdf(paste0("Sec_PM_ct6_Circlize_",Sys.Date(),".pdf"),width=7,height=7) 

        circos.par(track.height=0.01,cell.padding=c(0.001, 1.00, 0.001, 1.00))
        circos.initialize(node_df$ct_loci, x=node_df$order_id)

        ## track 1 and track 2
        dep_link_df = link_df %>% filter(from_sig_label=="sig") %>% filter(to_sig_label=="sig") 
        dep_from_link_df = dep_link_df %>% dplyr::select(from_ct_loci, from_order_id, from_genename)  %>% unique()
        dep_to_link_df = dep_link_df %>% dplyr::select(to_ct_loci, to_order_id, to_genename)  %>% unique()
        circos.labels(c(dep_from_link_df$from_ct_loci,dep_to_link_df$to_ct_loci), 
            x = c(dep_from_link_df$from_order_id,dep_to_link_df$to_order_id), 
            labels = c(dep_from_link_df$from_genename,dep_to_link_df$to_genename), 
            side = "outside",col="black", font=1,cex=1#,connection_height = mm_h(10)
                )

        ## track 3
        # cell types
        # Levels: 6
        set_track_gap(mm_h(1))
        circos.track(ylim = c(0, 1),track.height=0.01,bg.border = NA,track.index=3)
        
        sectors = ct_loci_order            
        sectors = factor(sectors, levels = sectors)
        sectors_text = gsub(sectors,pattern="_.*",replacement="")
        for (i in 1:length(ct_order)){
            j = which(sectors_text == ct_order[i])
            highlight.sector(sectors[j],track.index=3, col = ggsci::pal_npg("nrc",alpha = 1)(10)[i],text.vjust="5mm" ,#cex = 1, 
                            text = ct_order[i], text.col = "black", niceFacing = TRUE,facing = "bending.inside")

        }
        
        ## track 4 celltype vs Sec - PM
        # Levels: 12
        set_track_gap(mm_h(0.5)) # 2mm
        circos.track(ylim = c(0, 1),track.height=0.1,bg.border = NA,track.index=4)
        sectors_text = gsub(sectors,pattern=".*_",replacement="")
        sectors_text = gsub(sectors_text,pattern="Sec",replacement="Ligand")
        sectors_text = gsub(sectors_text,pattern="PM",replacement="Receptor")

        for (i in 1:length(sectors)){
            highlight.sector(sectors[i],  col = ct_loci_color[i], track.index=4,
                        text = sectors_text[i], cex = 0.8, text.col = "white", niceFacing = TRUE,facing = "bending.inside"
            )
        }

                                        
        ## track 5 proteins abundance
        # abundance / 20
        set_track_gap(mm_h(0.5)) # 2mm
        circos.track(ylim = c(0, 0.7), track.height = 0.1,bg.border = NA,track.index=5,
                    factors=sectors, panel.fun = function(x, y) {
            temp_df = node_df %>% dplyr::filter(ct_loci == CELL_META$sector.index)
        
            for(i in 1:dim(temp_df)[1]){
                # temp_gene = temp_df[i,"Gene"]
                temp_order = temp_df[i,"order_id"]
                temp_log2FC = temp_df[i,"abundance"]
                circos.barplot(temp_log2FC/20, temp_order, col = "#696969",bar_width=0.7,border = NA)
            }        

            
            circos.axis(h = "bottom", direction = "inside",
                        labels.facing = "reverse.clockwise",
                        labels = NULL,#major.tick.length = 10,
                        minor.ticks = 0,
                        # major.at = c(0:dim(temp_df)[1])
                        major.tick = FALSE
                        )        
        })

                           
        ## track 6 link between proteins
        set_track_gap(mm_h(0.2)) # 2mm

        for(i in 1:dim(link_df_order)[1]){
            temp_from_order_id = link_df_order[i,"from_order_id"]
            temp_to_order_id = link_df_order[i,"to_order_id"]

            temp_from_sector = link_df_order[i,"from_ct_loci"]
            temp_to_sector = link_df_order[i,"to_ct_loci"]

            temp_sig_order = link_df_order[i,"sig_order"]


            if(temp_sig_order == 1){#sig
                temp_from_color = ct_loci_link_color[which(temp_from_sector == sectors)]
                circos.link(sector.index1=temp_from_sector, point1 = temp_from_order_id,#c(temp_from_order_id,(temp_from_order_id+0.01)), 
            sector.index2=temp_to_sector, point2 = temp_to_order_id,#c(temp_to_order_id,temp_to_order_id+0.01), 
            lwd = 1, col=temp_from_color)                    

            }else{
                temp_from_color = "grey"
                circos.link(sector.index1=temp_from_sector, point1 = temp_from_order_id,#c(temp_from_order_id,(temp_from_order_id+0.01)), 
            sector.index2=temp_to_sector, point2 = temp_to_order_id,#c(temp_to_order_id,temp_to_order_id+0.01), 
            lwd = 0.001, col=temp_from_color)                    

            }

        }
            

    dev.off()

}
