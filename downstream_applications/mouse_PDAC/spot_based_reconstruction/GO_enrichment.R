library(dplyr)
library(ggplot2)
library(RColorBrewer)

library(clusterProfiler)
library(org.Mm.eg.db)
library(limma)
library(pheatmap)
library(ggpubr)


 # GO-BP analysis with DEP
{
    library(clusterProfiler)
    Obj_list = readRDS(file="01_data/MouseKPC_spatialdc_v7_ct10_v2_20240812.RDS") 

    meta_df = Obj_list$meta_df# %>% filter(CellType %in% Obj_list$ct3_order)                
    # data_df = Obj_list$data_df[,meta_df$SampleID]
    dep_df = Obj_list$dep_df %>% dplyr::filter(log2FC>log2(2) &fdr<0.05)        
    row.names(Obj_list$anno_df) = Obj_list$anno_df$pid
    dep_df$gene = Obj_list$anno_df[dep_df$pid,"gene"]

    write.csv(dep_df,"filter_ct10_impute_dep.csv")
    
        
    tran_df = bitr(dep_df$pids, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    plot_df = dep_df %>% merge(tran_df, by.y="UNIPROT", by.x="pids",all.x=F,all.y=F)    
    plot_df$celltype = factor(plot_df$celltype, levels=Obj_list$ct_order)    
    plot_df = plot_df %>% dplyr::arrange(celltype)
    
    # bg_df = Obj_list$data_df
    # row.names(bg_df) = gsub(row.names(bg_df),pattern=";.*",replacement="")
    # bg_pid_df = bitr(row.names(bg_df), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    
    formula_res <- compareCluster(data = plot_df,ENTREZID~celltype,fun="enrichGO", 
                                OrgDb = org.Mm.eg.db,minGSSize=3,maxGSSize=2000,
                                ont = "BP", pAdjustMethod = "BH",#universe = bg_pid_df$ENTREZID,
                                pvalueCutoff = 1,qvalueCutoff = 1,readable=TRUE
                                
    )
    formula_res_cutoff = formula_res
    formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$p.adjust<=0.05,]
    
    write.csv(x = formula_res_cutoff@compareClusterResult, file=paste0("analysis/GO/GOBP_",Sys.Date(),".csv"))  
    
    # # selcted GO items
    selected_GO = formula_res_cutoff@compareClusterResult
    selected_items = c("mitochondrial translation","mitochondrial gene expression","cell-cell adhesion mediated by cadherin", # PCC
        "blood vessel development","extracellular matrix organization","extracellular structure organization", #CAF
    "T cell costimulation","positive regulation of alpha-beta T cell proliferation","lymphocyte costimulation",#T4
    "regulation of B cell apoptotic process","B cell homeostasis","B cell apoptotic process",#B
    "myeloid leukocyte migration","neutrophil activation","neutrophil mediated immunity",#NEU
    "receptor-mediated endocytosis","negative regulation of leukocyte mediated immunity","antigen processing and presentation of peptide antigen via MHC class II"#DC
    
    )    
    formula_res_cutoff@compareClusterResult = selected_GO %>% filter(Description %in% selected_items)
    formula_res_cutoff@compareClusterResult$log10P = -log10(formula_res_cutoff@compareClusterResult$pvalue)
    p2=dotplot(formula_res_cutoff, label_format=50,showCategory=3,font.size=14,color="log10P",size="count") + 
        theme(panel.grid = element_blank(),axis.ticks.y = element_blank()) +
        scale_colour_gradientn(colours=colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "YlOrRd")[3:6])(30))

    # save dot plot
    pdf(file=paste0("analysis/GO/ct6_GOBP_",Sys.Date(),".pdf"),width = 8,height = 6)
    print(p2)
    dev.off() 
    
}


# GO-MF analysis with identited proteins
# 使用manual curated 的Sec-PM 数据库对重构的蛋白表达谱进行基本分析
# 1，鉴定到多少的Sec以及PM proteins
# 2，这些SEC 以及 PM的分子功能是什么样的

{
    rm(list=ls())        
    obj_list = readRDS(file = "01_data/MouseKPC_spatialdc_v7_ct10_v2_20240812.RDS")
    
    meta_df = obj_list$meta_df
    impute_df = obj_list$data_df # 2837
    
    impute_df$pid = row.names(impute_df)
    impute_df$genename = obj_list$anno_df[impute_df$pid,"gene"]        

    # 1，鉴定到多少的Sec以及PM proteins
    # Bar plot
    Sec_PM_db = read.delim("01_data/curated_mouse_LRdb_5145.csv",header=T,sep=",")
    
    Sec_db_pid = unique(Sec_PM_db[Sec_PM_db$loci=="Sec","Entry"])
    PM_db_pid = unique(Sec_PM_db[Sec_PM_db$loci=="PM","Entry"])
    
    length(intersect(impute_df$pid, Sec_db_pid)) #153
    length(intersect(impute_df$pid, PM_db_pid)) #610
    
    # 2, 这些Sec以及PM的功能是什么
    identified_PM_df = Sec_PM_db[which(Sec_PM_db$Entry  %in% impute_df$pid),]
    dim(identified_PM_df)#763
    tran_df = bitr(identified_PM_df$Entry, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    identified_PM_df = identified_PM_df %>% merge(tran_df, by.x="Entry", by.y="UNIPROT",all.x=F,all.y=F)

    formula_res <- compareCluster(data = identified_PM_df, ENTREZID~loci,fun="enrichGO", 
                                OrgDb = org.Mm.eg.db,minGSSize=3,maxGSSize=2000,
                                ont = "MF", pAdjustMethod = "BH",
                                pvalueCutoff = 1,qvalueCutoff = 1,#universe = bg_df$ENTREZID,
                                readable=TRUE
    )
    formula_res_cutoff = formula_res
    formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$p.adjust<=0.05,]
    
    write.csv(formula_res_cutoff@compareClusterResult, file="analysis/GO/Sec_PM_GOMF.csv")

    labeled_GO = read.csv("01_data/labeled_Sec_PM_GOMF.csv")        
    category = c("cell_adhesion", "transporter", "integrin", "ion_channel", "growth_factor", "peptidase", "protease", "ECM")
    category_df = data.frame() # 分类求和
    total_category_df = data.frame() # 存放每个ct 所有cg的num

    # 计算每个cluster 每个category的数量
    for (cl in unique(labeled_GO$Cluster)){
        # print(ct)
        temp_df = labeled_GO %>% filter(Cluster == cl)
        genes_ct = c()
        for (cg in category){
            temp_cg_df = temp_df %>% filter(Category == cg)
            genes_ct_cg = c()
            if(dim(temp_cg_df)[1]==1){
                genes_ct_cg = c(genes_ct_cg,strsplit(as.character(temp_cg_df["geneID"]),split = "/")[[1]])
            }else if (dim(temp_cg_df)[1] > 1) {
                for(i in 1:dim(temp_cg_df)[1]){                        
                genes_ct_cg = c(genes_ct_cg,strsplit(as.character(temp_cg_df[i,"geneID"]),split = "/")[[1]])
                }
            }
            # save category
            # print(paste0(cg,":",length(unique(genes_ct_cg))))
            cg_num_df = data.frame("num"=length(unique(genes_ct_cg)),"Category"=cg)
            cg_num_df$Cluster = cl
            category_df = rbind(category_df, cg_num_df)

            genes_ct = c(genes_ct, genes_ct_cg)
        }
            genes_ct_df = data.frame("num"=length(unique(genes_ct)), "Cluster"=cl)
            total_category_df = rbind(total_category_df, genes_ct_df)

            # save category of ct
            print(paste0(cl,":",length(unique(genes_ct))))
        
    }
    # begin plot
    head(category_df)
    # category_df$Cluster = factor(category_df$Cluster, levels=obj_list$ct3_order)
    category_df$Category = factor(category_df$Category, levels=category)
    # category_df$label = category_df$num        
    
    p1 = ggplot(category_df, aes(fill=Category, y=num, x=Cluster, label=num)) + 
        geom_bar(position="fill", stat="identity",width=0.5) + 
        # geom_bar(position="fill", stat="identity", aes(label = num)) + 
        geom_text(position=position_fill(vjust=0.5), size=6) +            
        theme_bw() + xlab("")+ylab("") +
        scale_fill_manual(values = ggsci::pal_npg()(length(category))) +
        theme(plot.background = element_blank(),panel.background = element_blank(),strip.background = element_blank(),
            legend.background = element_blank(),text = element_text(size=12),panel.grid = element_blank()) + coord_flip()
    
    pdf(file=paste0("analysis/GO/Sec_PM_GOMF_category_",Sys.Date(),".pdf"),width=7,height = 2)
    print(p1)
    dev.off()        

    write.csv(category_df, file=paste0("analysis/GO/Sec_PM_GOMF_category_",Sys.Date(),".csv"))


}

