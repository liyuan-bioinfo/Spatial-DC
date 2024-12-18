library(ggplot2)
library(survival)
library("survminer")


rm(list=ls())
setwd("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/MouseKPC2023/02_proteomic_profiles")

# data clean for proteins
{
    library(org.Hs.eg.db)
    rm(list=ls())
    data_df = read.delim("01_data/PDAC_CPTAC/PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt", check.names=FALSE)
    data_df$idx = gsub(data_df$idx,pattern="\\..*",replacement="")

    tran_df = clusterProfiler::bitr(data_df$idx, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
    data_df = data_df %>% merge(tran_df, by.x="idx",by.y="ENSEMBL", all.x=F,all.y=F)    

    meta_df = read.delim("01_data/PDAC_CPTAC/PDAC_survival.txt", check.names=FALSE) %>% dplyr::select(case_id,OS_days,OS_event)
    meta_df = na.omit(meta_df)
    dim(meta_df) # 97    

    filter_df = data_df[,c("SYMBOL",meta_df$case_id)]
    head(filter_df)

    write.csv(filter_df,file="01_data/PDAC_CPTAC/PDAC_protein_filter_97.csv")
    write.csv(meta_df,file="01_data/PDAC_CPTAC/PDAC_survival_filter_97.csv")

}

# Gene
# ITGA2 + COL8A1
{
    rm(list=ls())
    data_df = read.delim("01_data/PDAC_CPTAC/PDAC_gene_filter_97.csv", sep=",",row.names=1, check.names=FALSE)
    data_df = data_df %>% dplyr::filter(SYMBOL %in% c("COL8A1","ITGA2")) %>% na.omit()    
    row.names(data_df) = data_df$SYMBOL
    
    meta_df = read.delim("01_data/PDAC_CPTAC/PDAC_survival_filter_97.csv", sep=",",row.names=1, check.names=FALSE)
    data_df = na.omit(data_df[,meta_df$case_id]) 
    data_df = data_df %>% t() %>% as.data.frame()

    data_df = data_df %>% arrange(desc(COL8A1))
    data_df$rank1 = 1:dim(data_df)[1]
    head(data_df)

    data_df = data_df %>% arrange(desc(ITGA2))
    data_df$rank2 = 1:dim(data_df)[1]

    data_df$abundance = (data_df$rank1 + data_df$rank2)/2
    data_df$case_id = row.names(data_df)
    
    # merge
    cut_off = as.numeric(quantile(data_df$abundance,prob=0.5))    
    # high mean more abundance
    data_df$Group = "high"
    data_df[which(data_df$abundance >= cut_off), "Group"] = "low"
    
    # protein-1
    cut_off = as.numeric(quantile(data_df$COL8A1,prob=0.5))    
    # high mean more abundance
    data_df$Group_1 = "high"
    data_df[which(data_df$COL8A1 <= cut_off), "Group_1"] = "low"

    # protein-2
    cut_off = as.numeric(quantile(data_df$ITGA2,prob=0.5))    
    # high mean more abundance
    data_df$Group_2 = "high"
    data_df[which(data_df$ITGA2 <= cut_off), "Group_2"] = "low"


    plot_df = merge(data_df, meta_df, by="case_id",all=F)
    plot_df$OS_Years = plot_df$OS_days / 365


    filter_plot_df = plot_df %>% filter(Group_1 == Group_2)
    # filter_plot_df$Group_3 = paste0(filter_plot_df$Group_1,"_",filter_plot_df$Group_2)
            
    fit <- survfit(Surv(OS_Years, OS_event) ~ Group_1, data = filter_plot_df) # time; status ~ exp
    
    p1=ggsurvplot(fit,pval = TRUE,
            conf.int = FALSE,
            #risk.table.col = "strata", # Change risk table color by groups
            ggtheme = theme_bw(), # Change ggplot2 theme
            palette = c("#E7B800", "#2E9FDF")#c("#E7B800", "#2E9FDF","")#,
            # xlim = c(0, 600)
            )

    pdf(file = paste0("analysis/clinical_meaning/Survival_curve_COL8A1_ITGA2_intersect_gene_", Sys.Date(), ".pdf"),
        width = 5, height = 4)
    print(p1)
    dev.off()    

}
