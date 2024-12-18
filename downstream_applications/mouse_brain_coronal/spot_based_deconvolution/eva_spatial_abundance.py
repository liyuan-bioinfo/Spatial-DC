import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr,spearmanr
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import matplotlib as mat
import os
import sys
from scipy import stats
import warnings
from sklearn.preprocessing import MinMaxScaler,StandardScaler
warnings.filterwarnings("ignore")
import os

plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


# 绘制不同细胞之间的共定位特征
# 比较对比方法ct3的主要分布。

import cell2location
scaler = MinMaxScaler(feature_range=(0, 1))

os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/MouseBrain2022/01_cell_perc/")
sp_adata = sc.read_h5ad("01_data/MouseBrain2022_spot208_impute.h5ad")
sc.pp.normalize_total(sp_adata)
sc.pp.log1p(sp_adata)


project_dir = "03_benchmark_methods/exp_ct4/"
methods=["SpatialDC_reconstruct","Tangram","Cell2location","Stereoscope","DestVI","Seurat","spotlight","spatialdwls","CARD"]
# methods=["Seurat"]
ct_order = ["Neurons","Oligodendrocytes","Astrocytes"]

save_PCC = [] # average ct3
save_SPCC = []

for m in methods:
    print(f"{m}\n")
    if(m =="SpatialDC_initial"):
        temp_pred_path = f"{project_dir}/SpatialDC/initial.csv"
    elif(m =="SpatialDC_reconstruct"):
        temp_pred_path = f"{project_dir}/SpatialDC/reconstruct.csv"
    else:
        temp_pred_path = f"{project_dir}/{m}/{m}.csv"
        
    pred_df = pd.read_csv(temp_pred_path,index_col=0)                
    sp_adata.obs[ct_order]= scaler.fit_transform(pred_df[ct_order])            
    plt = cell2location.plt.plot_spatial(sp_adata,show_img=False,labels=ct_order,color=ct_order,circle_diameter=15.0,reorder_cmap=[1,3,4],max_color_quantile=1)

    # 计算mean PCC与SPCC    
    temp_PCC = [] #average markers
    temp_SPCC = []
    for ct in ct_order:
        if(ct == "Neurons"):
            pids = ["Q9JJV5", "P97441","Q8CC35"]

        elif(ct == "Oligodendrocytes"):
            pids = ["P04370", "P60202","P16330","Q60771","Q61885"]
        elif(ct == "Astrocytes"):
            pids = ["P28571","P31650", "Q9JM63"]
        if(np.sum(sp_adata.obs[ct]) == 0):
            temp_PCC.append(0)
            temp_SPCC.append(0)
            continue
        
        corr_matrix = pd.DataFrame(sp_adata.to_df()[pids]).corrwith(sp_adata.obs[ct])  
        corr_matrix = corr_matrix.dropna()
        temp_PCC.append(np.mean(corr_matrix))
        
        corr_matrix2 = pd.DataFrame(sp_adata.to_df()[pids]).corrwith(sp_adata.obs[ct],method="spearman")         
        corr_matrix2 = corr_matrix2.dropna()
        temp_SPCC.append(np.mean(corr_matrix2))
                                    
    save_PCC.append(round(np.mean(temp_PCC),3))
    save_SPCC.append(round(np.mean(temp_SPCC),3))
        
# save_PCC
    plt.show()
    # plt.savefig(f"figures/benchmark/{m}_ct3.pdf")
# results_df = pd.DataFrame({"PCC_CancerCell":save_PCC1, "SPCC_CancerCell":save_SPCC1, "PCC_ImmuneCell":save_PCC2, "SPCC_ImmuneCell":save_SPCC2},index=methods)
# results_df.to_csv("figures/benchmark/cellperc_corr_ct2.csv")
