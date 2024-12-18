# Compare the corr. between ground truth and predicted values. 
# Yuan
# 20241213

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
import matplotlib.colors as mcolors
scaler = MinMaxScaler(feature_range=(0, 1))
import cell2location

warnings.filterwarnings("ignore")

plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# benchmark of cell type composition as cell-type relative distribution
# Fig. 5d, PCC and SPCC metrics

os.chdir("")
scaler = MinMaxScaler(feature_range=(0, 1))
sp_adata = sc.read_h5ad("01_data/MousePDAC2023_spot108_2837.h5ad")

project_dir = "03_benchmark_methods/exp_ct10_impute_v2/" # save benchmark results here
methods = ["SpatialDC","Tangram","cell2location","Stereoscope","DestVI","Seurat","SPOTlight","SpatialDWLS","CARD"] # re-run cell2location

results = []
for m in methods:
    gd_df = pd.read_csv("01_data/gd.csv",index_col=0)        

    if(m =="SpatialDC_initial"):
        temp_pred_path = f"{project_dir}/SpatialDC/initial.csv"
    elif(m =="SpatialDC"):
        temp_pred_path = f"{project_dir}/SpatialDC/reconstruct.csv"
    else:
        temp_pred_path = f"{project_dir}/{m}/{m}.csv"
        
    pred_df = pd.read_csv(temp_pred_path, index_col=0)   

    ct_orders = ['PCC','Immune']    
    pred_df["Immune"] = pred_df[['T4', 'T8', 'Treg', 'B', 'NEU', 'MAC', 'MO', 'DC']].sum(axis=1).values
    
    for ct in ct_orders:
        pcc_corr = pearsonr(gd_df[ct], pred_df[ct])[0]
        spcc_corr = spearmanr(gd_df[ct], pred_df[ct])[0]
        results.append({
        'Method': m,
        'CellType': ct,
        'PCC': pcc_corr,
        'SPCC': spcc_corr
    })
    
    temp_adata = sp_adata.copy()
    temp_adata.obs[["PCC","Immune"]]= scaler.fit_transform(pred_df[["PCC","Immune"]])

    # Fig. 5d right panels    
    plt = cell2location.plt.plot_spatial(temp_adata,show_img=False,labels=["PCC"],color=["PCC"],circle_diameter=20.0,reorder_cmap=[1],max_color_quantile=1)
    plt.savefig(f"figures/benchmark_v20241127/spatial_map_{m}_PC.pdf")
    plt = cell2location.plt.plot_spatial(temp_adata,show_img=False,labels=["Immune"],color=["Immune"],circle_diameter=20.0,reorder_cmap=[3],max_color_quantile=1)
    plt.savefig(f"figures/benchmark_v20241127/spatial_map_{m}_PI.pdf")

# Fig. 5d, eva. metrics bellow predicted cell-type distribution
results_df = pd.DataFrame(results)
results_df.to_csv('figures/benchmark_v20241127/spatial_map_benchmark_correlation_results.csv', index=False)
