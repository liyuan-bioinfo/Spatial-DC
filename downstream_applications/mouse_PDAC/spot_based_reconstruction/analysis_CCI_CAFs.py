# -------------- 信号通讯分析-------------------
# COMMOT
import sys
import matplotlib.pyplot as plt
import seaborn as sns

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

# Commot
import gc
import ot
import pickle
from scipy import sparse
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt

import commot as ct

os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/MouseKPC2023/02_proteomic_profiles")

scp_adata = sc.read_h5ad("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/MouseKPC2023/02_proteomic_profiles/03_benchmark_methods/exp_ct10_impute_v2/SpatialDC/SpatialDCmodel_norm_log1p_scale_sample_10_100/reconstruct_not_norm.h5ad")

# set zero of cellperc<0
set_zero = scp_adata.obs["cellperc"] < 0.01
temp_arr = np.array(scp_adata.to_df())
temp_arr[set_zero,:]=0
scp_adata.X = temp_arr

# keep target cells, first CAF and then PC
ct_order = ["CAF"]
keep_index = scp_adata.var["gene"].isin(["Col8a1","Itga2"])
keep_ct = scp_adata.obs["celltype"].isin(ct_order)
target_adata = scp_adata[keep_ct,keep_index].copy()


df = target_adata.var
df['pids'] = df.apply(lambda row: str(row['gene']) +"_"+ str(row['pid']), axis=1)
target_adata.var.index = target_adata.var["gene"]
temp_index = [f"spot_{i}" for i in range(target_adata.shape[0])]
target_adata.obs.index = temp_index

temp_df = target_adata.to_df()

ct_order = ["PCC"]
keep_index = scp_adata.var["gene"].isin(["Itga2"])
keep_ct = scp_adata.obs["celltype"].isin(ct_order)
temp_df["Itga2"]  = scp_adata[keep_ct,keep_index].to_df().values

target_adata.X = np.array(temp_df)

# basic filter
# scp_adata = scp_adata[scp_adata.obs["cellperc"] > 0.01]
# sc.pp.filter_cells(scp_adata, min_counts=10)
# sc.pp.filter_genes(scp_adata, min_cells=3)


# set genes
# df = scp_adata.var
# df['pids'] = df.apply(lambda row: str(row['gene']) +"_"+ str(row['pid']), axis=1)
# scp_adata.var.index = scp_adata.var["gene"]

# prepare adata
target_adata.obs[["Y"]] = -target_adata.obs[["Y"]]
target_adata.obsm["spatial"] = np.array(target_adata.obs[["X","Y"]])
# sc.pp.normalize_total(target_adata)
# target_adata.raw = target_adata.copy()
sc.pp.log1p(target_adata)

# scaler = MinMaxScaler(feature_range=(0, 1))
# inputs = scaler.fit_transform(target_adata.to_df())
# target_adata.X = inputs

# 准备PPI
target_ppi_df = pd.DataFrame({0:'Col8a1',1:'Itga2',2:'sig_ppi',3:'Secreted Signaling'},index=[0,1,2,3])
ct.tl.spatial_communication(target_adata,
    database_name='selected_db', df_ligrec=target_ppi_df, dis_thr=500, heteromeric=True, pathway_sum=True) 

# 展示信号传递的方向
ct.tl.communication_direction(target_adata, database_name='selected_db', pathway_name='sig_ppi', k=5)

ct.pl.plot_cell_communication(target_adata, database_name='selected_db', pathway_name='sig_ppi',normalize_v=True,normalize_v_quantile=0.995,summary="sender",
grid_density=0.4,scale=1.5,ndsize=0,grid_thresh=10,grid_scale=10,filename="sender.pdf",plot_method="cell")

# ct.pl.plot_cell_communication(target_adata, database_name='selected_db', pathway_name='sig_ppi',normalize_v=True,normalize_v_quantile=0.995,clustering=None,summary="receiver",
# grid_density=0.4,scale=3,ndsize=50,grid_thresh=10,grid_scale=10)
ct.pl.plot_cell_communication(target_adata, database_name='selected_db', pathway_name='sig_ppi',normalize_v=True,normalize_v_quantile=0.995,clustering=None,summary="receiver",
grid_density=0.4,scale=1.5,ndsize=0,grid_thresh=10,grid_scale=10,filename="receiver.pdf",plot_method="cell")


# # 绘制每个蛋白的空间分布图
adata = target_adata.copy()
pts = adata.obsm['spatial']

adata.obsm['commot-selected_db-sum-sender'][['s-Col8a1-Itga2']]= scaler.fit_transform(adata.obsm['commot-selected_db-sum-sender'][['s-Col8a1-Itga2']])  
adata.obsm['commot-selected_db-sum-receiver'][['r-Col8a1-Itga2']]= scaler.fit_transform(adata.obsm['commot-selected_db-sum-receiver'][['r-Col8a1-Itga2']])  

s = adata.obsm['commot-selected_db-sum-sender']['s-Col8a1-Itga2']
r = adata.obsm['commot-selected_db-sum-receiver']['r-Col8a1-Itga2']


fig, ax1 = plt.subplots(nrows=1, ncols=2, figsize=(6, 2))

sc1 = ax1[0].scatter(pts[:,0], pts[:,1],s=105, edgecolors='none', facecolors='none',marker="s",
        c=s,cmap='coolwarm')
ax1[0].set_aspect('equal', 'box')
ax1[0].set_xticks([])
ax1[0].set_yticks([])

ax1[0].set_title("Sender signal")
ax1[0].axis('off')

sc2 = ax1[1].scatter(pts[:,0], pts[:,1],s=105, edgecolors='none', facecolors='none',marker="s",
        c=r,cmap='coolwarm')
ax1[1].set_aspect('equal', 'box')
ax1[1].set_xticks([])
ax1[1].set_yticks([])

ax1[1].set_title("Receiver signal")
ax1[1].axis('off')

fig.colorbar(sc1, ax=ax1[0])
fig.colorbar(sc2, ax=ax1[1])
plt.tight_layout()
# plt.show()
plt.savefig("figures/Sender_Receiver_Col8a2_Itga2.pdf")

# save table
adata.obsm['commot-selected_db-sum-sender'].to_csv("figures/Sender_Col8a2_Itga2.csv")
adata.obsm['commot-selected_db-sum-receiver'].to_csv("figures/Receiver_Col8a2_Itga2.csv")
