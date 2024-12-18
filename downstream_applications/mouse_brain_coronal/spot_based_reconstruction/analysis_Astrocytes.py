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

os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/MouseBrain2022/02_proteomic_profiles")
scp_adata = sc.read_h5ad("03_benchmark_methods/exp_ct4/SpatialDC/SpatialDCmodel_norm_log1p_scale_sample_10_100/reconstruct_not_norm.h5ad")

# set zero of cellperc<0
set_zero = scp_adata.obs["cellperc"] < 0.01
temp_arr = np.array(scp_adata.to_df())
temp_arr[set_zero,:]=0
scp_adata.X = temp_arr

sender_cell = "Astrocytes"
receiver_cell = "Neurons"
pathway_name = "EPHB"

scp_adata.var.index = scp_adata.var["gene"]
scp_adata.var_names = scp_adata.var_names.astype(str)
scp_adata.var_names_make_unique()

# 保留EPHB pathway
df_cellchat = ct.pp.ligand_receptor_database(species='mouse', signaling_type=None, database='CellChat')
df_cellchat = ct.pp.filter_lr_database(df_cellchat, scp_adata, min_cell_pct=0.05) #40
df_cellchat_filtered = df_cellchat[df_cellchat.iloc[:,2] == pathway_name]

sender_pid = list(set(df_cellchat_filtered.iloc[:,0]))
receiver_pid = list(set(df_cellchat_filtered.iloc[:,1]))

# for sender
keep_index = scp_adata.var["gene"].isin(sender_pid + receiver_pid)
keep_ct = scp_adata.obs["celltype"].isin([sender_cell])
target_adata = scp_adata[keep_ct,keep_index].copy()

target_adata.var.index = target_adata.var["gene"]
temp_df = target_adata.to_df()

# for receiver
keep_index = scp_adata.var["gene"].isin(receiver_pid)
keep_ct = scp_adata.obs["celltype"].isin([receiver_cell])
temp_df[receiver_pid]  = scp_adata[keep_ct,keep_index].to_df().values #replace

target_adata.X = np.array(temp_df)

# prepare adata
sc.pp.log1p(target_adata)
target_adata.obs[["Y"]] = -target_adata.obs[["Y"]]
target_adata.obsm["spatial"] = np.array(target_adata.obs[["X","Y"]])

temp_index = [f"spot_{i}" for i in range(target_adata.shape[0])]
target_adata.obs.index = temp_index


ct.tl.spatial_communication(target_adata,
    database_name='cellchat', df_ligrec=df_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True) 

# 展示信号传递的方向
ct.tl.communication_direction(target_adata, database_name='cellchat', pathway_name=pathway_name, k=5)

ct.pl.plot_cell_communication(target_adata, database_name='cellchat', pathway_name='EPHB',normalize_v=True,normalize_v_quantile=0.995,summary="sender",
grid_density=0.4,scale=1,ndsize=0,grid_thresh=10,grid_scale=10,filename=f"analysis/CCC/{sender_cell}_{receiver_cell}_EPHB_sender_0.5.pdf",plot_method="cell")

ct.pl.plot_cell_communication(target_adata, database_name='cellchat', pathway_name='EPHB',normalize_v=True,normalize_v_quantile=0.995,clustering=None,summary="receiver",
grid_density=0.4,scale=1,ndsize=0,grid_thresh=10,grid_scale=10,filename=f"analysis/CCC/{sender_cell}_{receiver_cell}_EPHC_receiver_0.5.pdf",plot_method="cell")

target_adata.obs[["Y"]] = -target_adata.obs[["Y"]]
target_adata.obsm["spatial"] = np.array(target_adata.obs[["X","Y"]])
target_adata.obs[f's-{pathway_name}'] = target_adata.obsm['commot-cellchat-sum-sender'][f's-{pathway_name}']
target_adata.obs[f'r-{pathway_name}'] = target_adata.obsm['commot-cellchat-sum-receiver'][f'r-{pathway_name}']
sc.pl.spatial(target_adata,spot_size=1,color=[f's-{pathway_name}',f'r-{pathway_name}'],cmap="coolwarm",save=f"{pathway_name}_strength.pdf")
