# Visualization for prediction cell-type composition using Spatial-DC and eight state-of-the-art deconvolution methods from mouse PDAC data
# Yuan
# 20241213
# Fig. 5f and Fig. 5g

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

import cell2location
import matplotlib.colors as mcolors

import warnings
from sklearn.preprocessing import MinMaxScaler,StandardScaler
warnings.filterwarnings("ignore")

plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ------------------------------------------------------------------------------------
os.chdir("")
sp_adata = sc.read_h5ad("00_raw/MousePDAC2023_impute_Spots108_protein3607.h5ad")

df = sp_adata.var
df['pids'] = df.apply(lambda row: str(row['gene']) +"_"+ str(row['pid']), axis=1)
sp_adata.var.index = sp_adata.var["pids"]

sc.pp.normalize_total(sp_adata)
sc.pp.log1p(sp_adata)

ct_perc_pred = pd.read_csv("01_cell_perc/03_benchmark_methods/exp_ct10_impute_v2/SpatialDC/reconstruct.csv",index_col=0)
ct_perc_pred["Immune"] = np.sum(ct_perc_pred[['T4', 'B', 'NEU', 'DC']],axis=1)

sp_adata.obs[["PCC","Immune"]]= scaler.fit_transform(ct_perc_pred[["PCC", "Immune"]])
pred_PCC = sp_adata.obs["PCC"]

selected_markers = ["Agr2_O88312","Tspan8_Q8R3G9","Clu_Q06890","Fn1_P11276"]
scale_df = sp_adata.to_df()[selected_markers]

# plot
# Fig. 5f
fig, axs = plt.subplots(1, 4, figsize=(12, 3))
for j, ct in enumerate(selected_markers):
    
    p = sp_adata.to_df()[ct].values
    g = pred_PCC

    values = np.vstack([g, p])
    kernel = stats.gaussian_kde(values)(values)
    axs[j].scatter(p, g, cmap="RdYlBu_r",c=kernel,marker='.')
    # axs[j].set_ylabel([])
    axs[j].set_xlabel("log1p(abundance)")
    axs[j].set_title(f"{ct}")

plt.tight_layout()
plt.savefig("01_cell_perc/analysis/vis_markers_scatter_v20241127.pdf")

# plot
# Fig. 5g
colors = ['white', 'grey','#ffb74d', 'red']
my_cmap = mcolors.LinearSegmentedColormap.from_list('mycmap', colors)
scale_df = sp_adata.to_df()[selected_markers]
fig, ax1 = plt.subplots(nrows=1, ncols=4, figsize=(12, 3))
spatial = pd.DataFrame(sp_adata.obsm['spatial'], index=sp_adata.obs_names,columns=["x","y"])
for i in range(4):
    selected_markers_order = selected_markers[i]
    pred_df_1 = scale_df[selected_markers_order].values
    sc_bar_1 = ax1[i].scatter(spatial['x'], -spatial['y'],s=190, edgecolors='none', facecolors='none',marker="s",
            c=pred_df_1,cmap=my_cmap)
    ax1[i].set_aspect('equal', 'box')
    ax1[i].set_xticks([])
    ax1[i].set_yticks([])

    ax1[i].set_title(selected_markers[i])
    ax1[i].axis('off')
    fig.colorbar(sc_bar_1, ax=ax1[i], orientation='vertical', pad=0.02, aspect=30, shrink=0.8)

plt.savefig("01_cell_perc/analysis/vis_markers_distriubtion_v20241127.pdf")