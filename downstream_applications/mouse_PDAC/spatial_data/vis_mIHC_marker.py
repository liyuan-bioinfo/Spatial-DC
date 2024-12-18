# Visualization of ground truth distribution for PC and PI
# 20241213
# Yuan
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

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
scaler = MinMaxScaler(feature_range=(0, 1))

import cell2location
scaler = MinMaxScaler(feature_range=(0, 1))

# Visualization of the cell type of ground truth distribution for PC and PI
# Fig. 5c middle and right

os.chdir("")
scaler = MinMaxScaler(feature_range=(0, 1))
sp_adata = sc.read_h5ad("MousePDAC2023_impute_Spots108_protein3607.h5ad")
gd_df = pd.read_csv("gd.csv",index_col=0)
sp_adata.obs = gd_df

sp_adata.obs[["PCC","Immune"]]= scaler.fit_transform(sp_adata.obs[["PCC", "Immune"]])            

# show the distribution of PCC and Immune
plt = cell2location.plt.plot_spatial(sp_adata,show_img=False,labels=["PCs"],color=["PCC"],circle_diameter=20.0,reorder_cmap=[1],max_color_quantile=1)
plt.savefig("figures/MousePDAC_gd_PCC.pdf") # middle panel

plt = cell2location.plt.plot_spatial(sp_adata,show_img=False,labels=["PIs"],color=["Immune"],circle_diameter=20.0,reorder_cmap=[3],max_color_quantile=1)
plt.savefig("figures/MousePDAC_gd_Immune.pdf") # right panel

