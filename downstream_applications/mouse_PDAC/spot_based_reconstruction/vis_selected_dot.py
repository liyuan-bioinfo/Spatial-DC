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
colors = ['navy',"white", 'red']#['black','navy',"orange", 'red']
my_cmap = mcolors.LinearSegmentedColormap.from_list('mycmap', colors)

os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/MouseKPC2023/02_proteomic_profiles")
scp_adata = sc.read_h5ad("03_benchmark_methods/exp_ct10_impute_v2/SpatialDC/SpatialDCmodel_norm_log1p_scale_sample_10_100/reconstruct_norm.h5ad")
sc.pp.filter_cells(scp_adata, min_counts=10)
sc.pp.filter_genes(scp_adata, min_cells=3)
scp_adata = scp_adata[scp_adata.obs["cellperc"] > 0.01]

scp_adata = scp_adata[scp_adata.obs["celltype"] != "Treg"]
scp_adata = scp_adata[scp_adata.obs["celltype"] != "MAC"]
scp_adata = scp_adata[scp_adata.obs["celltype"] != "T8"]
scp_adata = scp_adata[scp_adata.obs["celltype"] != "MO"]

df = scp_adata.var
df['pids'] = df.apply(lambda row: str(row['gene']) +"_"+ str(row['pid']), axis=1)
scp_adata.var.index = scp_adata.var["pids"]

sc.pp.normalize_total(scp_adata)
scp_adata.raw = scp_adata.copy()
sc.pp.log1p(scp_adata)

ct_order = ["PCC", "CAF", "T4", "B","NEU", "DC"]

scp_adata.obs["celltype"] = scp_adata.obs["celltype"].cat.set_categories(ct_order,ordered=False)

PPI_markers = ["Col8a1_Q00780","Itga2_Q62469"]


# sc.pl.stacked_violin(scp_adata, PPI_markers, groupby="celltype",figsize=[5,3],standard_scale="var",swap_axes="True",vmax=1,cmap=my_cmap)
sc.pl.dotplot(scp_adata, PPI_markers, groupby='celltype', dendrogram=True,standard_scale="var",swap_axes="True",cmap=my_cmap,save="PPI_selected_Col8a1_Itga2.pdf")
