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

# Panel-1, umap of single-cell spaital proteomics 
# cell-type proteomics data
os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/MouseKPC2023/02_proteomic_profiles")
scp_adata = sc.read_h5ad("03_benchmark_methods/exp_ct10_impute_v2/SpatialDC/SpatialDCmodel_norm_log1p_scale_sample_10_100/reconstruct_norm.h5ad")

sc.pp.filter_cells(scp_adata, min_counts=10)
sc.pp.filter_genes(scp_adata, min_cells=3)

scp_adata = scp_adata[scp_adata.obs["cellperc"] > 0.01]

sc.pp.normalize_total(scp_adata)
sc.pp.log1p(scp_adata)

ct_order = ["PCC", "CAF", "T4", "B","NEU", "DC"]

scp_adata.obs["celltype"] = scp_adata.obs["celltype"].astype("category")
scp_adata.obs["celltype"] = scp_adata.obs["celltype"].cat.set_categories(ct_order,ordered=False)

sc.pp.neighbors(scp_adata, n_pcs=20, n_neighbors=10)
sc.tl.umap(scp_adata)

my_palette = ["#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF","#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF"]
scp_adata.uns["celltype_colors"] = my_palette


sc.pl.umap(scp_adata, color=["celltype"], size=200,save="MouseKPC_Reconstruct_ct6_impute_umap.pdf")
# # sc.pl.umap(sc_adata, color=["celltype"], size=1000, save="MouseKPC_Reference_ct10_impute_umap.pdf")

# #保存为csv
# # 获取 X_umap 值
# umap_values = scp_adata.obsm["X_umap"]
# # 获取 celltype 的 obs 信息
# celltype_obs = scp_adata.obs["celltype"]
# # 将 X_umap 和 celltype 合并成一个 DataFrame
# df = pd.DataFrame(umap_values, columns=["UMAP1", "UMAP2"])
# df["celltype"] = celltype_obs.values

# df.to_csv("figures/umap_ct6.csv")
