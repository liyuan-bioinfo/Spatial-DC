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

# Panel-1, cell-type proteomics data
os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/MouseBrain2022/01_cell_perc/")

sc_adata = sc.read_h5ad("01_data/MouseBrain2022_CellType4_impute.h5ad")

sc.pp.normalize_total(sc_adata)
sc.pp.log1p(sc_adata)

sc.pp.neighbors(sc_adata, n_pcs=10, n_neighbors=5)
sc.tl.umap(sc_adata)

# sc.pl.umap(sc_adata, color=["celltype"], size=2000)
sc.pl.umap(sc_adata, color=["celltype"], size=3000, save="MouseBrain2022_impute_ct4_10734.pdf")

#保存为csv
# 获取 X_umap 值
umap_values = sc_adata.obsm["X_umap"]
# 获取 celltype 的 obs 信息
celltype_obs = sc_adata.obs["celltype"]
# 将 X_umap 和 celltype 合并成一个 DataFrame
df = pd.DataFrame(umap_values, columns=["UMAP1", "UMAP2"])
df["celltype"] = celltype_obs.values

df.to_csv("figures/umap_MouseBrain2022_impute_ct4_10734.csv")
