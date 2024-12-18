# Visualization of reference data
# Yuan
# 20241213
# Fig. 5b

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

# ----------------------------------
# Visualization of reference data
# Fig. 5b
os.chdir("")

sc_adata_df = pd.read_csv("ct10_adata_impute_v2.csv",sep=",")
sc_adata_df.index = sc_adata_df["pid"]

sc_meta_df = pd.read_csv("sc_adata_meta.csv",sep=",",index_col=0)
sc_meta_df = sc_meta_df[sc_meta_df.index.isin(sc_adata_df.columns)] # keep ordered meta_df

sc_adata_df = sc_adata_df[sc_meta_df.index] # 5788 * 50

# add var annotation
var_df = pd.read_csv("pid_gene_17207_mmu.csv")
sc_var_df = pd.DataFrame(index=sc_adata_df.index)

sc_var_df = pd.merge(sc_var_df,var_df,on="pid",how="left")
sc_var_df.index = sc_var_df["pid"]

sc_adata = sc.AnnData(X = csr_matrix(sc_adata_df.T),obs = sc_meta_df,var=sc_var_df)
sc_adata.write_h5ad("MousePDAC2023_impute_CellType10.h5ad")

ct_order = ['PCC', 'CAF', 'T4', 'T8', 'Treg', 'B', 
       'NEU', 'MAC', 'MO', 'DC']
sc_adata.obs["celltype"] = sc_adata.obs["celltype"].astype("category")
sc_adata.obs["celltype"] = sc_adata.obs["celltype"].cat.set_categories(ct_order,ordered=False)

# vis with UMAP
sc.pp.normalize_total(sc_adata)
sc.pp.log1p(sc_adata)

sc.pp.neighbors(sc_adata, n_pcs=20, n_neighbors=10)
sc.tl.umap(sc_adata)

sc.pl.umap(sc_adata, color=["celltype"], size=1000, save="_MousePDAC2023_impute_ct10_5788.pdf")
