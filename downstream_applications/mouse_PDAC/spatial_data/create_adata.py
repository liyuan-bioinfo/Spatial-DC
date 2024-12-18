# Create the Anndata object for spatial proteomics data
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


warnings.filterwarnings("ignore")

# create spatial proteomics adata
os.chdir("")

sp_adata_df = pd.read_csv("spatial_adata_impute.csv",sep=",",index_col=0)
sp_meta_df = pd.read_csv("sp_adata_meta.csv",sep=",",index_col=0)

# add var annotation
var_df = pd.read_csv("pid_gene_17207_mmu.csv")

sp_var_df = pd.merge(sp_adata_df,var_df,on="pid",how="left")
sp_var_df.index = sp_var_df["pid"]

sp_adata = sc.AnnData(X = csr_matrix(sp_adata_df.T),var=sp_var_df[["pid","gene","genes"]])
sp_adata.obs.index = sp_adata_df.columns
sp_adata.obs[["X","Y"]] = sp_meta_df[["x","y"]]
sp_adata.obsm["spatial"] = np.array(sp_meta_df[["x","y"]])

sp_adata.write_h5ad("MousePDAC2023_impute_Spots108_protein3607.h5ad")
