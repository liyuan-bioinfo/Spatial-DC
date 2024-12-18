# Keep intersected proteins for comparison of deconvolution methods
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
warnings.filterwarnings("ignore")
import os

# --------------------------------------------------------------------------
os.chdir("")
sc_adata = sc.read_h5ad("MousePDAC2023_impute_CellType10.h5ad")
sp_adata = sc.read_h5ad("MousePDAC2023_impute_Spots108_protein3607.h5ad")

intersected_pid = np.intersect1d(sc_adata.var_names,sp_adata.var_names)
print(intersected_pid.shape) # 2837

# save for benchmark study
sc_adata = sc_adata[:,intersected_pid].copy()
sp_adata = sp_adata[:,intersected_pid].copy()

del sc_adata.var
del sp_adata.var

sc_adata.write_h5ad("for_benchmark/MousePDAC2023_ct10_2837.h5ad")
sp_adata.write_h5ad("for_benchmark/MousePDAC2023_spot108_2837.h5ad")

# # for benchmark 
sc_adata.to_df().to_csv("for_benchmark/MousePDAC2023_ct10_2837_expression.csv")
sc_adata.obs.to_csv("for_benchmark/MousePDAC2023_ct10_meta.csv")
sp_adata.to_df().to_csv("for_benchmark/MousePDAC2023_spot108_2837_expression.csv")
sp_adata.obs[["X","Y"]].to_csv("for_benchmark/MousePDAC2023_spot108_meta.csv")