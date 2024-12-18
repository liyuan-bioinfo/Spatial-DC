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
# spatial proteomics distance analysis
# sp_adata = sc.read_h5ad("01_data/KPC2023_Spatial_Spots108_protein3607_add_anno.h5ad")

os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/MouseKPC2023/02_proteomic_profiles/")
scp_adata = sc.read_h5ad("03_benchmark_methods/exp_ct9_impute_v2/SpatialDC/SpatialDCmodel_norm_log1p_scale_sample_10_100/reconstruct_not_norm.h5ad")
# sc.pp.normalize_total(scp_adata)

ct_order = ["PCC","CAF","B","T4","NEU","DC"]
# scp_adata = scp_adata[scp_adata.obs["cellperc"] > 0.01]
sc.pp.log1p(scp_adata)
# sc.pp.filter_cells(scp_adata, min_counts=10)
# sc.pp.filter_genes(scp_adata, min_cells=3)

df = scp_adata.var
df['pids'] = df.apply(lambda row: str(row['gene']) +"_"+ str(row['pid']), axis=1)
scp_adata.var.index = scp_adata.var["pids"]

CAF_adata = scp_adata[scp_adata.obs["celltype"] == "CAF",:]

PCC_adata = scp_adata[scp_adata.obs["celltype"] == "PCC",:]

sp_adata = sc.read_h5ad("01_data/KPC2023_Spatial_Spots108_protein3607_add_anno_impute.h5ad")
sp_adata.obs["from"] = CAF_adata.to_df()["Col8a1_Q00780"].values
sp_adata.obs["to"] = PCC_adata.to_df()["Itga2_Q62469"].values

scaler = MinMaxScaler(feature_range=(0, 1))

spatial = pd.DataFrame(sp_adata.obsm['spatial'], index=sp_adata.obs_names,columns=["x","y"])
scale_df = sp_adata.obs[['from', 'to','from_to']]


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


x = scale_df['from']
y = scale_df['to']

fit = np.polyfit(x, y, 1)
fit_fn = np.poly1d(fit)

# 绘制散点图和拟合线
plt.plot(x, y, 'o',color="black")
plt.plot(x, fit_fn(x), '--k',color="red")
