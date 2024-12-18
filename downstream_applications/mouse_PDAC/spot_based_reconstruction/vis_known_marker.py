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

# Panel-2, known markers analysis
# this replace R to calculate ANOVA sig. 
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
# scp_adata.raw = scp_adata.copy()
# sc.pp.log1p(scp_adata)

ct_order = ["PCC", "CAF", "T4", "B","NEU", "DC"]

scp_adata.obs["celltype"] = scp_adata.obs["celltype"].astype("category")
scp_adata.obs["celltype"] = scp_adata.obs["celltype"].cat.set_categories(ct_order,ordered=False)


markers_list=["Krt7_Q9DCV7","Clu_Q06890","Krt19_P19001","Epcam_Q99JW5",#PCC
    "Ptgis_O35074","Col15a1_O35206","Col5a1_O88207","Tagln_P37804",#CAF
    "Cd5_P13379","Zap70_P43404",#T4
    "Cd74_P04441","Ighm_P01872",#B
    "S100a9_P31725","Mpo_P11247","Elane_Q3UP87",#NEU
    "Mrc1_Q61830"

]

scp_adata = scp_adata[:,markers_list].copy()
scp_adata.var.index = scp_adata.var["gene"]

markers_list=["Krt7","Clu","Krt19","Epcam",#PCC
    "Ptgis","Col15a1","Col5a1","Tagln",#CAF
    "Cd5","Zap70",#T4
    "Cd74","Ighm",#B
    "S100a9","Mpo","Elane",#NEU
    "Mrc1"

]

ax=sc.pl.matrixplot(scp_adata, markers_list,'celltype',return_fig=True,dendrogram=False,swap_axes="False",figsize=[5,4],cmap="Blues",standard_scale="var")
ax.edge_color="none"
# ax.show()
ax.savefig("figures/markers_ANOVA_selected_pheatmap.pdf")
# merged_df.to_csv("ct10_impute_v2_dep_anno.csv")
