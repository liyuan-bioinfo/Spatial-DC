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

# 展示四个主要区域的相对分布,使用高相关性marker的平均值
import cell2location
scaler = MinMaxScaler(feature_range=(0, 1))

os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/MouseBrain2022/01_cell_perc/")
sp_adata = sc.read_h5ad("01_data/MouseBrain2022_spot208_impute.h5ad")
sc.pp.normalize_total(sp_adata)
sc.pp.log1p(sp_adata)

pids = ["Q9JJV5", "P97441","Q8CC35"] #Cortex
sp_adata.obs["Cortex"]= np.mean(sp_adata.to_df()[pids],axis=1)
   
pids = ["P51830", "Q8VHW2","Q0VE82","P23818","Q61097","P62748","Q9QVP9","Q3UH99"] #Hip
# pids = ["Q8VHW2","P23818"] #Hip
sp_adata.obs["Hip"]= np.mean(sp_adata.to_df()[pids],axis=1)            

pids = ["Q3TVA9", "P28867","Q91YE8"] #Tha 
# pids = ["Q3TVA9", "Q91YE8"] #Tha
sp_adata.obs["Thalamus"]= np.mean(sp_adata.to_df()[pids],axis=1)            

pids = ["P12961", "Q99P58","Q9QXV0"] #Hypo
# pids = ["P12961", "Q9QXV0"] #Hypo
sp_adata.obs["Hypo"]= np.mean(sp_adata.to_df()[pids],axis=1)            

sp_adata.obs[["Cortex", "Hip", "Thalamus", "Hypo"]]= scaler.fit_transform(sp_adata.obs[["Cortex", "Hip", "Thalamus", "Hypo"]])
# show the co-localization of interacted cells
plt = cell2location.plt.plot_spatial(sp_adata,show_img=False,labels=["Cortex", "Hip", "Thalamus", "Hypo"],
                                     color=["Cortex", "Hip", "Thalamus", "Hypo"],circle_diameter=14.0,reorder_cmap=[1,3,4,0],max_color_quantile=1)
# plt.savefig("figures/region4_distribution.pdf")
plt.show()

# sp_adata.obs[["Cortex", "Hip", "Thalamus", "Hypo"]].to_csv("figures/region4_distribution.csv")
