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


# 展示Spatial-DC预测出来的细胞比例分布
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
scaler = MinMaxScaler(feature_range=(0, 1))

os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/MouseBrain2022/01_cell_perc/")
colors = ['#001f3f', '#545b62', '#7fc97f', '#ffb74d', 'red']
my_cmap = mcolors.LinearSegmentedColormap.from_list('mycmap', colors)

pids = ["P55088","Q8CC35","P16330"]
pids_gene = ["Aqp4","Synpo","Cnp"]

ct_order = ["Astrocytes","Neurons","Oligodendrocytes"]

# scale spatial data
sp_adata = sc.read_h5ad("01_data/MouseBrain2022_spot208_impute.h5ad")
sc.pp.normalize_total(sp_adata)
sc.pp.log1p(sp_adata)
sc.pp.scale(sp_adata)

# scale pred cell-perc
pred_df = pd.read_csv("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/MouseBrain2022/01_cell_perc/03_benchmark_methods/exp_ct4/SpatialDC/SpatialDCmodel_norm_log1p_scale_sample_10_100/reconstruct.csv",index_col=0)
pred_df[pred_df.columns]= scaler.fit_transform(pred_df[pred_df.columns])

# vis-1
sp_adata.obs = pred_df
sc.pl.spatial(sp_adata, spot_size=1,color=ct_order,cmap="magma",save="SpatialDC_ct3.pdf")

# vis-2
sc.pl.spatial(sp_adata, spot_size=1,color=pids,cmap=my_cmap,vmin=-1.5,vmax=1.5,title=pids_gene,save="ct3_markers.pdf")

# vis-3
sp_adata = sc.read_h5ad("01_data/MouseBrain2022_spot208_impute.h5ad")
sc.pp.normalize_total(sp_adata)
sc.pp.log1p(sp_adata)

fig, axs = plt.subplots(1, 3, figsize=(9, 3))

for j, ct in enumerate(pids):
    
    p = sp_adata.to_df()[ct].values
    g = pred_df[ct_order[j]].values

    temp_pcc = '%.3f' % pearsonr(p,g)[0]
    
    values = np.vstack([g, p])
    kernel = stats.gaussian_kde(values)(values)
    axs[j].scatter(p, g, cmap="RdYlBu_r",c=kernel,marker='.')
    # axs[j].set_ylabel([])
    axs[j].set_xlabel("log1p(abundance)")
    axs[j].set_title(f"{pids_gene[j]}_{temp_pcc}")
    # axs[j].set_yticks([])
#         axs[j, i].set_xlim(0,1)
#         axs[j, i].set_ylim(0,1)
plt.tight_layout()
# plt.show()
plt.savefig("figures/scatter_ct3_markers.pdf", dpi=300, bbox_inches='tight', pad_inches=0.2)
