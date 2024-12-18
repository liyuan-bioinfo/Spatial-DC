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

sc.pl.spatial(sp_adata, spot_size=1,color=pids,cmap=my_cmap,vmin=-1.5,vmax=1.5,title=pids_gene,save="ct3_markers.pdf")
