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
from sklearn.cluster import KMeans
warnings.filterwarnings("ignore")


plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ----------------------------------
# create spatial proteomics adata
os.chdir("")

sp_adata = sc.read_h5ad("MousePDAC2023_impute_Spots108_protein3607.h5ad")

# Clustering of spatial data
# Fig. S10e left
sc.pp.normalize_total(sp_adata)
sc.pp.log1p(sp_adata)

kmeans = KMeans(n_clusters=2, random_state=42).fit(np.array(sp_adata.to_df()))
sp_adata.obs["true_labels"] = kmeans.labels_.astype(str)

sc.pl.spatial(sp_adata, color=['true_labels'],size=1,spot_size=1,save='_MousePDAC_kmeans2.pdf')

