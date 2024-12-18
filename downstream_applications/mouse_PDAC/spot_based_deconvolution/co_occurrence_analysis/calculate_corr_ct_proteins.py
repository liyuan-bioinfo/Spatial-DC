# obtain the corr between proteins and ct perc
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

plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# --------------------------------------------------------
# Calculate the Pearson Corr. Coef between predicted cell-type composition and original cell-type proteomic profiles
# for Fig. 5e-g
os.chdir("")

sp_adata = sc.read_h5ad("00_raw/MousePDAC2023_impute_Spots108_protein3607.h5ad")
sc.pp.normalize_total(sp_adata)
sc.pp.log1p(sp_adata)

df = sp_adata.var
df['pids'] = df.apply(lambda row: str(row['gene']) +"_"+ str(row['pid']), axis=1)
sp_adata.var.index = sp_adata.var["pids"]
gene_expressions = sp_adata.to_df()

ct_perc_pred = pd.read_csv("01_cell_perc/03_benchmark_methods/exp_ct10_impute_v2/SpatialDC/reconstruct.csv",index_col=0)

high_corr_genes_df = pd.DataFrame(columns=['CellType', 'Protein', 'Corr'])
for ct in ct_order:    
    correlations = gene_expressions.corrwith(ct_perc_pred[ct],method='pearson')      
    highly_correlated_genes = correlations  

    for gene in highly_correlated_genes.index.tolist():        
        high_corr_genes_df = pd.concat([high_corr_genes_df, pd.DataFrame([{'CellType': ct, 'Corr': highly_correlated_genes[gene], 'Protein': gene}])], ignore_index=True)

high_corr_genes_df.to_csv("01_cell_perc/analysis/spatial_corr_proteins_SpatialDC_v20241127.csv")

# conduct analysis and visualization with R scripts
