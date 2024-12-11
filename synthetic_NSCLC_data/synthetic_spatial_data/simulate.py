import time
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl 
import matplotlib
from scipy.sparse import csr_matrix
import anndata as anndata

os.chdir("NSCLC2019/01_cell_perc/01_data/simulations")

# This grid method extended from (Li, et al. Nat.Methods. 2022.) [https://www.nature.com/articles/s41592-022-01480-9]
def Simulated(adata, CoordinateXlable, CoordinateYlable, window):
    if os.path.exists(save_dir):
        print ('The output file is in ' + save_dir)
    else:
        os.mkdir(save_dir)

    spatial_rna = adata.to_df()
    spatial_meta = adata.obs
    spatial_loc = adata.obsm.to_df()
    
    combined_spot = []
    combined_spot_loc = []
    window=window
    c = 0
    for x in np.arange((spatial_loc[CoordinateXlable].min()//window),spatial_loc[CoordinateXlable].max()//window+1):
        for y in np.arange((spatial_loc[CoordinateYlable].min()//window),spatial_loc[CoordinateYlable].max()//window+1):
            tmp_loc = spatial_loc[(x*window < spatial_loc[CoordinateXlable]) & (spatial_loc[CoordinateXlable] < (x+1)*window) & (y*window < spatial_loc[CoordinateYlable]) & (spatial_loc[CoordinateYlable] < (y+1)*window)]
            if len(tmp_loc) > 0:
                c += 1
                combined_spot_loc.append([x,y])
                combined_spot.append(tmp_loc.index.to_list())
            
    combined_cell_counts = pd.DataFrame([len(s) for s in combined_spot],columns=['cell_count'])    
    print ('The simulated spot has cells with ' + str(combined_cell_counts.min()[0]) + ' to ' + str(combined_cell_counts.max()[0]))
    combined_spot_loc = pd.DataFrame(combined_spot_loc, columns=['x','y'])

    combined_spot_exp = []
    combined_cell_type_exp = []

    for s in combined_spot:
        spot_exp = spatial_rna.loc[s,:].sum(axis=0).values
        combined_spot_exp.append(spot_exp)

        spot_exp_celltype_meta = spatial_meta.loc[s,:] # obtain the index of selected spots
        cell_type_exp = []
        for cell_type in np.unique(spatial_meta['celltype']): # obtain the index of celltypes from selcted spots [spot => cell types]
            spot_celltype_exp = spatial_rna.loc[s,:].loc[spot_exp_celltype_meta["celltype"] == cell_type].sum(axis=0).values
            cell_type_exp.append(spot_celltype_exp)
        combined_cell_type_exp.append(cell_type_exp)

    combined_spot_exp = pd.DataFrame(combined_spot_exp, columns=spatial_rna.columns)
    combined_cell_type_exp = pd.DataFrame(combined_cell_type_exp, columns=[f'{cell_type}' for cell_type in np.unique(spatial_meta['celltype'])])

    combined_spot_clusters = pd.DataFrame(np.zeros((len(combined_spot_loc.index),len(np.unique(spatial_meta['celltype'])))),columns=np.unique(spatial_meta['celltype']))
    for i,c in enumerate(combined_spot):
        for clt in spatial_meta.loc[c,'celltype']:
            combined_spot_clusters.loc[i,clt] += 1

    print ('The simulated spot has size ' + str(combined_spot_clusters.shape[0]))
    return(combined_spot_exp, combined_spot_loc, combined_spot_clusters, combined_cell_counts,  combined_cell_type_exp)

save_dir = "synthetic_cellnum_noise_Lung5_add_gd_scp2019"
dataset_name = "Simu"

def add_gaussian_noise_by_all_v2(matrix, noise_ratio):
    log1p_matrix = np.log1p(matrix)        
    matrix_sd = np.std(log1p_matrix)        
    log1p_noise_matrix = np.random.normal(0, matrix_sd * noise_ratio, matrix.shape)  

    output_matrix = (np.exp(log1p_noise_matrix + log1p_matrix)-1).clip(0)
    return(output_matrix)

for seed in [0]:
    for cells in [10,15,20]: # varying cell densities

        sp_adata = sc.read_h5ad(f'NSCLC_ct3_scp_Lung5_R1_impute.h5ad')
        
        combined_spot_exp, combined_spot_loc, combined_spot_clusters, combined_cell_counts, combined_cell_type_exp = Simulated(adata=sp_adata,CoordinateXlable="spatial1", 
                                                                                    CoordinateYlable="spatial2", window=cells) # generate spatial proteomics without external noise

        for noise in [0,0.25,0.5,0.75,1.0]:  # varying external Gaussian noise
        
            combined_spot_exp2 = add_gaussian_noise_by_all_v2(combined_spot_exp, noise_ratio=noise)
            temp_sp_adata = ad.AnnData(csr_matrix(combined_spot_exp2))
            temp_sp_adata.var_names = combined_spot_exp2.columns.values # add protein name

            combined_spot_clusters = (combined_spot_clusters.T/combined_spot_clusters.sum(axis=1)).T

            random_spots = temp_sp_adata.shape[0]
            simu_index = ["spot_"+str(i) for i in list(range(1,random_spots+1))] # set the simulated spot index
            combined_spot_clusters.index = simu_index
            temp_sp_adata.obs.index = simu_index            

            temp_sp_adata.uns["celltype_gd"] = combined_spot_clusters
            temp_sp_adata.uns["combined_cell_counts"] = combined_cell_counts
            temp_sp_adata.obsm["spatial"]=combined_spot_loc.values
                        
            temp_sp_adata.write_h5ad(f"{save_dir}/{dataset_name}_seed{seed}_cells{str(cells)}_noise{int(noise*100)}.h5ad") # generate synthetic spatial proteomics with external varying Gaussian noise
            temp_sp_adata.uns['celltype_gd'].to_csv(f"{save_dir}/{dataset_name}_seed{seed}_cells{str(cells)}_noise{int(noise*100)}_gd.csv") # corresponding ground truth of cell-type composition of each spot
            temp_sp_adata.uns['combined_cell_counts'].to_csv(f"{save_dir}/{dataset_name}_seed{seed}_cells{str(cells)}_noise{int(noise*100)}_combined_cell_counts.csv") # summed cell counts of each spot

            temp_spatial_df = pd.DataFrame(temp_sp_adata.obsm["spatial"],index=temp_sp_adata.obs.index,columns=["x","y"])
            temp_spatial_df.to_csv(f"{save_dir}/{dataset_name}_seed{seed}_cells{str(cells)}_noise{int(noise*100)}_loc.csv") # syntheitc spots' spatial coordinates
            temp_sp_adata.to_df().to_csv(f"{save_dir}/{dataset_name}_seed{seed}_cells{str(cells)}_noise{int(noise*100)}_expression.csv") # synthetic spots' proteomic profiles

            # obtain the ground truth for assessment of reconstruction proteomic profiles
            combined_cell_type_exp.index = simu_index
            output_combined_cell_type_exp = []
            output_combined_cell_type_scaled_exp = []
            output_combined_cell_type_gd = []
            output_combined_cell_type_exp_index = []
            output_combined_cell_type_exp_index_celltype = []
            for spot in combined_cell_type_exp.index:
                for ct in combined_cell_type_exp.columns:
                    spot_exp = combined_cell_type_exp.loc[spot,ct]
                    spot_gd = temp_sp_adata.uns['celltype_gd'].loc[spot,ct]
                    spot_scaled_exp = spot_exp / (spot_gd + 1e-6)
                    
                    output_combined_cell_type_gd.append(spot_gd.tolist())
                    output_combined_cell_type_exp.append(spot_exp.tolist())
                    output_combined_cell_type_scaled_exp.append(spot_scaled_exp.tolist())
                    output_combined_cell_type_exp_index.append(f"{spot}_{ct}")
                    output_combined_cell_type_exp_index_celltype.append(f"{ct}")

            combined_cell_type_adata = anndata.AnnData(X=csr_matrix(np.array(output_combined_cell_type_exp)))
            combined_cell_type_adata.obs.index = output_combined_cell_type_exp_index
            combined_cell_type_adata.obs["cellperc"] = output_combined_cell_type_gd
            combined_cell_type_adata.obs["celltype"] = output_combined_cell_type_exp_index_celltype
            combined_cell_type_adata.var = temp_sp_adata.var
            combined_cell_type_adata.write_h5ad(f"{save_dir}/{dataset_name}_seed{seed}_cells{str(cells)}_noise{int(noise*100)}_combined_celltype_exp.h5ad") # obtain cell-type proteomic profiles of each spot, served as ground truth for reconstructed proteomics profiles
            
            combined_cell_type_scaled_adata = combined_cell_type_adata.copy()
            combined_cell_type_scaled_adata.X = csr_matrix(np.array(output_combined_cell_type_scaled_exp))
            combined_cell_type_scaled_adata.write_h5ad(f"{save_dir}/{dataset_name}_seed{seed}_cells{str(cells)}_noise{int(noise*100)}_combined_celltype_scaled_exp.h5ad") # scale the cell-type proteomic profiles of each spot from ground truth, which used for statistical analysis
