from sklearn.metrics import adjusted_rand_score
from scipy.spatial.distance import jensenshannon

from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
# Panel- 
# 基于细胞比例进行聚类，看看与基于count 聚类的区别
# panel-4 umap of predicted cell-type perc
os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/MouseBrain2022/01_cell_perc/")
ct_order = ['Astrocytes', 'Microglia', 'Neurons', 'Oligodendrocytes']

sp_adata = sc.read_h5ad("01_data/MouseBrain2022_spot208_impute.h5ad")

sc.pp.normalize_total(sp_adata)
sc.pp.log1p(sp_adata)

kmeans = KMeans(n_clusters=2, random_state=42).fit(np.array(sp_adata.to_df()))
sp_adata.obs["true_labels"] = kmeans.labels_.astype(str)

methods=["SpatialDC_reconstruct","Tangram","Cell2location","Stereoscope","DestVI","CARD","spotlight","spatialdwls","Seurat"]

save_ARI = []
save_Purity = []
for m in methods:
    if(m =="SpatialDC_initial"):
        temp_pred_path = f"{project_dir}/SpatialDC/initial.csv"
    elif(m =="SpatialDC_reconstruct"):
        temp_pred_path = f"{project_dir}/SpatialDC/reconstruct.csv"
    else:
        temp_pred_path = f"{project_dir}/{m}/{m}.csv"

    ct_perc_pred = pd.read_csv(temp_pred_path,index_col=0)                    
        
    ct_perc_adata = ad.AnnData(X=csr_matrix(np.array(ct_perc_pred[ct_order]).clip(0)))    
    ct_perc_adata.obsm["spatial"] = sp_adata.obsm["spatial"]
#     sc.pp.filter_genes(ct_perc_adata, min_cells=1)
    kmeans = KMeans(n_clusters=2, random_state=42).fit(np.array(ct_perc_pred[ct_order]).clip(0))

    sp_adata.obs[f"{m}_labels"] = kmeans.labels_.astype(str)
  
    ari = adjusted_rand_score(sp_adata.obs['true_labels'].values, sp_adata.obs[f"{m}_labels"])    
   
    
    # 计算Purity值
    y_pred = sp_adata.obs[f"{m}_labels"].values.astype(int)
    y_true = sp_adata.obs['true_labels'].values.astype(int)
    N = len(y_pred)
    purity = 0
    for k in np.unique(y_pred):
        idx = (y_pred == k)
        labels = y_true[idx]
        count = np.bincount(labels)
        purity += np.max(count)
    purity /= N

    print(f'{m}_{round(ari,3)}_{round(purity,3)}')
    save_ARI.append(round(ari,3))
    save_Purity.append(round(purity,3))
#     sc.pl.spatial(sp_adata, color=[f"{m}_labels"],size=1,spot_size=1)    
#     # 计算Silhouette Coefficient
#     score = silhouette_score(np.array(ct_perc_pred[ct_order]),sp_adata.obs["leiden"].values, metric='euclidean')
#     print(f'{m}_Silhouette Coefficient:', score)

sc.pl.spatial(sp_adata, color=sp_adata.obs.columns[3:],size=1,title=methods,spot_size=1,palette=['#ff7f0e','#1f77b4'],save="comparision_kmean_cluster.pdf")    
