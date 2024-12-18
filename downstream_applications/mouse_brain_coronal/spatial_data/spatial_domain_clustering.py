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
