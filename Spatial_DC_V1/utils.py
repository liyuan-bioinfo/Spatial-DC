import sys
import os
import pandas as pd
import scanpy as sc
import numpy as np
from scipy.stats import pearsonr,spearmanr
from sklearn.metrics import mean_squared_error

import matplotlib.pyplot as plt

import torch
import torch.nn as nn
import torch.nn.functional as F

import warnings
warnings.filterwarnings("ignore")
import glob
from sklearn.preprocessing import scale,MinMaxScaler, StandardScaler
import anndata
from scipy.sparse import csr_matrix
from scipy.stats import norm
from tqdm import tqdm
import math
import random
import time

import torch_geometric
from torch_geometric.nn import VGAE
from torch_geometric.nn import VGAE, GCNConv
from torch_geometric.utils import negative_sampling, remove_self_loops, add_self_loops

# for higher resolutions    
from scipy.spatial.distance import cdist
from scipy.stats import median_abs_deviation
class scale_dataset(torch.utils.data.Dataset):

    def __init__(self, inputs, labels,cell_matrix):
        super(scale_dataset, self).__init__()

        # normalized total
        # target_sum = np.median(np.sum(inputs,axis=1))
        # total_counts = np.sum(inputs, axis=1)
        # inputs = inputs / total_counts[:, np.newaxis]
        # inputs = inputs * target_sum

        inputs = np.log1p(inputs) # log2 normalization                        

        inputs = scale(inputs,axis=1) #sample
        
        self.inputs = inputs

        self.cell_matrix = cell_matrix
        self.labels = labels

    def __getitem__(self, index):
        _input = self.inputs[index]
        _label = self.labels[index]        
        cell_matrix = self.cell_matrix[index]        

        return _input, _label,cell_matrix

    def __len__(self):
        return len(self.labels)   



class Utils():

    def generate_vectors(X, N, min_sum, max_sum):
        """
        Generates X series of vectors with a random sum within the range [min_sum, max_sum].
        Each vector has size N.
        """
        vectors = []
        for i in range(X):
            # initialize the vector with zeros
            vec = np.zeros(N, dtype=int)
            # generate a random sum within the range [min_sum, max_sum]
            total_sum = np.random.randint(min_sum, max_sum+1)
            # repeat until the sum of the vector is equal to total_sum
            while np.sum(vec) < total_sum:
                # generate a random index in the range [0, N-1]
                j = np.random.randint(0, N-1)
                # generate a random value between 1 and the remaining sum
                val = np.random.randint(1, total_sum-np.sum(vec)+1)
                # update the vector at the selected index
                vec[j] += val
            vectors.append(vec)
        return vectors    
    def generate_simulated_data(sc_data, cell_type_obs="celltype",samplenum=10):

        # print(sc_data)
        if isinstance(sc_data.X, np.ndarray):
            pass
        else:
            sc_data.X = sc_data.X.toarray() #73 * 10466
        
        sc_data = pd.DataFrame(sc_data.X, index=sc_data.obs[cell_type_obs], columns=sc_data.var.index)

        sc_data[cell_type_obs] = sc_data.index
        sc_data.index = range(len(sc_data))
        celltype_groups = sc_data.groupby(cell_type_obs).groups
        num_celltype = len(sc_data[cell_type_obs].value_counts())
        
        sc_data.drop(columns=cell_type_obs, inplace=True)
        
        sc_data = anndata.AnnData(sc_data)
        
        sc_data = sc_data.X
        sc_data = np.ascontiguousarray(sc_data, dtype=np.float32)

        for key, value in celltype_groups.items():
            celltype_groups[key] = np.array(value)
                
        # History
        prop = np.random.dirichlet(np.ones(num_celltype), samplenum)
        n = np.random.randint(10,100,size=prop.shape[0]).reshape(-1,1)
        cell_num = np.floor(n * prop) # 可以增加超参增加单个spot组成的细胞数量        
        prop = cell_num / np.sum(cell_num, axis=1).reshape(-1, 1)

        # update in 20240803
        # cell_num = Utils.generate_vectors(X=samplenum, N=num_celltype, min_sum=10, max_sum=100)
        # prop = cell_num / np.sum(cell_num, axis=1).reshape(-1, 1)


        sample = np.zeros((prop.shape[0], sc_data.shape[1])) # 1000 * 253
        allcellname = celltype_groups.keys() # dict_keys, 19
        


        cell_matrix_intensity = np.zeros((prop.shape[0], sc_data.shape[1] *  num_celltype))
        cell_matrix_weight = np.zeros((prop.shape[0], sc_data.shape[1] *  num_celltype))

        # add cell-type matrix
        for i, sample_prop in tqdm(enumerate(cell_num)): # 遍历每个spot 对应的cell perc
            
            
            temp_subspot_intensity = np.zeros((num_celltype, sc_data.shape[1])) # celltypes * proteins
            

            for j, cellname in enumerate(allcellname): # 遍历cell types
                select_index = np.random.choice(celltype_groups[cellname], size=int(sample_prop[j]), replace=True) # 从cell-type抽取该数量的cell
                celltype_intensity = sc_data[select_index].sum(axis=0)
                sample[i] += celltype_intensity# mixed
                
                
                temp_subspot_intensity[j] = celltype_intensity # intensity        
                    
            # print(temp_subspot_weight.sum(axis=0).shape)
            cell_matrix_intensity[i] = temp_subspot_intensity.T.reshape(1,sc_data.shape[1] *  num_celltype)


            temp_abundance = temp_subspot_intensity.sum(axis=0)
            temp_abundance = np.where(temp_abundance == 0, 1, temp_abundance) #防止 /0

            temp_scale = temp_subspot_intensity / temp_abundance #celltypes * proteins, scale to 1
            cell_matrix_weight[i] =  temp_scale.T.reshape(1,sc_data.shape[1] *  num_celltype)
        
        _spot_data = sample # spots * proteins
        _spot_type = prop # spots * celltypes
        ct_protein_intensity = cell_matrix_intensity
        ct_protein_weight = cell_matrix_weight
        return(_spot_data, _spot_type, ct_protein_intensity,ct_protein_weight)

    def construct_knn(sp_adata,k_graph):
        adata_X = sp_adata.X.toarray()        
        # adata_X = (adata_X - adata_X.mean(0)) / (adata_X.std(0) + 1e-10)
        adata_X = np.log1p(adata_X)
        adata_X = scale(adata_X,axis=1) # by sample
        # adata_X = scale(adata_X,axis=0)
        gene_mat = torch.Tensor(adata_X)
        del adata_X
        # if pca_dim > 0 and gene_mat.shape[1] > pca_dim:  # PCA
        #     u, s, v = torch.pca_lowrank(gene_mat, pca_dim)
        #     gene_mat = torch.matmul(gene_mat, v)

        # Construct spatial graph.        

        # gene_mat = torch.Tensor(adata.X / adata.X.sum(-1))
        cell_coo = torch.Tensor(sp_adata.obsm['spatial'])
        data = torch_geometric.data.Data(x=gene_mat, pos=cell_coo)

        data = torch_geometric.transforms.KNNGraph(k=k_graph, loop=True)(data)
        data = torch_geometric.transforms.Distance()(data)
        data.z = torch.Tensor(np.array(sp_adata.uns["protein_weight_pred"]))
        data.y = torch.Tensor(np.array(sp_adata.uns["cellperc_initial"]))
        # print(data.y)
        
        data.edge_weight = 1 - data.edge_attr[:,0]
        return(data)

    def prepare_dataloader(spot_data, spot_type,protein_weight,test_size,batch_size):
        # Set training and validation sets.
        num_spots = spot_data.shape[0]
        
        if test_size < 1:
            test_size = int(num_spots * test_size)
        
        _train_set = scale_dataset(spot_data[test_size:], spot_type[test_size:], protein_weight[test_size:])
        _valid_set = scale_dataset(spot_data[:test_size], spot_type[:test_size], protein_weight[:test_size])
        # print('  Train: %d spots,  Valid: %d spots' % (len(_train_set), len(_valid_set)))

        _train_loader = torch.utils.data.DataLoader(_train_set, batch_size=batch_size, shuffle=True, pin_memory=True)
        _valid_loader = torch.utils.data.DataLoader(_valid_set, batch_size=batch_size * 2, shuffle=False, pin_memory=True)
        return _train_loader, _valid_loader

    # 给定坐标轴[x,y], 产生对应个数的坐标
    def getHighresCords(Cords, numCell):
        ED = cdist(Cords, Cords)
        n = Cords.shape[0]
        dis = np.zeros(n)
        for i in range(n):
            dis[i] = np.min(np.delete(ED[i], i))
        min_distance = np.median(dis) / 2
        np.random.seed(20240201)
        Cords_new = []
        for i in range(n):

            def getPointsWithinSquare(Cords, i, numCell, min_distance):
                minX = min_distance
                minY = min_distance
                rectangular_x = np.random.uniform(low=Cords[i, 0] - minX, high=Cords[i, 0] + minX, size=int(numCell[i]))
                rectangular_y = np.random.uniform(low=Cords[i, 1] - minY, high=Cords[i, 1] + minY, size=int(numCell[i]))#change this
                df = np.column_stack((rectangular_x, rectangular_y))
                return df

            df = getPointsWithinSquare(Cords, i, numCell, min_distance)
            df = np.column_stack((df, np.repeat(f"{Cords[i, 0]}x{Cords[i, 1]}", df.shape[0])))
            df = np.column_stack((df, np.repeat(Cords[i, 0], df.shape[0])))
            df = np.column_stack((df, np.repeat(Cords[i, 1], df.shape[0])))
            Cords_new.append(df)
        Cords_new = np.concatenate(Cords_new, axis=0)        
        return Cords_new        