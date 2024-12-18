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

import warnings
warnings.filterwarnings("ignore")
# ---------------------------------------------------
class scale_dataset(torch.utils.data.Dataset):

    def __init__(self, inputs, labels,cell_matrix):
        super(scale_dataset, self).__init__()

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
  
    def construct_knn(sp_adata,k_graph):
        adata_X = sp_adata.X.toarray()        
        # adata_X = (adata_X - adata_X.mean(0)) / (adata_X.std(0) + 1e-10)
        adata_X = np.log1p(adata_X)
        adata_X = scale(adata_X,axis=1) # by sample
        # adata_X = scale(adata_X,axis=0)
        gene_mat = torch.Tensor(adata_X)
        del adata_X

        # Construct spatial graph.        
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
   
