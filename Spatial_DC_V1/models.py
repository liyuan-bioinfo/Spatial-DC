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

random.seed(0)
os.environ['PYTHONHASHSEED'] = '0'
np.random.seed(0)
torch.manual_seed(0)
torch.cuda.manual_seed(0)
torch.cuda.manual_seed_all(0)
torch.backends.cudnn.benchmark = False
torch.backends.cudnn.deterministic = True
torch.use_deterministic_algorithms(True)
os.environ['CUBLAS_WORKSPACE_CONFIG']=':4096:8'

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

from .utils import Utils
import anndata as ad

# Use DNN to train a distriubtion model, which initial predict cell-type abundance of each spot
class DNN(nn.Module):
    state = {} # save model params

    def __init__(self,protein_num, ct_num, drop_rate=0.5):
        # print("initial Distribution Model")
        super(DNN, self).__init__()        

        self.protein_num = protein_num
        self.ct_num = ct_num

        self.encoder = nn.Sequential(
            # Layer 1
            nn.Dropout(drop_rate),            
            nn.Linear(self.protein_num, self.protein_num // 2),    
            nn.CELU(),

            # Layer 2
            nn.Dropout(drop_rate),            
            nn.Linear(self.protein_num // 2, self.protein_num // 4),    
            nn.CELU()

        ) 

        self.classifier = nn.Sequential(            
                                 
            nn.Linear(self.protein_num // 4, self.ct_num),    
            nn.ReLU()            
        )        

        self.deconv = nn.Sequential(            
                        
            nn.Linear(self.protein_num // 4, self.protein_num * self.ct_num),                
            nn.ReLU()
        )        


    def forward(self, x):
        z_cat = self.encoder(x)

        # cell type perc
        ct_perc = self.classifier(z_cat)
        temp_sum = ct_perc.sum(axis=1) + 1e-6
        # temp_sum = torch.where(temp_sum == 0, 1, temp_sum) # divide zero
        ct_perc = (ct_perc.T/temp_sum).T        

        # protein weight
        protein_weight = self.deconv(z_cat)
        protein_weight = protein_weight.reshape(z_cat.shape[0]*self.protein_num,self.ct_num) # scale
        
        temp_sum = protein_weight.sum(axis=1) + 1e-6
        # temp_sum = torch.where(temp_sum == 0, 1, temp_sum) # divide zero
        protein_weight = (protein_weight.T/temp_sum).T
        protein_weight = protein_weight.reshape(z_cat.shape[0], self.protein_num * self.ct_num)
        
        return ct_perc, protein_weight
    


class GCN(nn.Module):

    state = {}
    def __init__(self, protein_num, ct_num):
        super(GCN, self).__init__()
        self.protein_num = protein_num
        self.ct_num = ct_num
        self.gae_dim, self.dae_dim, self.feat_dim = [32, 20], [100, 20], 64
        self.fcat_dim = self.dae_dim[1] + self.gae_dim[1]
        self.drop_rate = 0.5

        self.encoder = nn.Sequential(            
            nn.Dropout(self.drop_rate),
            nn.Linear(self.protein_num, self.dae_dim[0]),
            nn.CELU(),                 

            nn.Dropout(self.drop_rate),              
            nn.Linear(self.dae_dim[0], self.dae_dim[1]),
            nn.CELU(),
            
        )
        self.decoder = nn.Linear(self.fcat_dim, self.protein_num)
        self.vgae = VGAE(GraphEncoder(self.dae_dim[1], self.gae_dim[0], self.gae_dim[1]))
        
        self.feat_fc_g = nn.Sequential(
            # nn.Dropout(self.drop_rate),
            nn.Linear(self.fcat_dim, self.feat_dim),
            nn.ReLU()
        )
        # self.classifier = nn.Linear(self.fcat_dim, self.ct_num)
        self.classifier = nn.Sequential(
            # nn.Dropout(self.drop_rate),
            nn.Linear(self.fcat_dim, self.ct_num),
            nn.ReLU()
        )

        self.deconv = nn.Sequential(            
            
            # nn.Dropout(self.drop_rate),            
            nn.Linear(self.fcat_dim,self.protein_num * self.ct_num),
            nn.ReLU()
        )


    def forward(self, x, edge_index, edge_weight):

        feat_x = self.encoder(x) #[5577, 254] to [5577, 20]        
        feat_g = self.vgae.encode(feat_x, edge_index, edge_weight) #[5577, 20] to [5577, 20], add spatial info
        feat = torch.cat([feat_x, feat_g], 1) # [5577, 40]
        
        # task-1 cell-abundance prediction
        ct_perc = self.classifier(feat) # [5577, 40] -> [5577, 19]                        
        temp_sum = ct_perc.sum(axis=1) + 1e-6
        # temp_sum = torch.where(temp_sum == 0, 1, temp_sum)
        ct_perc = (ct_perc.T/temp_sum).T                

        # task-2, protein weight
        protein_weight = self.deconv(feat)
        protein_weight = protein_weight.reshape(feat.shape[0]*self.protein_num,self.ct_num)#protein spot * celltype

        temp_sum = protein_weight.sum(axis=1) + 1e-6
        # temp_sum = torch.where(temp_sum == 0, 1, temp_sum)
        protein_weight = (protein_weight.T/temp_sum).T
        protein_weight = protein_weight.reshape(feat.shape[0], self.protein_num * self.ct_num)#protein spot * celltype

        # feature encoder-decoder matrix
        x_dec = self.decoder(feat) # [5577, 64] -> [5577, 254]
        
        feat_g = self.feat_fc_g(feat) # [5577, 40] -> [5577, 64]
        gae_loss = self.recon_loss(feat_g, edge_weight, edge_index) + (1 / len(x)) * self.vgae.kl_loss()
        
        return ct_perc, x_dec, gae_loss,protein_weight

    # 计算GAE Loss
    def recon_loss(self, z, edge_weight, pos_edge_index, neg_edge_index=None): # edge weight
        pos_dec = self.vgae.decoder(z, pos_edge_index, sigmoid=True)
        

        # pos_loss = F.binary_cross_entropy_with_logits(pos_dec, edge_weight)

        pos_loss = nn.MSELoss()(pos_dec, edge_weight)
        # Make sure self-loops in positive samples.
        pos_edge_index, _ = remove_self_loops(pos_edge_index)
        pos_edge_index, _ = add_self_loops(pos_edge_index)
        if neg_edge_index is None:
            neg_edge_index = negative_sampling(pos_edge_index, z.size(0))

        neg_dec = self.vgae.decoder(z, neg_edge_index, sigmoid=True)
        # neg_loss = nn.Sigmoid()(neg_dec).mean()
        neg_loss = nn.MSELoss()(neg_dec, torch.zeros_like(edge_weight))
        # neg_loss = -F.logsigmoid(-neg_dec).mean()
        return pos_loss + neg_loss    
    
class GraphEncoder(nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels):
        super(GraphEncoder, self).__init__()
        self.gc_feat = GCNConv(in_channels, hidden_channels)
        self.gc_mean = GCNConv(hidden_channels, out_channels)
        self.gc_logstd = GCNConv(hidden_channels, out_channels)

    def forward(self, x, edge_index, edge_weight):
        x = self.gc_feat(x, edge_index, edge_weight).relu()
        # print(x.shape) # [5615, 32]
        mean = self.gc_mean(x, edge_index, edge_weight)
        logstd = self.gc_logstd(x, edge_index, edge_weight)

        return mean, logstd

class SpatialDC:

    def __init__(self,sc_adata,sp_adata,print_info=True,ct_obs="celltype"):

        self.celltype_obs = ct_obs

        intersect = np.intersect1d(sc_adata.var_names, sp_adata.var_names)
        self.sc_adata = sc_adata[:, intersect].copy()
        self.sp_adata = sp_adata[:, intersect].copy()
        self.obs_index = sp_adata.obs.index.values    
        self.obs_num = len(self.obs_index)

        self.protein_num = self.sc_adata.shape[1]
        self.celltype_num = len(self.sc_adata.obs[self.celltype_obs].cat.categories.tolist())
        self.celltype_order = self.sc_adata.obs[self.celltype_obs].cat.categories.tolist()

        # use for DNN model
        self.distriubtion_model = None

        self.print_info = print_info


    def __str__(self):
        return f"Intersected sc_adata shape: {self.sc_adata.shape}\n" \
            f"Intersected sp_adata shape: {self.sp_adata.shape}\n" \
            f"obs_num: {self.obs_num}\n" \
            f"protein_num: {self.protein_num}\n" \
            f"celltype_num: {self.celltype_num}\n" \
            f"celltype_order: {self.celltype_order}"

    def setup_distribution_model(self, test_size=0.1, batch_size=128, spot_num=10000, epochs=10, lr=0.001):

        # batch_size: 128 => lr 0.001
        # batch_size = 128
        # factor = int(batch_size / 128) - 1
        # lr = 0.001 * (np.sqrt(2))**factor        

        # prepare simulated data                
        spot_data, spot_type,_,protein_weight = Utils.generate_simulated_data(sc_data=self.sc_adata,cell_type_obs=self.celltype_obs, 
                                                                    samplenum=spot_num)        
        
        self.train_loader, self.valid_loader = Utils.prepare_dataloader(spot_data=spot_data, spot_type=spot_type,protein_weight=protein_weight,test_size=test_size,batch_size=batch_size) # create dataload of sc        

        # register model object
        self.distriubtion_model = DNN(protein_num=self.protein_num,ct_num=self.celltype_num).to(device)
        self.distriubtion_model.state = {"lr":lr,"lr_decay":0.95,"lr_step":1,"weight_decay":1e-5,"start_epoch":1,"epochs":epochs,"spot_num":spot_num}


    def train_distribution_model(self):

        # set optimizer for weight_decay and scheduler for lr_decay
        lr = self.distriubtion_model.state["lr"]
        lr_step = self.distriubtion_model.state["lr_step"]
        lr_decay = self.distriubtion_model.state["lr_decay"]
        weight_decay = self.distriubtion_model.state["weight_decay"]
        start_epoch = self.distriubtion_model.state["start_epoch"]
        epochs = self.distriubtion_model.state["epochs"]


        params = []
        for key, value in dict(self.distriubtion_model.named_parameters()).items():
            if 'bias' in key:
                params.append({'params': [value], 'lr': lr, 'weight_decay': 0})
            else:
                params.append({'params': [value], 'lr': lr, 'weight_decay': weight_decay,'betas':[0.9,0.999]})
            
        optimizer = torch.optim.Adam(params)
        
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, lr_step, gamma=lr_decay)

        TRAIN_LOSS = []
        TRAIN_X_RECON_LOSS = []
        TRAIN_X_PROP_LOSS = []

        VALID_LOSS = []
        VALID_X_RECON_LOSS = []
        VALID_X_PROP_LOSS = []

        VALID_METRIC = []

        criterion = torch.nn.L1Loss()
        # criterion = torch.nn.MSELoss()#
        for epoch in range(start_epoch, epochs + 1):

            # ************** Train ************** #
            if(self.print_info):
                print('[Epoch %2d][Lr  %.e]' % (epoch, scheduler.get_last_lr()[0]))
            self.distriubtion_model.train()
            train_loss = 0
            train_x_recon_loss = 0
            train_cell_prop = 0

            # Train model.
            for batch_idx, (protein_intensity, ct_perc, protein_weight) in enumerate(self.train_loader):
                                
                protein_intensity  = protein_intensity.float().to(device)
                ct_perc = ct_perc.float().to(device)
                protein_weight = protein_weight.float().to(device)
                

                # Forward propagation.
                with torch.cuda.amp.autocast():
                
                    pred_cell_perc, pred_protein_weight = self.distriubtion_model(protein_intensity)                    

                    loss = criterion(pred_cell_perc, ct_perc) + criterion(pred_protein_weight, protein_weight)                    
                                        
                train_loss += loss.item() 

                # if math.isnan(train_loss):
                #     epochs = epoch-1                                        
                #     break

                # Backward propagation.
                optimizer.zero_grad()
                loss.backward()
                
                optimizer.step()

                if(self.print_info):
                    print('[Epoch %2d][Iter %4d/ %4d] Loss: %.5f \r'
                        % (epoch, batch_idx + 1, len(self.train_loader), train_loss / (batch_idx + 1)), end='')


            # if math.isnan(train_loss):
            #     epochs = epoch-1
            #     break
            # Update global variables.
            TRAIN_LOSS.append(train_loss / len(self.train_loader))            

            # ************** Valid ************** #
            self.distriubtion_model.eval()
            valid_loss = 0

            # Validate model.
            with torch.no_grad():
                for batch_idx, (protein_intensity, cell_perc,protein_weight) in enumerate(self.valid_loader):
                    protein_intensity = protein_intensity.float().to(device)
                    cell_perc = cell_perc.float().to(device)
                    protein_weight = protein_weight.float().to(device)
                    
                    pred_cell_perc,pred_protein_weight = self.distriubtion_model(protein_intensity)
                    
                    valid_loss = criterion(cell_perc, pred_cell_perc) + criterion(protein_weight, pred_protein_weight)
                    
            
            VALID_LOSS.append(valid_loss / len(self.valid_loader))


            if(self.print_info):
                print('[Epoch %2d][Train] Loss: %.5f [Valid] Loss: %.5f ' % (epoch, TRAIN_LOSS[-1], VALID_LOSS[-1]))

            self.distriubtion_model.state.update({
                'model': self.distriubtion_model.state_dict(),
                'optimizer': optimizer.state_dict(),
                'epoch': epoch,
                'train_loss': TRAIN_LOSS,
                'valid_loss': VALID_LOSS,
                'lr': scheduler.get_last_lr()[0],
            })
            scheduler.step()


    def save_distribution_model(self, save_model_path):


        torch.save(self.distriubtion_model.state, save_model_path)
                
        # train_loss = self.distriubtion_model.state["train_loss"]
        # valid_loss = self.distriubtion_model.state["valid_loss"]
        # epochs = self.distriubtion_model.state["epochs"]
        # spot_num = self.distriubtion_model.state["spot_num"]
        
        # plt.plot(np.arange(0,epochs), train_loss)
        # plt.xlabel('Epochs')
        # plt.ylabel('Loss')
        # plt.title(f'Simu_Spot:[{spot_num}] | Train Loss:[{"%.4f" % train_loss[-1]}] | Valid Loss:[{"%.4f" % valid_loss[-1]}]')
        # plt.savefig(f"{save_model_path}.png")
        # plt.close()        

    def load_distribution_model(self, load_model_path):
        self.distriubtion_model = DNN(protein_num=self.protein_num,ct_num=self.celltype_num).to(device)

        checkpoint = torch.load(load_model_path, map_location=torch.device('cpu'))
        self.distriubtion_model.load_state_dict(checkpoint['model'])        
        
    # transfer trained model to spatial proteomics data
    # Prepare simulated data sets
    # output: cell perc and protein weight
    def transfer_distribution_model(self):

        self.distriubtion_model.eval()
        inputs = self.sp_adata.X.toarray() # 2492 * 205
        
        # inputs = scale(inputs,axis=1) # sample
        # Normalized total
        # target_sum = np.median(np.sum(inputs,axis=1))
        # total_counts = np.sum(inputs, axis=1)
        # inputs = inputs / total_counts[:, np.newaxis]
        # inputs = inputs * target_sum        

        inputs = np.log1p(inputs)
        inputs = scale(inputs,axis=1) #sample
        
        inputs = torch.from_numpy(inputs)

        inputs  = inputs.float().to(device)
        ct_perc,protein_weight = self.distriubtion_model(inputs)
        ct_perc = ct_perc.detach().cpu().numpy()
        protein_weight = protein_weight.detach().cpu().numpy()

        results = [] # meta_df

        for i in self.sp_adata.var_names:
            for ct in self.celltype_order:
                temp_ind = i + "_" + ct
                temp_ct = ct
                temp_p = i
                result = {"ID":temp_ind,"celltype":ct,"protein":i}
                results.append(result)
        results = pd.DataFrame(results)
        results.index = results["ID"]        
        
        self.sp_adata.uns["cellperc_initial"] = pd.DataFrame(ct_perc,columns=self.celltype_order,index=self.obs_index)
        self.sp_adata.uns["protein_weight_pred"] = pd.DataFrame(protein_weight.T,index=results.index,columns=self.obs_index)

        return(self.sp_adata)
        

    # reconstruct initial predicted cell-type abundance when spatial coordinate available(x, y)
    def setup_reconstruction_model(self, k_graph=30, epochs=100, w_cls=30, w_rec=30):

        knn_data = Utils.construct_knn(sp_adata=self.sp_adata,k_graph=k_graph)
        self.reconstruction_model = GCN(self.protein_num, self.celltype_num).to(device)
        self.reconstruction_model.state = {"knn_data":knn_data,"lr":1e-3,"lr_decay":1,"lr_step":1,"w_cls":w_cls,"w_dae":1,"w_gae":1,"w_rec":w_rec,"kd_T":1,"epochs":epochs,"weight_decay":1e-5}

    def reconstruct(self):
        knn_data = self.reconstruction_model.state["knn_data"]
        lr = self.reconstruction_model.state["lr"]
        lr_decay = self.reconstruction_model.state["lr_decay"]
        lr_step = self.reconstruction_model.state["lr_step"]
        w_cls = self.reconstruction_model.state["w_cls"]
        w_dae = self.reconstruction_model.state["w_dae"]
        w_gae = self.reconstruction_model.state["w_gae"]
        w_rec = self.reconstruction_model.state["w_rec"]

        kd_T = self.reconstruction_model.state["kd_T"]
        epochs = self.reconstruction_model.state["epochs"]                
        weight_decay = self.reconstruction_model.state["weight_decay"]                

        # prepare data for training
        knn_data = knn_data.to(device, non_blocking=True)        
        inputs, targets, protein_weight = knn_data.x, knn_data.y, knn_data.z.T
        edge_index  = knn_data.edge_index
        edge_weight = knn_data.edge_weight

        criterion = nn.L1Loss()

        params = []
        for key, value in dict(self.reconstruction_model.named_parameters()).items():
            if 'bias' in key:
                params.append({'params': [value], 'lr': lr, 'weight_decay': 0})
            else:
                params.append({'params': [value], 'lr': lr, 'weight_decay': weight_decay,'betas':[0.9,0.999]})
            
        optimizer = torch.optim.Adam(params)

        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, lr_step, gamma=lr_decay)        
        
        self.reconstruction_model.train()
        
        TRAIN_LOSS = []
        for epoch in range(1, epochs + 1):
            # Forward propagation.
            with torch.cuda.amp.autocast():

                outputs, x_dec, gae_loss,pred_protein_weight = self.reconstruction_model(inputs, edge_index, edge_weight)
                dae_loss = F.mse_loss(x_dec, inputs) # MSELoss, 内部固定写. 
                                
                cls_loss = criterion(outputs, targets)

                rec_loss = criterion(pred_protein_weight, protein_weight)
                
                loss = w_cls * cls_loss + w_rec * rec_loss + w_dae * dae_loss + w_gae * gae_loss 
            
            train_loss = loss.item()

            # Backward propagation.
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            scheduler.step()

            TRAIN_LOSS.append(train_loss)
            if(self.print_info):            
                print('  [Epoch %3d] Loss: %.5f' % (epoch, train_loss))
            
        self.reconstruction_model.state.update({"train_loss":TRAIN_LOSS})            

        # reconstruct
        self.reconstruction_model.eval()
        ct_perc, _, _,protein_weight = self.reconstruction_model(inputs, edge_index, edge_weight)

        ct_perc = ct_perc.detach().cpu().numpy()
        
        self.sp_adata.uns["cellperc_reconstruct"] = pd.DataFrame(ct_perc,columns=self.celltype_order,index=self.obs_index)

        protein_weight = protein_weight.detach().cpu().numpy()        
        re_cell_metrix = pd.DataFrame(protein_weight.T, index=self.sp_adata.uns["protein_weight_pred"].index, columns=self.sp_adata.uns["protein_weight_pred"].columns)
        self.sp_adata.uns["protein_weight_reconstruct"] = re_cell_metrix

        return(self.sp_adata)

    def save_reconstruction_model(self, save_model_path):

        torch.save(self.reconstruction_model.state, save_model_path)
                
        # train_loss = self.reconstruction_model.state["train_loss"]
        # epochs = self.reconstruction_model.state["epochs"]
                
        # plt.plot(np.arange(0,epochs), train_loss)
        # plt.xlabel('Epochs')
        # plt.ylabel('Loss')
        # plt.title(f'Train Loss:[{"%.4f" % train_loss[-1]}]')
        # plt.savefig(f"{save_model_path}.png")
        # plt.close()     

    # extended function-1
    # purify cell-type expression profiles from each spot
    # downstream tasks: novel celltype markers; far or near distance cells; protein protein interactions
    def purify_spots(self, norm=True, spatial_info=False):        
        sp_purified_adata = self.sp_adata.copy()

        cellperc_pred = sp_purified_adata.uns["cellperc_reconstruct"]
        protein_weight_pred = sp_purified_adata.uns["protein_weight_reconstruct"]        
        # cellperc_pred = sp_purified_adata.uns["cellperc_initial"]
        # protein_weight_pred = sp_purified_adata.uns["protein_weight_pred"]        

        
        # before normalized
        purified_weight_arr = protein_weight_pred.T # spot * nfeatures, protein weight        
        purified_intensity_arr = np.repeat(np.array(sp_purified_adata.to_df()), self.celltype_num, axis=1)
        purified_intensity_arr = np.multiply(purified_weight_arr, purified_intensity_arr) # spot * nfeatures, protein intensity
    

        # # normalized with cell perc
        # cellperc_pred_reshape = np.array(cellperc_pred).reshape(spot_num * celltype_num,1) #nspot * features
        if(norm):
            cellperc_pred_reshape = np.where(cellperc_pred == 0, 1e20, cellperc_pred) #防止 /0
            cellperc_pred_reshape = np.tile(np.array(cellperc_pred_reshape), self.protein_num)
            purified_intensity_arr = purified_intensity_arr / cellperc_pred_reshape

        # array to h5ad
        results = [] # meta_df
        for i in sp_purified_adata.var_names:
            for ct in self.celltype_order:
                temp_ind = i + "_" + ct
                temp_ct = ct
                temp_p = i
                result = {"ID":temp_ind,"celltype":ct,"protein":i}
                results.append(result)
        results = pd.DataFrame(results)
        results.index = results["ID"]

        purified_ct_adata = sc.AnnData(X=csr_matrix(purified_intensity_arr))
        purified_ct_adata.var = results
        purified_ct_adata.obs = sp_purified_adata.obs

        # print(purified_ct_adata.to_df())
        results = []
        for ct in self.celltype_order:
            temp_ind = (purified_ct_adata.var["ID"][purified_ct_adata.var["celltype"] == ct]).index.values
            results = np.append(results,temp_ind)

        purified_ct_adata = purified_ct_adata[:,results]
        # print(results)
        purified_ct_adata_reshape = np.array(purified_ct_adata.to_df()).reshape(self.obs_num * self.celltype_num, self.protein_num)

        results = [] # meta_df
        if(spatial_info):
            
            spatial_df = pd.DataFrame(sp_purified_adata.obsm["spatial"], index=sp_purified_adata.obs.index, columns=["X","Y"])
            for i in sp_purified_adata.obs.index:
                for ct in self.celltype_order:

                    temp_ind = i + "_" + ct
                    temp_ct = ct

                    temp_x = spatial_df["X"][i]
                    temp_y = spatial_df["Y"][i]
                    result = {"ID":temp_ind,"celltype":ct,"X":temp_x,"Y":temp_y}
                    results.append(result)
            results = pd.DataFrame(results)
            results.index = results["ID"]
        else:
            for i in sp_purified_adata.obs.index:
                for ct in self.celltype_order:

                    temp_ind = i + "_" + ct
                    temp_ct = ct

                    result = {"ID":temp_ind,"celltype":ct}
                    results.append(result)
            results = pd.DataFrame(results)
            results.index = results["ID"]

        purified_ct_adata_final = sc.AnnData(X=csr_matrix(purified_ct_adata_reshape))
        purified_ct_adata_final.obs = results
        # purified_ct_adata_final.obsm["spatial"] = np.array(results[["X","Y"]])
        purified_ct_adata_final.var = sp_purified_adata.var

        # purified_ct_adata_final.uns["cellperc_reconstruct"] = np.array(cellperc_pred)
        purified_ct_adata_final.obs["cellperc"] = np.array(cellperc_pred).reshape(self.obs_num * self.celltype_num,1)
                

        return(purified_ct_adata_final)

    # extended function-2
    # enhance the spot's resolution with 10 cells for each spot
    # downstream tasks: cell-cell communications
    def enhance_resolution(self,cellnum_in_spot=10):
        purified_adata = self.purify_spots()

        cellperc_pred = self.sp_adata.uns["cellperc_reconstruct"]

        enhanced_adata = None           
        prop = np.array(cellperc_pred) # 把细胞比例转换成细胞数量 
        cell_nums = np.floor(cellnum_in_spot * prop) # 可以增加超参增加单个spot组成的细胞数量

        for i, spot_cell_num in enumerate(cell_nums): #每个spot中的细胞数量
            for j, ct_cell_num in enumerate(spot_cell_num): #每个细胞的数量

                spot_index = self.obs_index[i] # 每个spot的index
                ct_index = self.celltype_order[j] # 每个ct的index

                if (ct_cell_num>0):
                    expand_index = [spot_index+"_"+ct_index for x in range(int(ct_cell_num))]
                    sub_adata = purified_adata[expand_index,:] # 使用直接复制的策略, 需要改成平均值的策略
                    
                    # sub_adata.X = sub_adata / ct_cell_num # 均分到每个细胞类型
                    if enhanced_adata is None:
                        enhanced_adata = sub_adata.copy()        
                    else:
                        enhanced_adata = ad.AnnData.concatenate(enhanced_adata, sub_adata, join='outer', batch_key='batch').copy()

        Cords = np.array(self.sp_adata.obsm["spatial"]) #spatial cords
        # 重新设置index; 重新分配spatial-grid          
        Cords_new = Utils.getHighresCords(Cords=Cords, numCell=np.sum(cell_nums,axis=1))

        enhanced_adata.obsm["spatial"] = Cords_new[:,[0,1]].astype('float64')

        return(enhanced_adata)

