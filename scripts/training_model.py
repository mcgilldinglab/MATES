import torch
import torch.nn.functional as F
from torch import nn
from torch import optim
from torch.utils.data import DataLoader, Dataset
from torch.autograd import Variable
import scipy
from scipy import sparse
import os
from os import walk
import pandas as pd
import numpy as np
from os.path import exists
import pickle
import datetime
import math
from pathlib import Path
from AutoEncoder import AutoEncoder
from MLP import MultiLayerPerceptron
from MLP import MLP_loss
import matplotlib.pyplot as plt
from os.path import join
import sys

##settings
torch.manual_seed(3407)
device = torch.device('cuda:0')
torch.cuda.empty_cache()
torch.cuda.memory_allocated()
BATCH_SIZE = 4096
AE_LR = 1e-4
EPOCHS = 200
USE_GPU = False
if USE_GPU:
    gpu_status = torch.cuda.is_available()
else:
    gpu_status = False
cur_path = os.getcwd()
    
##load data    
file_name = sys.argv[1]
bin_size = sys.argv[2]
prop = sys.argv[3]
with open('./'+file_name) as file:
    sample_list = file.readlines()
for i in range(len(sample_list)):
    sample_list[i] = sample_list[i][:-1]
for sample in sample_list:
    barcodes_file = cur_path + '/STAR_Solo/' + sample + '/'+sample+'_Solo.out'+ "/Gene/filtered/barcodes.tsv"
    if Path(barcodes_file).is_file():
        with open(barcodes_file, "r") as fh:
            barcodes = [l.rstrip() for l in fh.readlines()]
    
    path = cur_path + '/MU_Stats/'+sample
    
    if not os.path.isdir(path + '/seperate_training_'+bin_size + '_' + prop):
        os.mkdir(path + '/seperate_training_'+bin_size + '_' + prop)
        os.mkdir(path + '/seperate_training_'+bin_size + '_' + prop+'/AE_pretrain')
        os.mkdir(path + '/seperate_training_'+bin_size + '_' + prop+'/MLP')
          
    p1= path + '/Unique_TE_train_'+bin_size + '_' + prop + '.npz'
    unique_vec_matrix = scipy.sparse.load_npz(p1)
    p2 = path + '/Unique_BATCH_train_'+bin_size + '_' + prop+ '.npz'
    unique_TE_matrix = scipy.sparse.load_npz(p2)
    with open(path + '/Unique_selected_meta_'+bin_size + '_' + prop+'.pkl', 'rb') as f:
        unique_meta = pickle.load(f)
    TE_train = unique_vec_matrix.toarray()
    i = TE_train.shape[0]
    Batch_train = unique_TE_matrix.toarray().reshape(i,)
    p3 = path + '/Multi_TE_train_'+bin_size + '_' + prop+'.npz'
    p4 = path + '/Multi_Batch_train_'+bin_size + '_' + prop+'.npz'
    p5 = path + '/Multi_Region_train_'+bin_size + '_' + prop +'.npz'
    
    TE_FAM_NUMBER = len(unique_meta)
    MLP_TE_train = scipy.sparse.load_npz(p3).toarray()
    MLP_Batch_train = scipy.sparse.load_npz(p4).toarray().reshape(i,)
    MLP_Region_train = scipy.sparse.load_npz(p5).toarray()
    with open(path + '/Multi_meta_train_'+bin_size + '_' + prop+'.pkl', 'rb') as f:
        MLP_meta_train=pickle.load(f)
    
    path_dir = join(path, 'seperate_training_'+bin_size + '_' + prop)      

##AE Net Pretrain
    AENet = AutoEncoder(2001, 128, TE_FAM_NUMBER)
    ##optimizer  & loss
    optimizer = torch.optim.Adam(AENet.parameters(), lr=AE_LR)
    loss_f = nn.MSELoss()
          
    PATH = join(path_dir,'AE_pretrain/AE_state.pth')
    torch.save(AENet.state_dict(), PATH)
    train_data = []
    for i in range(len(TE_train)):
        train_data.append([TE_train[i],  Batch_train[i]])
    train_loader = DataLoader(train_data, BATCH_SIZE, shuffle = True, drop_last=True)

    TE_data = torch.Tensor(BATCH_SIZE,1*2001)
    TE_data = Variable(TE_data)
    BATCH_data = torch.Tensor(BATCH_SIZE,1*1)
    BATCH_data = Variable(BATCH_data)

    ##load to device
    if torch.cuda.is_available():
        AENet = AENet.to(device)
        TE_data = TE_data.to(device)
        BATCH_data = BATCH_data.to(device)
    loss_list = []
    AE_log = open(join(path_dir,'AE_pretrain_loss.txt'), 'w')
    for epoch in range(EPOCHS):
        loss_val = []
        starttime = datetime.datetime.now()
        AENet.train()

        if epoch+1 == EPOCHS:
            tmp_input = scipy.sparse.coo_matrix((0, 0))
            hidden_info = scipy.sparse.coo_matrix((0, 0))
            Batch_Info = np.array([])

        for step, (TE_vecs, Batch_ids) in enumerate(train_loader):

                TE_data.data.copy_(TE_vecs)

                Batch_info  = Batch_ids.clone().detach().view(BATCH_SIZE,1)
                BATCH_data.data.copy_(Batch_info)
                hidden, reconstruct = AENet(TE_data*1e6, BATCH_data, BATCH_SIZE)
                loss = loss_f(reconstruct, TE_data*1e6)
                if epoch+1 == EPOCHS:
                    Batch_Info = np.append(Batch_Info,(BATCH_data.cpu().detach().numpy().reshape(1,BATCH_SIZE)))
                    hidden_sparse = sparse.csr_matrix(hidden.cpu().detach().numpy())
                    hidden_info = scipy.sparse.vstack([hidden_info,hidden_sparse])
                    tmp_input_sparse = sparse.csr_matrix(TE_data.cpu().detach().numpy())
                    tmp_input = scipy.sparse.vstack([tmp_input,tmp_input_sparse])
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                if step % 100 == 0:
                    AENet.eval()
                loss_val.append(loss.item())

    #     compute the epoch training loss
        loss_list.append(np.mean(loss_val))
        endtime = datetime.datetime.now()
        if (epoch+1) % 10 == 0:
            p = join(path_dir, 'AE_pretrain/AE_'+str((epoch+1))+'_pretrained.pt')
            torch.save(AENet,p)
        # display the epoch training loss
        print("epoch : {}/{}, loss = {:.6f}, takes : {} seconds".format(epoch + 1, EPOCHS, np.mean(loss_val), endtime-starttime), file = AE_log)

    AE_log.close()
    

## AE get embeddings
    AENet = AutoEncoder(2001, 128, TE_FAM_NUMBER)
    AENet = torch.load(join(path_dir, 'AE_pretrain/AE_200_pretrained.pt'))
             
    MLP_meta_train_tmp = MLP_meta_train.tolist()
    MLP_train_data = []
    for i in range(len(MLP_TE_train)):
        MLP_train_data.append([MLP_TE_train[i],  MLP_Batch_train[i], MLP_Region_train[i],MLP_meta_train_tmp[i]])

    MLP_train_loader = DataLoader(MLP_train_data, BATCH_SIZE, shuffle = True, drop_last=True)

    MLP_TE_data = torch.Tensor(BATCH_SIZE,1*2001)
    MLP_TE_data = Variable(MLP_TE_data)
    MLP_BATCH_data = torch.Tensor(BATCH_SIZE,1*1)
    MLP_BATCH_data = Variable(MLP_BATCH_data)
    MLP_Region_data = torch.Tensor(BATCH_SIZE,1*5)
    MLP_Region_data = Variable(MLP_Region_data)
    MLP_meta_data = torch.Tensor(BATCH_SIZE,1*2)
    MLP_meta_data = Variable(MLP_meta_data)
    ##load to device
    if torch.cuda.is_available():
        AENet = AENet.to(device)
        MLP_TE_data = MLP_TE_data.to(device)
        MLP_BATCH_data = MLP_BATCH_data.to(device)
        MLP_Region_data = MLP_Region_data.to(device)
    #optimizer  & loss
    loss_f = nn.MSELoss()
    optimizer = torch.optim.Adam(AENet.parameters(), lr=AE_LR)
                
    EPOCHS = 1
    for epoch in range(EPOCHS):
        MLP_loss_val = []
        MLP_meta_data = [(),()]
        loss_val = []
        starttime = datetime.datetime.now()
        AENet.train()
        if (epoch == EPOCHS-1):
            problist=[]
            metalist=np.array([])
            Embedding_Info = scipy.sparse.coo_matrix((0, 0))
            Batch_Info = np.array([])
            Region_Info = []
            hidden_info = scipy.sparse.coo_matrix((0, 0))
            Meta_Data = [[],[],[]]
        for step, (TE_vecs, Batch_ids, TE_region, metainfo) in enumerate(MLP_train_loader):
                MLP_TE_data.data.copy_(TE_vecs)
                MLP_Region_data.data.copy_(TE_region)
                Batch_info  = Batch_ids.clone().detach().view(BATCH_SIZE,1)
                # Batch_info  = torch.tensor(Batch_ids).view(BATCH_SIZE,1)
                MLP_BATCH_data.data.copy_(Batch_info)
                embedding, reconstruct = AENet(MLP_TE_data*1000000, MLP_BATCH_data, BATCH_SIZE)
                loss = loss_f(reconstruct, MLP_TE_data*1000000)
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                if step % 100 == 0:
                    AENet.eval()
                if (epoch == EPOCHS-1):
                    Meta_Data[0] = Meta_Data[0]+(list(metainfo[0]))
                    Meta_Data[1] = Meta_Data[1]+(list(metainfo[1]))
                    Meta_Data[2] = Meta_Data[2]+(list(metainfo[2]))

                    Batch_Info = np.append(Batch_Info,(Batch_info.cpu().detach().numpy().reshape(1,BATCH_SIZE)))
                    Region_Info = Region_Info+(MLP_Region_data.cpu().detach().numpy().tolist())
                    hidden_sparse = sparse.csr_matrix(embedding.cpu().detach().numpy())
                    hidden_info = scipy.sparse.vstack([hidden_info,hidden_sparse])
#Prepare MLP Train Data
    MLP_meta_data_tmp = []
    for i in range(len(Meta_Data[0])):
        MLP_meta_data_tmp.append((Meta_Data[0][i], Meta_Data[1][i],Meta_Data[2][i]))
    MLP_TE_trained = hidden_info.toarray()
    MLP_Batch_trained = Batch_Info
    MLP_Region_trained = np.array(Region_Info)
    MLP_meta_trained = MLP_meta_data_tmp
    MLP_LR = 1e-6        
    MLP_trained_data = []
    for i in range(len(MLP_TE_trained)):
        MLP_trained_data.append([MLP_TE_trained[i],  MLP_Batch_trained[i], MLP_Region_trained[i],MLP_meta_trained[i]])

    MLP_trained_loader = DataLoader(MLP_trained_data, BATCH_SIZE, shuffle = True, drop_last=True)

    embeddings = torch.Tensor(BATCH_SIZE,1*128)
    embeddings = Variable(embeddings)
    MLP_BATCH_data_2 = torch.Tensor(BATCH_SIZE,1*1)
    MLP_BATCH_data_2 = Variable(MLP_BATCH_data_2)
    MLP_Region_data_2 = torch.Tensor(BATCH_SIZE,1*5)
    MLP_Region_data_2 = Variable(MLP_Region_data_2)
    MLP_meta_data_2 = torch.Tensor(BATCH_SIZE,1*2)
    MLP_meta_data_2 = Variable(MLP_meta_data_2)

##MLP
    MLP = MultiLayerPerceptron(TE_FAM_NUMBER,128)   
    PATH = join(path_dir,'MLP/MLP_state_250.pth')
    torch.save(MLP.state_dict(), PATH)
    MLP_model = torch.load(PATH)
                
    MLP_log = open(join(path_dir,'MLP_train_loss.txt'), 'w')
    ##load to device
    if torch.cuda.is_available():
        MLP = MLP.to(device)
        embeddings = embeddings.to(device)
        MLP_BATCH_data_2 = MLP_BATCH_data_2.to(device)
        MLP_Region_data_2= MLP_Region_data_2.to(device)
    #optimizer  & loss
    MLP_optimizer = torch.optim.Adam(MLP.parameters(), lr=MLP_LR)
    MLP_Loss_list = []
    criterion = MLP_loss()
    EPOCHS =200
    
    for epoch in range(EPOCHS):
        MLP_Loss_val = np.array([])
        starttime = datetime.datetime.now()
        MLP.train()
        if epoch+1 == EPOCHS:
            problist=np.array([])
            metalist=np.array([])
            Batch_Info = np.array([])
            Region_Info = []
            Meta_Data = [[],[],[]]
        for step, (embedding, Batch_ids, TE_region, metainfo) in enumerate(MLP_trained_loader):
            embeddings.data.copy_(embedding)
            Batch_info  = Batch_ids.clone().detach().view(BATCH_SIZE,1)
            MLP_BATCH_data_2.data.copy_(Batch_info)
            MLP_Region_data_2.data.copy_(TE_region)
            
            alpha = MLP(embeddings, MLP_BATCH_data_2, BATCH_SIZE)
            
            MLP_Loss = criterion(alpha, MLP_Region_data_2, BATCH_SIZE, bin_size)
            MLP_optimizer.zero_grad()
            MLP_Loss.mean().backward()
            MLP_optimizer.step()
            if epoch+1 == EPOCHS:
                Meta_Data[0] = Meta_Data[0]+(list(metainfo[0]))
                Meta_Data[1] = Meta_Data[1]+(list(metainfo[1]))
                Meta_Data[2] = Meta_Data[2]+(list(metainfo[2]))
                Batch_Info = np.append(Batch_Info,(Batch_info.cpu().detach().numpy().reshape(1,BATCH_SIZE)))
                Region_Info = Region_Info+(MLP_Region_data_2.cpu().detach().numpy().tolist())
                problist = np.append(problist, alpha.cpu().detach().numpy().reshape(BATCH_SIZE))
            MLP_Loss_val = np.append(MLP_Loss_val,MLP_Loss.cpu().detach())
        MLP_Loss_list.append(np.mean(MLP_Loss_val))
        endtime = datetime.datetime.now()
        if (epoch+1) % 10 == 0:
            p = join(path_dir,'MLP/MLP_'+str(epoch+1)+'_250.pt')
            torch.save(MLP,p)
        print("epoch : {}/{}, loss = {:.6f}, takes : {} seconds".format(epoch + 1, EPOCHS, np.mean(MLP_Loss_val), endtime-starttime), file = MLP_log)
    MLP_log.close()
    
    ##plt
    # figure, axis = plt.subplots(2)
  
    # axis[0].plot(MLP_Loss_list)
    # axis[0].set_title("MLP Loss")
    # tmp=[]
    # for i in problist:
    #     tmp.append(round(i*10)/10)
    # axis[1]=pd.Series(tmp).value_counts(sort=True).plot(kind='bar')
    # axis[1].set_title('MLP Dist')
    # plt.savefig(path_dir+'/MLP_dist_plot.png')
    print('Finised training for sample:' + sample)
