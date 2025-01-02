from pathlib import Path
import pandas as pd

import pickle
import torch
from torch import nn
from torch.utils.data import DataLoader
from torch.autograd import Variable
import scipy
from scipy import sparse
import os
import numpy as np
import datetime
from .AutoEncoder import AutoEncoder
from .MLP import MultiLayerPerceptron
from .MLP import MLP_loss
import matplotlib.pyplot as plt
from os.path import join
from tqdm import tqdm
import matplotlib.pyplot as plt

def pretrain_AE(EPOCHS, bin_size, prop, BATCH_SIZE, device, AE_LR, TE_FAM_NUMBER,
                TE_train, Batch_train, data_mode, sample=None):
    if data_mode == '10X':
        path_dir = join(os.getcwd(), 'training_'+ str(bin_size) + '_' + str(prop))    
        path_dir = join(path_dir, sample)  
    elif data_mode == 'Smart_seq':
        path_dir = join(os.getcwd(), 'training_'+ str(bin_size) + '_' + str(prop))
    AENet = AutoEncoder(2001, 128, TE_FAM_NUMBER, device)
    ##optimizer  & loss
    optimizer = torch.optim.Adam(AENet.parameters(), lr=AE_LR)
    loss_f = nn.MSELoss()
          
    PATH = join(path_dir,'AE_pretrain/AE_state.pth')
    torch.save(AENet.state_dict(), PATH)
    train_data = []
    for i in range(len(TE_train)):
        train_data.append([TE_train[i],  Batch_train[i]])
    train_loader = DataLoader(train_data, BATCH_SIZE, shuffle = True, drop_last=False)

    
    ##load to device
    if torch.cuda.is_available():
        AENet = AENet.to(device)
        
    loss_list = []
    AE_log = open(join(path_dir,'AE_pretrain_loss.txt'), 'w')
    with tqdm(total = EPOCHS) as pbar:
        for epoch in range(EPOCHS):
            loss_val = []
            starttime = datetime.datetime.now()
            AENet.train()

            if epoch+1 == EPOCHS:
                tmp_input = scipy.sparse.coo_matrix((0, 0))
                hidden_info = scipy.sparse.coo_matrix((0, 0))
                Batch_Info = np.array([])

            for step, (TE_vecs, Batch_ids) in enumerate(train_loader):
                    TE_data = torch.Tensor(len(Batch_ids),1*2001)
                    TE_data = Variable(TE_data)
                    BATCH_data = torch.Tensor(len(Batch_ids),1*1)
                    BATCH_data = Variable(BATCH_data)
                    if torch.cuda.is_available():
                        TE_data = TE_data.to(device)
                        BATCH_data = BATCH_data.to(device)
                    TE_data.data.copy_(TE_vecs)
                    Batch_info  = Batch_ids.clone().detach().view(len(Batch_ids),1)
                    BATCH_data.data.copy_(Batch_info)
                    hidden, reconstruct = AENet(TE_data*1e6, BATCH_data, len(Batch_ids), device)
                    loss = loss_f(reconstruct, TE_data*1e6)
                    if epoch+1 == EPOCHS:
                        Batch_Info = np.append(Batch_Info,(BATCH_data.cpu().detach().numpy().reshape(1,len(Batch_ids))))
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
            if (epoch+1) % 10 == 0 or epoch+1 == EPOCHS:
                p = join(path_dir, 'AE_pretrain/AE_'+str((epoch+1))+'_pretrained.pt')
                torch.save(AENet,p)
            # display the epoch training loss
            print("epoch : {}/{}, loss = {:.6f}, takes : {} seconds".format(epoch + 1, EPOCHS, np.mean(loss_val), endtime-starttime), file = AE_log)
            pbar.update(1)
    AE_log.close()

## AE get embeddings                                                                 
def get_AE_embedding(data_mode, bin_size, prop, BATCH_SIZE, device, AE_LR, TE_FAM_NUMBER,
                     MLP_TE_train, MLP_Batch_train, MLP_Region_train, MLP_meta_train, EPOCHS, sample=None):
    if data_mode == '10X':
        path_dir = join(os.getcwd(), 'training_'+ str(bin_size) + '_' + str(prop))    
        path_dir = join(path_dir, sample)  
    elif data_mode == 'Smart_seq':
        path_dir = join(os.getcwd(), 'training_'+ str(bin_size) + '_' + str(prop))
    AENet = AutoEncoder(2001, 128, TE_FAM_NUMBER, device)
    AENet = torch.load(join(path_dir, 'AE_pretrain/AE_'+str(EPOCHS)+'_pretrained.pt'))
             
    MLP_meta_train_tmp = MLP_meta_train.tolist()
    MLP_train_data = []
    for i in range(len(MLP_TE_train)):
        MLP_train_data.append([MLP_TE_train[i],  MLP_Batch_train[i], MLP_Region_train[i],MLP_meta_train_tmp[i]])

    MLP_train_loader = DataLoader(MLP_train_data, BATCH_SIZE, shuffle = True, drop_last=False)

    

    # ##load to device
    # if torch.cuda.is_available():
    #     AENet = AENet.to(device)
    #     MLP_TE_data = MLP_TE_data.to(device)
    #     MLP_BATCH_data = MLP_BATCH_data.to(device)
    #     MLP_Region_data = MLP_Region_data.to(device)

    AENet.to(device)

    #optimizer  & loss
    loss_f = nn.MSELoss()
    optimizer = torch.optim.Adam(AENet.parameters(), lr=AE_LR)

    for epoch in range(1):
        MLP_meta_data = [(),()]
        AENet.train()
        Batch_Info = np.array([])
        Region_Info = []
        hidden_info = scipy.sparse.coo_matrix((0, 0))
        Meta_Data = [[],[],[]]
        for step, (TE_vecs, Batch_ids, TE_region, metainfo) in enumerate(MLP_train_loader):
                MLP_TE_data = torch.Tensor(len(Batch_ids),1*2001)
                MLP_TE_data = Variable(MLP_TE_data)
                MLP_BATCH_data = torch.Tensor(len(Batch_ids),1*1)
                MLP_BATCH_data = Variable(MLP_BATCH_data)
                MLP_Region_data = torch.Tensor(len(Batch_ids),1*5)
                MLP_Region_data = Variable(MLP_Region_data)
                if torch.cuda.is_available():
                    MLP_TE_data = MLP_TE_data.to(device)
                    MLP_BATCH_data = MLP_BATCH_data.to(device)
                    MLP_Region_data = MLP_Region_data.to(device)
                MLP_TE_data.data.copy_(TE_vecs)
                MLP_Region_data.data.copy_(TE_region)
                Batch_info  = Batch_ids.clone().detach().view(len(Batch_ids),1)
                MLP_BATCH_data.data.copy_(Batch_info)
                embedding, reconstruct = AENet(MLP_TE_data*1000000, MLP_BATCH_data, len(Batch_ids), device)
                loss = loss_f(reconstruct, MLP_TE_data*1000000)
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                if step % 100 == 0:
                    AENet.eval()
                if (epoch == 0):
                    Meta_Data[0] = Meta_Data[0]+(list(metainfo[0]))
                    Meta_Data[1] = Meta_Data[1]+(list(metainfo[1]))
                    Meta_Data[2] = Meta_Data[2]+(list(metainfo[2]))

                    Batch_Info = np.append(Batch_Info,(Batch_info.cpu().detach().numpy().reshape(1,len(Batch_ids))))
                    Region_Info = Region_Info+(MLP_Region_data.cpu().detach().numpy().tolist())
                    hidden_sparse = sparse.csr_matrix(embedding.cpu().detach().numpy())
                    hidden_info = scipy.sparse.vstack([hidden_info,hidden_sparse])
    return Meta_Data, hidden_info, Batch_Info, Region_Info

def get_MLP_input(BATCH_SIZE, Meta_Data, hidden_info, Batch_Info, Region_Info):
    #Prepare MLP Train Data
    MLP_meta_data_tmp = []
    for i in range(len(Meta_Data[0])):
        MLP_meta_data_tmp.append((Meta_Data[0][i], Meta_Data[1][i],Meta_Data[2][i]))
    MLP_TE_trained = hidden_info.toarray()
    MLP_Batch_trained = Batch_Info
    MLP_Region_trained = np.array(Region_Info)
    MLP_meta_trained = MLP_meta_data_tmp       
    MLP_trained_data = []
    for i in range(len(MLP_TE_trained)):
        MLP_trained_data.append([MLP_TE_trained[i],  MLP_Batch_trained[i], MLP_Region_trained[i],MLP_meta_trained[i]])
#     print(len(MLP_TE_trained))
    MLP_trained_loader = DataLoader(MLP_trained_data, BATCH_SIZE, shuffle = True, drop_last=False)



    return MLP_trained_loader

def training_MLP(EPOCHS, bin_size, prop, device, BATCH_SIZE, TE_FAM_NUMBER, MLP_LR, MLP_trained_loader, data_mode, sample=None):
    if data_mode == '10X':
        path_dir = join(os.getcwd(), 'training_'+str(bin_size) + '_' + str(prop))    
        path_dir = join(path_dir, sample)  
    elif data_mode == 'Smart_seq':
        path_dir = join(os.getcwd(), 'training_'+str(bin_size) + '_' + str(prop))
    MLP = MultiLayerPerceptron(TE_FAM_NUMBER,128, device)   
    # PATH = join(path_dir,'MLP/MLP_state.pth')
    # torch.save(MLP.state_dict(), PATH)
    # MLP = torch.load(PATH)
                
    MLP_log = open(join(path_dir,'MLP_train_loss.txt'), 'w')
    ##load to device
    # if torch.cuda.is_available():
    #     MLP = MLP.to(device)
    #     embeddings = embeddings.to(device)
    #     MLP_BATCH_data_2 = MLP_BATCH_data_2.to(device)
    #     MLP_Region_data_2= MLP_Region_data_2.to(device)
    
    #optimizer  & loss
    MLP_optimizer = torch.optim.Adam(MLP.parameters(), lr=MLP_LR)
    MLP_Loss_list = []
    criterion = MLP_loss()
    with tqdm(total = EPOCHS) as pbar:
        for epoch in range(EPOCHS):
            MLP_Loss_val = np.array([])
            starttime = datetime.datetime.now()
            MLP.train()
            if epoch+1 == EPOCHS:
                problist=np.array([])
                Batch_Info = np.array([])
                Region_Info = []
                Meta_Data = [[],[],[]]
            for step, (embedding, Batch_ids, TE_region, metainfo) in enumerate(MLP_trained_loader):
                embeddings = torch.Tensor(len(Batch_ids),1*128)
                embeddings = Variable(embeddings)
                MLP_BATCH_data_2 = torch.Tensor(len(Batch_ids),1*1)
                MLP_BATCH_data_2 = Variable(MLP_BATCH_data_2)
                MLP_Region_data_2 = torch.Tensor(len(Batch_ids),1*5)
                MLP_Region_data_2 = Variable(MLP_Region_data_2)
                MLP_meta_data_2 = torch.Tensor(len(Batch_ids),1*2)
                MLP_meta_data_2 = Variable(MLP_meta_data_2)
                MLP = MLP.to(device)
                embeddings = embeddings.to(device)
                MLP_BATCH_data_2 = MLP_BATCH_data_2.to(device)
                MLP_Region_data_2= MLP_Region_data_2.to(device)
                embeddings.data.copy_(embedding)
                Batch_info  = Batch_ids.clone().detach().view(len(Batch_ids),1)
                MLP_BATCH_data_2.data.copy_(Batch_info)
                MLP_Region_data_2.data.copy_(TE_region)
                
                alpha = MLP(embeddings, MLP_BATCH_data_2, len(Batch_ids), device)

                Loss = criterion(alpha, MLP_Region_data_2, len(Batch_ids),bin_size)
                MLP_optimizer.zero_grad()
                Loss.mean().backward()
                MLP_optimizer.step()
                if epoch+1 == EPOCHS:
                    Meta_Data[0] = Meta_Data[0]+(list(metainfo[0]))
                    Meta_Data[1] = Meta_Data[1]+(list(metainfo[1]))
                    Meta_Data[2] = Meta_Data[2]+(list(metainfo[2]))
                    Batch_Info = np.append(Batch_Info,(Batch_info.cpu().detach().numpy().reshape(1,len(Batch_ids))))
                    Region_Info = Region_Info+(MLP_Region_data_2.cpu().detach().numpy().tolist())
                    problist = np.append(problist, alpha.cpu().detach().numpy().reshape(len(Batch_ids)))
                MLP_Loss_val = np.append(MLP_Loss_val,Loss.cpu().detach())
            MLP_Loss_list.append(np.mean(MLP_Loss_val))
            endtime = datetime.datetime.now()
            if (epoch+1) % 10 == 0 or epoch+1 == EPOCHS:
                p = join(path_dir,'MLP/MLP_'+str(epoch+1)+'_'+str(EPOCHS)+'.pt')
                torch.save(MLP,p)
            print("epoch : {}/{}, loss = {:.6f}, takes : {} seconds".format(epoch + 1, EPOCHS, np.mean(MLP_Loss_val), endtime-starttime), file = MLP_log)
            pbar.update(1)
    MLP_log.close()
    ##plt
    figure, axis = plt.subplots(2)
  
    axis[0].plot(MLP_Loss_list)
    axis[0].set_title("MLP Loss")
    tmp=[]
    for i in problist:
        tmp.append(round(i*10)/10)
    axis[1]=pd.Series(tmp).value_counts(sort=True).plot(kind='bar')
    axis[1].set_title('MLP Dist')
    plt.savefig(path_dir+'/MLP_dist_plot.png')
    

def MATES_train(data_mode, file_name, bin_size, prop, BATCH_SIZE= 4096, AE_LR = 1e-4, MLP_LR = 1e-6, 
                 AE_EPOCHS = 200, MLP_EPOCHS = 200, DEVICE = 'cude:0'):
    cur_path = os.getcwd()
    torch.manual_seed(3407)
    BIN_SIZE = str(bin_size)
    PROP = str(prop)

    def check_cuda_device(device='cuda:0'):
        if device == 'cpu':
            print("Running on CPU.")
        return
    
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA is not available.")

        if device not in [f'cuda:{i}' for i in range(torch.cuda.device_count())]:
            raise RuntimeError(f"CUDA device '{device}' is not available.")

        print(f"Device '{device}' is available.")

    
    try:
        check_cuda_device(DEVICE) 
    except RuntimeError as e:
        print(e)
    
    DEVICE = torch.device(DEVICE)
    print('Data Mode: ',data_mode)
    print('AE Settings:  Epoch: {:6d}, Learning Rate: {:.6f}'.format(AE_EPOCHS,AE_LR))
    print('MLP Settings: Epoch: {:6d}, Learning Rate: {:.6f}'.format(MLP_EPOCHS,MLP_LR))
    print('Batch Size: {:6d}'.format(BATCH_SIZE))
    print('Searching Bin Size: {:6d}'.format(bin_size))
    print('Dominate Proportion: {:6d}'.format(prop))
    
    if data_mode == 'Smart_seq':
        print("Loading training data...")
        path_dir = join(os.getcwd(), 'training_'+BIN_SIZE + '_' + PROP)
        if not os.path.isdir(path_dir):
            os.mkdir(path_dir)
            os.mkdir(path_dir + '/AE_pretrain')
            os.mkdir(path_dir + '/MLP')

        path = cur_path + '/MU_Stats'  
        with open(join(path, 'total_unique_TE_reads.txt'), 'r') as f:
            total_unique_TE_reads = int(f.read())
        with open(join(path, 'total_multi_TE_reads.txt'), 'r') as f:
            total_multi_TE_reads = int(f.read())
        if  total_unique_TE_reads + total_multi_TE_reads == 0:
            raise ValueError("The provided bam files don't have enough reads mapped to TE loci.")
        elif total_unique_TE_reads > 0 and total_multi_TE_reads == 0:
            #warning
            print("**Warning**: The provided bam files don't have enough multi-mapping TE reads!\n**Warning**: Skip model training, will quantify unique TE only.")
            return
        p1= path + '/Unique_TE_train_'+BIN_SIZE + '_' + PROP + '.npz'
        unique_vec_matrix = scipy.sparse.load_npz(p1)
        p2 = path + '/Unique_BATCH_train_'+BIN_SIZE + '_' + PROP+ '.npz'
        unique_TE_matrix = scipy.sparse.load_npz(p2)
        with open(path + '/Unique_selected_meta_'+BIN_SIZE + '_' + PROP+'.pkl', 'rb') as f:
            unique_meta = pickle.load(f)
        TE_train = unique_vec_matrix.toarray()
        i = TE_train.shape[0]
        try:
            Batch_train = unique_TE_matrix.toarray().reshape(i,)
        except:
            with open(join(path, 'total_unique_TE_reads.txt'), 'r') as f:
                total_unique_TE_reads = int(f.read())
            with open(join(path, 'total_multi_TE_reads.txt'), 'r') as f:
                total_multi_TE_reads = int(f.read())
            if  total_unique_TE_reads + total_multi_TE_reads == 0:
                raise ValueError("The provided bam files don't have enough reads mapped to TE loci.")
            elif total_unique_TE_reads > 0 and total_multi_TE_reads == 0:
                #warning
                print("**Warning**: The provided bam files don't have enough multi-mapping TE reads!\n**Warning**: Skip model training, will quantify unique TE only.")
                return
            else:
                raise RuntimeError("Failed to load training data.")
        p3 = path + '/Multi_TE_train_'+BIN_SIZE + '_' + PROP+'.npz'
        p4 = path + '/Multi_Batch_train_'+BIN_SIZE + '_' + PROP+'.npz'
        p5 = path + '/Multi_Region_train_'+BIN_SIZE + '_' + PROP +'.npz'
        
        TE_FAM_NUMBER = len(unique_meta)
        MLP_TE_train = scipy.sparse.load_npz(p3).toarray()
        MLP_Batch_train = scipy.sparse.load_npz(p4).toarray().reshape(i,)
        MLP_Region_train = scipy.sparse.load_npz(p5).toarray()
        with open(path + '/Multi_meta_train_'+BIN_SIZE + '_' + PROP+'.pkl', 'rb') as f:
            MLP_meta_train=pickle.load(f)

        print("Training model...")
        pretrain_AE(AE_EPOCHS, bin_size, prop, BATCH_SIZE, DEVICE, AE_LR,TE_FAM_NUMBER,
                        TE_train, Batch_train, data_mode)
            
        Meta_Data, hidden_info, Batch_Info, Region_Info = get_AE_embedding(data_mode, bin_size, prop, 
                                                                        BATCH_SIZE, DEVICE, AE_LR,TE_FAM_NUMBER,
                                                                        MLP_TE_train, MLP_Batch_train, 
                                                                        MLP_Region_train, MLP_meta_train, AE_EPOCHS)
        MLP_trained_loader = get_MLP_input(BATCH_SIZE, Meta_Data, hidden_info, Batch_Info, Region_Info)
        training_MLP(MLP_EPOCHS,bin_size, prop, DEVICE, BATCH_SIZE, TE_FAM_NUMBER, MLP_LR, MLP_trained_loader, data_mode)

        print("Finish training model.")

    elif data_mode == '10X':
        sample=file_name
        print("Loading training data for " + sample + '...')
        path_dir = join(os.getcwd(), 'training_'+BIN_SIZE + '_' + PROP)
        if not os.path.isdir(path_dir):
            os.mkdir(path_dir)
        if not os.path.isdir(join(path_dir, sample)):
            os.mkdir(join(path_dir, sample))
        if not os.path.isdir(join(path_dir, sample) + '/AE_pretrain'):
            os.mkdir(join(path_dir, sample) + '/AE_pretrain')
        if not os.path.isdir(join(path_dir, sample) + '/MLP'):
            os.mkdir(join(path_dir, sample) + '/MLP')

        path = cur_path + '/MU_Stats/'+sample   
        with open(join(path, 'total_unique_TE_reads.txt'), 'r') as f:
            total_unique_TE_reads = int(f.read())
        with open(join(path, 'total_multi_TE_reads.txt'), 'r') as f:
            total_multi_TE_reads = int(f.read())
        if  total_unique_TE_reads + total_multi_TE_reads == 0:
            raise ValueError("The provided bam files don't have enough reads mapped to TE loci.")
        elif total_unique_TE_reads > 0 and total_multi_TE_reads == 0:
            #warning
            print(f"**Warning**: The provided bam files don't have enough multi-mapping TE reads in sample: {sample}!\n**Warning**: Skip model training for sample: {sample}, will quantify unique TE only.")
            return   
        p1= path + '/Unique_TE_train_'+BIN_SIZE + '_' + PROP + '.npz'
        unique_vec_matrix = scipy.sparse.load_npz(p1)
        p2 = path + '/Unique_BATCH_train_'+BIN_SIZE + '_' + PROP+ '.npz'
        unique_TE_matrix = scipy.sparse.load_npz(p2)
        with open(path + '/Unique_selected_meta_'+BIN_SIZE + '_' + PROP+'.pkl', 'rb') as f:
            unique_meta = pickle.load(f)
        TE_train = unique_vec_matrix.toarray()
        i = TE_train.shape[0]
        try:
            Batch_train = unique_TE_matrix.toarray().reshape(i,)
        except:
            with open(join(path, 'total_unique_TE_reads.txt'), 'r') as f:
                total_unique_TE_reads = int(f.read())
            with open(join(path, 'total_multi_TE_reads.txt'), 'r') as f:
                total_multi_TE_reads = int(f.read())
            if  total_unique_TE_reads + total_multi_TE_reads == 0:
                raise ValueError("The provided bam files don't have enough reads mapped to TE loci.")
            elif total_unique_TE_reads > 0 and total_multi_TE_reads == 0:
                #warning
                print(f"**Warning**: The provided bam files don't have enough multi-mapping TE reads! in sample: {sample}.\n**Warning**: Skip model training for sample: {sample}, will quantify unique TE only.")
                return
            else:
                raise ValueError("Failed to load training data!")
        p3 = path + '/Multi_TE_train_'+BIN_SIZE + '_' + PROP+'.npz'
        p4 = path + '/Multi_Batch_train_'+BIN_SIZE + '_' + PROP+'.npz'
        p5 = path + '/Multi_Region_train_'+BIN_SIZE + '_' + PROP +'.npz'
        
        TE_FAM_NUMBER = len(unique_meta)
        MLP_TE_train = scipy.sparse.load_npz(p3).toarray()
        MLP_Batch_train = scipy.sparse.load_npz(p4).toarray().reshape(i,)
        MLP_Region_train = scipy.sparse.load_npz(p5).toarray()
        with open(path + '/Multi_meta_train_'+BIN_SIZE + '_' + PROP+'.pkl', 'rb') as f:
            MLP_meta_train=pickle.load(f)

        print("Training model for " + sample + '...')
        pretrain_AE(AE_EPOCHS, bin_size, prop, BATCH_SIZE, DEVICE, AE_LR,TE_FAM_NUMBER,
                    TE_train, Batch_train, data_mode, sample)
        
        Meta_Data, hidden_info, Batch_Info, Region_Info = get_AE_embedding(data_mode, bin_size, prop, 
                                                                        BATCH_SIZE, DEVICE, AE_LR,TE_FAM_NUMBER,
                                                                        MLP_TE_train, MLP_Batch_train, 
                                                                        MLP_Region_train, MLP_meta_train, AE_EPOCHS, sample)
        
        
        MLP_trained_loader = get_MLP_input(BATCH_SIZE, Meta_Data, hidden_info, Batch_Info, Region_Info)
        training_MLP(MLP_EPOCHS,bin_size, prop, DEVICE, BATCH_SIZE, TE_FAM_NUMBER, MLP_LR, MLP_trained_loader, data_mode, sample)

        print("Finish training model for " + sample + '.')